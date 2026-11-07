"""review: one sheet per recording for looking at its flags, the frames chosen from the tracking itself.

    python -I tools/video_qc/review.py <manifest.csv> <list.csv> <out_dir> [--draw COLUMN] [--id-columns "a b c"]
                                       [--pairs FILE] [--only REGEX]

The manifest is described in core.py; the optional columns nose and tail name those landmarks (default nose, tailbase;
body as there). list.csv names the recordings: id (or --id-columns: these columns joined with spaces make the manifest
id), flags (separated by ';'), and optionally the flags to draw, in the column --draw names (default 'draw'). When the
list has that column, a row with it empty is skipped; without it every flag is drawn. A note in brackets after a flag
is ignored ("corner (tracking changed since its decision)"). For duplicate, the column pair names the other recording
(several separated by ';') and pair_offset where its frames match (its frame = this one + offset, default 0); --pairs
reads them from a table with columns pair ("a = b") and offset_frames (b's frame t matches a's frame t + offset), as
validation/exp9_boris/scripts/31_light_checks.R writes it.

Writes <out_dir>/<id>.jpg (id with spaces and other characters as '_', as the other stills) and review_sheets.csv:
per recording the sheet, the flags drawn and the frames each part shows (tracking rows, counted from 0; video frame =
row). Without --only, sheets an earlier run listed there and this run did not write are deleted.

A sheet: the id, all flags, the flags drawn, the video and the scale; a legend; the clip (first and last frame and four
between, every tracked landmark drawn); then per flag drawn:
  length                 the first and last frames; with a clip longer than expected_s also where the protocol time
                         ends counted from the first frame and where it starts counted back from the last
  corner                 the frames in which an arena landmark (corners and outline) is farthest from its median,
                         whole and zoomed, the medians marked
  camera_shift           the 10 s in which the arena landmarks together sit farthest from their medians: a zoom on
                         each corner then, and the 10 s with the least shift for comparison
  nose_tail_order        moving frames (body faster than 10 cm/s) in which the body does not move toward the nose
                         (step . (nose - tail) <= 0), the middle of the longest such runs; arrow: body motion
  pose_spikes            the largest whole-skeleton spikes (nose, body and tail more than 1 cm off the midpoint of the
                         rows either side) with the frame before and after
  objects                the object points over the median of 15 frames, a zoom on each, and the frames in which a
                         point is farthest from its median
  duplicate              this clip and the other's at matched frames, each with its own tracking
  landmark_missing       frames from the longest gaps of the body landmark tracked least
  filled_gap             the two longest straight-line fills of the body point (a gap the tracking file bridged):
                         the last row before, the middle and the first row after each, whole frames
  arena_shape, arm_width the median outline with the length of each side
  mouse_size             moving frames with a short, a typical and a long nose to tail length, at one scale
  off_maze, oft_dropped  the frames in which the body is farthest outside the median outline (else the largest jumps)
  video:<check>:<flag>   as <flag> for the frames check; the objects, sync and clip checks draw their own stills
  any other flag         the clip only (readouts and summaries)

Landmarks: the arena (corners and outline), fixed objects (the manifest's objects, and any other landmark that stays
within 5% of the arena's size of its median in 95% of the frames, e.g. cups), nose, body and tail, and the other body
points. Dots are tracked positions in the frame shown, crosses medians over the recording. Without a scale cm limits
use --px-per-cm. Test: tests/test_video_qc_review.py.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
import traceback
import warnings

import cv2
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core  # noqa: E402
from objects import COLOURS as OBJECT_COLOURS  # noqa: E402

TILE = 320            # tile width, px
PER_ROW = 6
GAP = 2               # px between tiles
WIDTH = PER_ROW * (TILE + GAP)
FONT = cv2.FONT_HERSHEY_SIMPLEX
NOSE, BODY, TAIL, OTHER = (0, 0, 255), (0, 200, 0), (255, 90, 0), (170, 170, 170)
ARENA, MEDIAN, MOTION, MARK = (255, 255, 255), (255, 0, 255), (0, 165, 255), (255, 255, 255)
MOVING_CM_S = 10      # moving: body faster than this
STRAIGHT_PX = 0.001   # filled_gap: rows whose point steps exactly as in the rows either side (a bridged gap)
SPIKE_CM = 1.0        # pose spike: nose, body and tail this far off the midpoint of the rows either side
OFF_CM = 2.0          # corner and objects: a point this far from its median is off
SHIFT_S = 10          # camera shift: median over this many seconds
FIXED = 0.05          # a landmark within this share of the arena's size of its median in 95% of frames is fixed
STILL_CHECKS = {"objects", "sync", "clip"}   # video checks that draw their own stills


# ---- one recording --------------------------------------------------------------------------------------------------

class Recording:
    """Tracking, video and landmark roles of one manifest row; decoded frames are kept while it lives."""

    def __init__(self, rec, opts):
        self.rec, self.id = rec, rec["id"]
        self.track, self.rows = core.read_tracking(rec["tracking"])
        self.video = v = core.Video(core.field(rec, "video"), opts.ffmpeg, opts.ffprobe)
        self.fps = v.fps if np.isfinite(v.fps) and v.fps > 0 else 30.0
        self.video_frames = v.nb_frames or (int(round(v.duration * self.fps)) if np.isfinite(v.duration) else self.rows)
        self.n = min(self.rows, self.video_frames)   # rows with a frame
        scale = core.cm_per_px(rec, self.track)
        self.scaled = scale is not None
        self.cm = scale or 1 / opts.px_per_cm     # cm per video pixel
        self.corners = [k for k in core.names(rec, "corners", "tl tr br bl") if k in self.track]
        outline = core.names(rec, "outline") or core.names(rec, "corners", "tl tr br bl")
        self.outline = [k for k in outline if k in self.track]
        self.arena = list(dict.fromkeys(self.outline + self.corners))
        self.nose, self.body, self.tail = (core.field(rec, k, d) for k, d in (("nose", "nose"), ("body", "bodycentre"),
                                                                                ("tail", "tailbase")))
        self.poly = core.arena_polygon(self.track, self.outline)
        size = np.sqrt(core.polygon_area(self.poly)) if self.poly is not None else max(v.width, v.height)
        key = (self.nose, self.body, self.tail)
        self.objects = [k for k in core.names(rec, "objects") if k in self.track]
        for k in self.track:
            if k not in self.arena and k not in self.objects and k not in key:
                d = self.from_median(k)
                if np.isfinite(d).any() and np.nanpercentile(d, 95) < FIXED * size:
                    self.objects.append(k)
        self.others = [k for k in self.track if k not in self.arena and k not in self.objects and k not in key]
        self._frames = {}

    def xy(self, k):
        """(rows, 2) positions of landmark k, NaN where untracked or absent."""
        if k not in self.track:
            return np.full((self.rows, 2), np.nan)
        return np.column_stack(self.track[k]).astype(float)

    def at(self, k, i):
        if k not in self.track or not 0 <= i < self.rows:
            return np.array([np.nan, np.nan])
        return np.array([self.track[k][0][i], self.track[k][1][i]], float)

    def median(self, k):
        return core.median_point(self.track, k)

    def from_median(self, k):
        m = self.median(k)
        return np.full(self.rows, np.nan) if m is None else np.hypot(*(self.xy(k) - m).T)

    def frame(self, i):
        i = int(min(max(i, 0), max(self.n - 1, 0)))
        if i not in self._frames:
            img = self.video.frame(i)
            if img is None:
                img = np.full((self.video.height or TILE, self.video.width or TILE, 3), 40, np.uint8)
                cv2.putText(img, "no frame decoded", (10, 40), FONT, 1, (255, 255, 255), 2, cv2.LINE_AA)
            self._frames[i] = img
        return self._frames[i]

    def clock(self, i):
        m, s = divmod(round(i / self.fps, 1), 60)
        return f"{int(m)}:{s:04.1f}"


# ---- drawing --------------------------------------------------------------------------------------------------------

def ok(p):
    return p is not None and bool(np.isfinite(p).all())


class Tile:
    """An image scaled to a tile: the whole of it, or a square zoom of +-half video px around centre."""

    def __init__(self, img, centre=None, half=None, width=TILE):
        self.zoom = centre is not None and ok(np.asarray(centre, float))
        if self.zoom:
            self.z = width / (2 * half)
            self.origin = np.asarray(centre, float) - half
            warp = np.float32([[self.z, 0, -self.origin[0] * self.z], [0, self.z, -self.origin[1] * self.z]])
            self.img = cv2.warpAffine(img, warp, (width, width), flags=cv2.INTER_LINEAR, borderValue=(50, 50, 50))
        else:
            self.z, self.origin = width / img.shape[1], np.zeros(2)
            self.img = cv2.resize(img, (width, max(1, int(round(img.shape[0] * self.z)))), interpolation=cv2.INTER_AREA)

    def pt(self, p):
        q = np.clip(np.round((np.asarray(p, float) - self.origin) * self.z), -9999, 9999)
        return int(q[0]), int(q[1])


def put(img, text, at, colour=(255, 255, 255), scale=0.4):
    """Text with a dark rim, readable on any background."""
    cv2.putText(img, text, at, FONT, scale, (0, 0, 0), 3, cv2.LINE_AA)
    cv2.putText(img, text, at, FONT, scale, colour, 1, cv2.LINE_AA)


def caption(img, lines, scale=0.42):
    h = 6 + 16 * len(lines)
    cv2.rectangle(img, (0, 0), (img.shape[1], h), (0, 0, 0), -1)
    for k, line in enumerate(lines):
        cv2.putText(img, line, (5, 15 + 16 * k), FONT, scale, (255, 255, 255), 1, cv2.LINE_AA)
    return img


def draw_marks(t, r, i, highlight=(), medians=False):
    """Every tracked landmark of row i; with medians also the median outline and positions."""
    rk, ro = (5, 2) if t.zoom else (3, 1)
    if medians:
        med = [r.median(k) for k in r.outline]
        if len(med) >= 3 and all(m is not None for m in med):
            cv2.polylines(t.img, [np.array([t.pt(m) for m in med], np.int32)], True, MEDIAN, 1, cv2.LINE_AA)
        for k in r.arena + r.objects:
            m = r.median(k)
            if m is not None:
                cv2.drawMarker(t.img, t.pt(m), MEDIAN, cv2.MARKER_CROSS, 16 if t.zoom else 10, 1, cv2.LINE_AA)
    ring = [r.at(k, i) for k in r.outline]
    if len(ring) >= 3 and all(ok(p) for p in ring):
        cv2.polylines(t.img, [np.array([t.pt(p) for p in ring], np.int32)], True, ARENA, 1, cv2.LINE_AA)
    for k in r.arena:
        p = r.at(k, i)
        if ok(p):
            cv2.circle(t.img, t.pt(p), ro + 1, ARENA, -1, cv2.LINE_AA)
            if t.zoom or k in r.corners:
                x, y = t.pt(p)
                put(t.img, k, (x + 5, y - 5), ARENA, 0.4 if t.zoom else 0.33)
    for j, k in enumerate(r.objects):
        p, col = r.at(k, i), OBJECT_COLOURS[j % len(OBJECT_COLOURS)]
        if ok(p):
            x, y = t.pt(p)
            cv2.circle(t.img, (x, y), rk, col, -1, cv2.LINE_AA)
            put(t.img, k, (x + 6, y - 6), col, 0.4 if t.zoom else 0.33)
    for k in r.others:
        p = r.at(k, i)
        if ok(p):
            cv2.circle(t.img, t.pt(p), ro, OTHER, -1, cv2.LINE_AA)
    chain = [r.at(k, i) for k in (r.nose, r.body, r.tail)]
    for a, b in zip(chain, chain[1:]):
        if ok(a) and ok(b):
            cv2.line(t.img, t.pt(a), t.pt(b), OTHER, 1, cv2.LINE_AA)
    for p, col in zip(chain, (NOSE, BODY, TAIL)):
        if ok(p):
            cv2.circle(t.img, t.pt(p), rk, col, -1, cv2.LINE_AA)
    for k in highlight:
        p = r.at(k, i)
        if ok(p):
            x, y = t.pt(p)
            cv2.circle(t.img, (x, y), rk + 7, MARK, 1, cv2.LINE_AA)


def scale_bar(t, cm):
    """A bar of a round number of cm in the lower left corner of a zoom."""
    per_cm = t.z / cm
    length = next((L for L in (0.5, 1, 2, 5, 10, 20, 50) if L * per_cm >= 0.15 * t.img.shape[1]), 50)
    y, x1 = t.img.shape[0] - 10, 10 + int(round(length * per_cm))
    cv2.line(t.img, (10, y), (x1, y), (0, 0, 0), 4)
    cv2.line(t.img, (10, y), (x1, y), (255, 255, 255), 2)
    put(t.img, f"{length:g} cm", (10, y - 6))


def tile(r, i, note="", centre=None, half=None, highlight=(), medians=False, arrow=None):
    """Frame i with its tracking, whole or zoomed, captioned with frame number and time."""
    i = int(i)
    t = Tile(r.frame(i), centre, half)
    draw_marks(t, r, i, highlight, medians)
    if arrow is not None and ok(arrow[0]) and ok(arrow[1]):
        cv2.arrowedLine(t.img, t.pt(arrow[0]), t.pt(arrow[1]), MOTION, 2, cv2.LINE_AA, tipLength=0.3)
    if t.zoom:
        scale_bar(t, r.cm)
    return caption(t.img, [f"frame {i}  {r.clock(i)}"] + ([note] if note else []))


def wrap(text, width, scale):
    lines, cur = [], ""
    for word in text.split(" "):
        t = f"{cur} {word}".strip()
        if cur and cv2.getTextSize(t, FONT, scale, 1)[0][0] > width:
            lines.append(cur)
            cur = word
        else:
            cur = t
    return lines + [cur]


class Sheet:
    """Text and rows of tiles, stacked."""

    def __init__(self):
        self.blocks = []

    def text(self, text, scale=0.5, colour=(255, 255, 255)):
        lines = wrap(text, WIDTH - 20, scale)
        step = int(30 * scale) + 6
        img = np.zeros((8 + step * len(lines), WIDTH, 3), np.uint8)
        for k, line in enumerate(lines):
            cv2.putText(img, line, (10, 2 + step * (k + 1)), FONT, scale, colour, 1, cv2.LINE_AA)
        self.blocks.append(img)

    def tiles(self, tiles):
        for k in range(0, len(tiles), PER_ROW):
            row = tiles[k:k + PER_ROW]
            h = max(t.shape[0] for t in row)
            pad = cv2.BORDER_CONSTANT
            img = np.hstack([cv2.copyMakeBorder(t, 0, h - t.shape[0], 0, GAP, pad, value=0) for t in row])
            self.blocks.append(cv2.copyMakeBorder(img, 0, GAP, 0, max(0, WIDTH - img.shape[1]), pad, value=0))

    def legend(self, items):
        rows, x = [np.zeros((24, WIDTH, 3), np.uint8)], 10
        for kind, col, name in items:
            need = 24 + cv2.getTextSize(name, FONT, 0.45, 1)[0][0] + 20
            if x + need > WIDTH:
                rows.append(np.zeros((24, WIDTH, 3), np.uint8))
                x = 10
            img = rows[-1]
            if kind == "dot":
                cv2.circle(img, (x + 7, 12), 5, col, -1, cv2.LINE_AA)
            elif kind == "cross":
                cv2.drawMarker(img, (x + 7, 12), col, cv2.MARKER_CROSS, 14, 2)
            elif kind == "ring":
                cv2.circle(img, (x + 7, 12), 7, col, 1, cv2.LINE_AA)
            elif kind == "arrow":
                cv2.arrowedLine(img, (x, 12), (x + 16, 12), col, 2, cv2.LINE_AA, tipLength=0.4)
            else:
                cv2.line(img, (x, 12), (x + 16, 12), col, 2)
            cv2.putText(img, name, (x + 24, 17), FONT, 0.45, (255, 255, 255), 1, cv2.LINE_AA)
            x += need
        self.blocks += rows

    def image(self):
        return np.vstack([b[:, :WIDTH] for b in self.blocks])


# ---- choosing frames ------------------------------------------------------------------------------------------------

def pick(score, k, gap):
    """Rows of the k largest finite scores, at least gap rows apart, in time order."""
    s = np.where(np.isfinite(score), score, -np.inf)
    out = []
    for i in np.argsort(-s, kind="stable"):
        if len(out) == k or not np.isfinite(s[i]):
            break
        if all(abs(int(i) - j) >= gap for j in out):
            out.append(int(i))
    return sorted(out)


def runs(mask):
    """(start, length) of the runs of True."""
    idx = np.flatnonzero(mask)
    if not len(idx):
        return []
    return [(int(p[0]), len(p)) for p in np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1)]


def off_midpoint(m):
    """Distance of each row from the midpoint of the rows either side."""
    d = np.full(len(m), np.nan)
    d[1:-1] = np.hypot(*(m[1:-1] - (m[:-2] + m[2:]) / 2).T)
    return d


def motion(r):
    """Body step per row (px), speed (cm/s) and the moving rows."""
    b = r.xy(r.body)[:r.n]
    step = np.vstack([[np.nan, np.nan], np.diff(b, axis=0)])
    speed = np.hypot(*step.T) * r.fps * r.cm
    with np.errstate(invalid="ignore"):
        return step, speed, np.isfinite(speed) & (speed > MOVING_CM_S)


def body_length(r, rows=None):
    d = np.hypot(*(r.xy(r.nose) - r.xy(r.tail))[:r.n].T)
    d = d[rows] if rows is not None else d
    return float(np.nanmedian(d)) if np.isfinite(d).any() else 5 / r.cm


# ---- per flag -------------------------------------------------------------------------------------------------------
# Each part writes a title and tiles to the sheet and returns the rows it shows.

def part_length(r, s, ctx):
    sec, n = int(round(r.fps)), r.n
    exp = core.number(r.rec, "expected_s")
    dur = r.video_frames / r.fps
    shown = {0: "first frame", min(2 * sec, n - 1): "2 s in", max(n - 1 - 2 * sec, 0): "2 s before the last",
             n - 1: "last frame"}
    text = f"{ctx['flag']}: the clip runs {dur:.1f} s"
    if exp:
        text += f", {dur - exp:+.1f} s against the protocol's {exp:g} s"
        e = int(round(exp * r.fps))
        if e < n:
            shown[e - 1] = f"{exp:g} s after the first frame"
            shown[n - e] = f"{exp:g} s before the last frame"
    s.text(text + ". Does it hold the whole trial, from after the animal is put in to before it is taken out? With a "
           "longer clip also where the protocol time ends counted from the first frame, and starts counted back from "
           "the last.")
    rows = sorted(shown)
    s.tiles([tile(r, i, shown[i]) for i in rows])
    return rows


def part_corner(r, s, ctx):
    if not r.arena:
        s.text(f"{ctx['flag']}: no arena landmarks in the tracking")
        return []
    d = {k: r.from_median(k)[:r.n] for k in r.arena}
    d = {k: v for k, v in d.items() if np.isfinite(v).any()}
    q95 = {k: np.nanpercentile(v, 95) * r.cm for k, v in d.items()}
    off = {k: np.mean(v[np.isfinite(v)] * r.cm > OFF_CM) for k, v in d.items()}
    far = {k: np.nanmax(v) * r.cm for k, v in d.items()}
    order = sorted(d, key=lambda k: (-q95[k], -off[k], -far[k]))
    if not order:
        s.text(f"{ctx['flag']}: the arena landmarks are not tracked")
        return []
    # the landmark the flag is about (largest 95th percentile), and one frame of another that jumps far
    extra = [k for k in order[1:] if far[k] > OFF_CM][:1]
    plan = [(order[0], 3 - len(extra))] + [(k, 1) for k in extra]
    s.text(f"{ctx['flag']}: arena landmarks off their median (magenta: median outline and positions; white: tracked "
           "in the frame). 95th percentile distance from the median: " +
           ", ".join(f"{k} {q95[k]:.2f} cm" for k in order[:4]) + f"; frames more than {OFF_CM:g} cm off: " +
           ", ".join(f"{k} {off[k]:.1%}" for k in order[:4]) + "; farthest: " +
           ", ".join(f"{k} {far[k]:.1f} cm" for k in order[:4]) + ". The frames in which " +
           " and ".join(k for k, _ in plan) + (" are" if len(plan) > 1 else " is") + " farthest off: whole, and "
           "zoomed on the median (cross) and, when far from it, on the tracked point (ring).")
    rows, tiles = [], []
    for k, n in plan:
        for i in pick(np.where(d[k] * r.cm > 0.3, d[k], np.nan), n, int(round(r.fps))):
            med, p = r.median(k), r.at(k, i)
            note = f"{k} {d[k][i] * r.cm:.1f} cm from its median"
            rows.append(i)
            tiles.append(tile(r, i, note, highlight=[k], medians=True))
            if d[k][i] * r.cm < 8:
                tiles.append(tile(r, i, note, (med + p) / 2, max(0.7 * d[k][i] + 2 / r.cm, 4 / r.cm), [k], True))
            else:
                tiles += [tile(r, i, f"{k}: median", med, 4 / r.cm, [k], True),
                          tile(r, i, f"{k}: tracked", p, 4 / r.cm, [k], True)]
    if not rows:
        s.text("No arena landmark is more than 0.3 cm off its median in any frame.", 0.5)
    s.tiles(tiles)
    return rows


def part_shift(r, s, ctx):
    meds = [(k, r.median(k)) for k in r.arena]
    meds = [(k, m) for k, m in meds if m is not None]
    if not meds:
        s.text(f"{ctx['flag']}: no arena landmarks in the tracking")
        return []
    offs = np.stack([r.xy(k)[:r.n] - m for k, m in meds])            # landmarks, rows, xy
    size = max(1, int(round(SHIFT_S * r.fps)))
    bins = [(a, min(a + size, r.n)) for a in range(0, r.n, size)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        shift = np.array([np.hypot(*np.nanmedian(np.nanmedian(offs[:, a:b], axis=1), axis=0)) for a, b in bins])
    if not np.isfinite(shift).any():
        s.text(f"{ctx['flag']}: the arena landmarks are not tracked")
        return []
    hi, lo = int(np.nanargmax(shift)), int(np.nanargmin(shift))
    i, j = (bins[hi][0] + bins[hi][1]) // 2, (bins[lo][0] + bins[lo][1]) // 2
    s.text(f"{ctx['flag']}: the arena landmarks moving together. Their median offset from the median position over "
           f"{SHIFT_S} s is largest at {r.clock(bins[hi][0])}-{r.clock(bins[hi][1])}, {shift[hi] * r.cm:.2f} cm "
           f"(first {SHIFT_S} s {shift[0] * r.cm:.2f} cm, last {shift[-1] * r.cm:.2f} cm, least {shift[lo] * r.cm:.2f} "
           f"cm at {r.clock(bins[lo][0])}). The middle frame of those {SHIFT_S} s, a zoom on each corner (magenta "
           "cross: median), and the middle frame of the least shifted 10 s. When the picture moved too, the arena's "
           "edges sit beside the magenta outline.")
    half = max(3 / r.cm, 4 * shift[hi])
    tiles = [tile(r, i, f"largest shift, {shift[hi] * r.cm:.2f} cm", medians=True)]
    for k in (r.corners or r.arena)[:4]:
        tiles.append(tile(r, i, k, centre=r.median(k), half=half, medians=True))
    tiles.append(tile(r, j, f"least shift, {shift[lo] * r.cm:.2f} cm", medians=True))
    s.tiles(tiles)
    return [i, j]


def part_nose_tail(r, s, ctx):
    if r.nose not in r.track or r.tail not in r.track:
        s.text(f"{ctx['flag']}: no {r.nose} or {r.tail} in the tracking")
        return []
    nose, tail, body = r.xy(r.nose)[:r.n], r.xy(r.tail)[:r.n], r.xy(r.body)[:r.n]
    step, speed, moving = motion(r)
    ahead = np.sum(step * (nose - tail), axis=1)
    valid = moving & np.isfinite(ahead)
    back = valid & (ahead <= 0)
    longest = sorted(runs(back), key=lambda ab: -ab[1])[:PER_ROW]
    rows = sorted(a + n // 2 for a, n in longest)
    share = float(np.mean(ahead[valid] > 0)) if valid.any() else float("nan")
    s.text(f"{ctx['flag']}: the body moves toward the nose in {share:.1%} of {int(valid.sum())} moving frames (body "
           f"faster than {MOVING_CM_S} cm/s); in {int(back.sum())} it does not. The middle frames of the longest such "
           "runs; arrow: body from 3 frames before to 3 after. A swapped nose and tail base puts the red point on "
           "the tail; backing out and turning are real.")
    half = max(body_length(r, moving), 4 / r.cm)
    run_of = {a + n // 2: n for a, n in longest}
    s.tiles([tile(r, i, f"{speed[i]:.0f} cm/s, nose behind ({run_of[i]} rows)", centre=body[i], half=half,
                  arrow=(body[max(i - 3, 0)], body[min(i + 3, r.n - 1)])) for i in rows])
    return rows


def part_spikes(r, s, ctx):
    nose, body, tail = (r.xy(k)[:r.n] for k in (r.nose, r.body, r.tail))
    dn, db, dt = off_midpoint(nose), off_midpoint(body), off_midpoint(tail)
    lim = SPIKE_CM / r.cm
    with np.errstate(invalid="ignore"):
        spike = (dn > lim) & (db > lim) & (dt > lim)
    size = np.where(spike, np.minimum(np.minimum(dn, db), dt), np.nan)
    path = np.nansum(np.hypot(*np.diff(body, axis=0).T))
    smooth = body.copy()
    ix = np.flatnonzero(spike)
    smooth[ix] = (body[ix - 1] + body[ix + 1]) / 2
    added = 100 * (path - np.nansum(np.hypot(*np.diff(smooth, axis=0).T))) / path if path > 0 else 0.0
    rows = pick(size, 4, int(round(r.fps)))
    s.text(f"{ctx['flag']}: {int(spike.sum())} rows in which nose, body and tail base all sit more than {SPIKE_CM:g} "
           f"cm off the midpoint of the rows either side; they add {added:.2f}% to the body's path. The largest, each "
           "with the row before and after: a pose read from another frame jumps away from the animal and back.")
    half = max(body_length(r), 4 / r.cm)
    tiles = []
    for k, i in enumerate(rows):
        c = np.nanmean(body[max(i - 1, 0):i + 2], axis=0)
        h = max(half, 1.5 * size[i])
        tiles += [tile(r, i - 1, "row before", c, h), tile(r, i, f"spike {k + 1}: {size[i] * r.cm:.1f} cm off", c, h),
                  tile(r, i + 1, "row after", c, h)]
    s.tiles(tiles)
    return [j for i in rows for j in (i - 1, i, i + 1)]


def part_objects(r, s, ctx):
    if not r.objects:
        s.text(f"{ctx['flag']}: no fixed objects in the manifest or the tracking")
        return []
    bg = r.video.background(15, gray=False)
    near = {o: float(np.mean(d[np.isfinite(d)] * r.cm <= OFF_CM)) for o in r.objects
            for d in [r.from_median(o)[:r.n]] if np.isfinite(d).any()}
    together = []
    for j, a in enumerate(r.objects):
        for b in r.objects[j + 1:]:
            d = np.hypot(*(r.xy(a) - r.xy(b))[:r.n].T)
            if np.isfinite(d).any():
                together.append(f"{a} and {b} {np.mean(d[np.isfinite(d)] * r.cm <= OFF_CM):.1%}")
    s.text(f"{ctx['flag']}: the object points over the median of 15 frames (dots: every tracked position, cross: the "
           f"median), a zoom on each, and the frames in which a point is farthest from its median. Within {OFF_CM:g} "
           "cm of the median: " + ", ".join(f"{o} {v:.1%}" for o, v in near.items()) +
           (f"; within {OFF_CM:g} cm of each other: " + ", ".join(together) if together else "") + ".")
    tiles = []
    if bg is not None:
        for centre in [None] + [r.median(o) for o in r.objects]:
            t = Tile(bg, centre, 6 / r.cm)
            for j, o in enumerate(r.objects):
                col = OBJECT_COLOURS[j % len(OBJECT_COLOURS)]
                x, y = r.xy(o)[:r.n].T
                step = max(1, len(x) // 600)
                for p in zip(x[::step], y[::step]):
                    if ok(np.array(p)):
                        cv2.circle(t.img, t.pt(p), 2, col, -1)
                m = r.median(o)
                if m is not None:
                    cv2.drawMarker(t.img, t.pt(m), col, cv2.MARKER_CROSS, 18, 2)
                    put(t.img, o, (t.pt(m)[0] + 8, t.pt(m)[1] - 8), col)
            if t.zoom:
                scale_bar(t, r.cm)
            tiles.append(caption(t.img, ["median of 15 frames"]))
    rows = []
    for o in r.objects:
        d = r.from_median(o)[:r.n]
        for i in pick(np.where(d * r.cm > OFF_CM, d, np.nan), 2, int(round(r.fps))):
            tiles.append(tile(r, i, f"{o} {d[i] * r.cm:.1f} cm from its median", highlight=[o], medians=True))
            rows.append(i)
    s.tiles(tiles)
    return sorted(rows)


def part_pair(r, s, ctx):
    if not ctx["pairs"]:
        s.text(f"{ctx['flag']}: the other recording is not given (list column pair, or --pairs)")
        return []
    rows = []
    for other, off in ctx["pairs"]:
        orec = ctx["manifest"].get(other)
        if orec is None:
            s.text(f"{ctx['flag']}: {other} is not in the manifest")
            continue
        o = Recording(orec, ctx["opts"])
        lo, hi = max(0, -off), min(r.n, o.n - off)
        if hi - lo < 3:
            s.text(f"{ctx['flag']}: {other} does not overlap this recording at offset {off}")
            continue
        a, b = r.xy(r.body)[lo:hi], o.xy(o.body)[lo + off:hi + off]
        gap = np.hypot(*(a - b).T)
        gap = np.nanmedian(gap) * r.cm if np.isfinite(gap).any() else float("nan")
        picks = np.linspace(lo, hi - 1, 5).round().astype(int)[1:-1]
        s.text(f"{ctx['flag']}: this recording (left) and {other} (right) at matched frames (its frame = this one "
               f"{off:+d}); their body points are {gap:.2f} cm apart at the median. The same picture in both: one "
               "trial cut twice. Different pictures: the tracking belongs to one of them (it misses the animal in the "
               "other).")
        tiles = []
        for i in picks:
            tiles += [tile(r, i, r.id), tile(o, i + off, other)]
        s.tiles(tiles)
        rows += [int(i) for i in picks]
    return rows


def part_missing(r, s, ctx):
    shares = {k: float(np.isfinite(r.xy(k)[:r.n]).all(axis=1).mean()) for k in r.track
              if k not in r.arena and k not in r.objects}
    low = sorted((k for k in shares if shares[k] < 0.95), key=lambda k: (k not in (r.nose, r.body, r.tail), shares[k]))
    some = [k for k in low if shares[k] > 0]
    never = [k for k in low if shares[k] == 0]
    text = f"{ctx['flag']}: " + (", ".join(f"{k} tracked in {shares[k]:.1%} of frames" for k in some) or
                                  "every tracked landmark is in at least 95% of frames")
    if never:
        text += "; never tracked: " + " ".join(never)
    if not some:
        s.text(text + ".")
        return []
    k = some[0]
    gaps = sorted(runs(~np.isfinite(r.xy(k)[:r.n]).all(axis=1)), key=lambda ab: -ab[1])[:PER_ROW]
    rows = sorted(a + n // 2 for a, n in gaps)
    length = {a + n // 2: n for a, n in gaps}
    s.text(text + f". Frames from the longest gaps of {k}, around the body.")
    body, half = r.xy(r.body), max(1.5 * body_length(r), 5 / r.cm)
    s.tiles([tile(r, i, f"{k} untracked for {length[i] / r.fps:.1f} s", centre=body[i] if ok(body[i]) else None,
                  half=half) for i in rows])
    return rows


def part_filled(r, s, ctx):
    body = r.xy(r.body)[:r.n]
    d2 = np.full(len(body), np.inf)
    with np.errstate(invalid="ignore"):
        d2[1:-1] = np.hypot(*(body[2:] - 2 * body[1:-1] + body[:-2]).T)
    filled = d2 < STRAIGHT_PX
    fills = sorted(runs(filled), key=lambda ab: -ab[1])[:2]
    text = (f"{ctx['flag']}: {filled.sum() / r.fps:.1f} s ({filled.mean():.1%}) of the {r.body} track lie on "
            "straight-line fills: a gap the tracking file bridged, every point of it invented")
    if not fills:
        s.text(text + ".")
        return []
    s.text(text + ". The longest, each with the last row before it, its middle and the first row after it: what does "
           "the video show while the track is filled?")
    rows, tiles = [], []
    for a, n in sorted(fills):
        for i, note in ((a - 1, "last row before"), (a + n // 2, f"middle of a {n / r.fps:.1f} s fill"),
                        (a + n, "first row after")):
            i = int(min(max(i, 0), r.n - 1))
            rows.append(i)
            tiles.append(tile(r, i, note))
    s.tiles(tiles)
    return rows


def part_geometry(r, s, ctx):
    if r.poly is None:
        s.text(f"{ctx['flag']}: the arena outline is not tracked")
        return []
    i = r.n // 2
    t = Tile(r.frame(i), width=2 * TILE)
    draw_marks(t, r, i, medians=True)
    sides = []
    for a, b, ka, kb in zip(r.poly, np.roll(r.poly, -1, axis=0), r.outline, r.outline[1:] + r.outline[:1]):
        length = np.hypot(*(a - b)) * r.cm
        sides.append(f"{ka}-{kb} {length:.1f}")
        x, y = t.pt((a + b) / 2)
        put(t.img, f"{length:.1f}", (x - 14, y + 4), MEDIAN, 0.45)
    how = (f"from arena_cm {core.field(r.rec, 'arena_cm')}" if core.field(r.rec, "arena_cm") else
           "from scale_cm_per_px") if r.scaled else "no scale: assumed"
    s.text(f"{ctx['flag']}: the median outline (magenta) over the middle frame with the length of each side in cm "
           f"({how}): " + ", ".join(sides) + ".")
    s.tiles([caption(t.img, [f"frame {i}  {r.clock(i)}  median outline"])])
    return [i]


def part_size(r, s, ctx):
    if r.nose not in r.track or r.tail not in r.track:
        s.text(f"{ctx['flag']}: no {r.nose} or {r.tail} in the tracking")
        return []
    length = np.hypot(*(r.xy(r.nose) - r.xy(r.tail))[:r.n].T)
    _, _, moving = motion(r)
    use = moving & np.isfinite(length)
    if use.sum() < 10:
        s.text(f"{ctx['flag']}: the animal hardly moves (faster than {MOVING_CM_S} cm/s in {int(use.sum())} frames)")
        return []
    q = np.percentile(length[use], [5, 50, 95])
    s.text(f"{ctx['flag']}: nose to tail base while moving (body faster than {MOVING_CM_S} cm/s, {int(use.sum())} "
           f"frames): median {q[1] * r.cm:.2f} cm, 5th and 95th percentiles {q[0] * r.cm:.2f} and {q[2] * r.cm:.2f} "
           "cm. Two frames at each, at one scale; the flag compares the median with the test, sex and batch.")
    body, half = r.xy(r.body), max(1.2 * q[2], 4 / r.cm)
    rows, tiles = [], []
    for v, name in zip(q, ("short", "typical", "long")):
        for i in pick(np.where(use, -np.abs(length - v), np.nan), 2, int(round(r.fps))):
            rows.append(i)
            tiles.append(tile(r, i, f"{name}: {length[i] * r.cm:.1f} cm", centre=body[i], half=half))
    s.tiles(tiles)
    return rows


def part_outside(r, s, ctx):
    if r.poly is None:
        s.text(f"{ctx['flag']}: the arena outline is not tracked")
        return []
    body = r.xy(r.body)[:r.n]
    cont = r.poly.astype(np.float32).reshape(-1, 1, 2)
    out = np.full(r.n, np.nan)
    for i in np.flatnonzero(np.isfinite(body).all(axis=1)):
        out[i] = -cv2.pointPolygonTest(cont, (float(body[i, 0]), float(body[i, 1])), True)
    with np.errstate(invalid="ignore"):
        outside = out > 0
    half = max(2 * body_length(r), 6 / r.cm)
    if outside.any():
        rows = pick(np.where(outside, out, np.nan), 3, int(round(r.fps)))
        s.text(f"{ctx['flag']}: the body point is outside the median outline (magenta) in {int(outside.sum())} frames "
               f"({outside.sum() / r.fps:.1f} s), at most {np.nanmax(out) * r.cm:.1f} cm. The frames farthest out, "
               "whole and zoomed.")
        tiles = []
        for i in rows:
            note = f"{out[i] * r.cm:.1f} cm outside"
            tiles += [tile(r, i, note, medians=True), tile(r, i, note, centre=body[i], half=half, medians=True)]
    else:
        jump = np.r_[np.nan, np.hypot(*np.diff(body, axis=0).T)]
        rows = pick(jump, 3, int(round(r.fps)))
        s.text(f"{ctx['flag']}: the body point never leaves the median outline. The largest jumps from one row to "
               "the next, the row before and the row of the jump.")
        tiles = []
        for i in rows:
            c = np.nanmean(body[max(i - 1, 0):i + 1], axis=0)
            h = max(half, jump[i])
            tiles += [tile(r, i - 1, "row before", c, h), tile(r, i, f"jump {jump[i] * r.cm:.1f} cm", c, h)]
        rows = [j for i in rows for j in (i - 1, i)]
    s.tiles(tiles)
    return rows


PARTS = {"length": part_length, "corner": part_corner, "camera_shift": part_shift, "nose_tail_order": part_nose_tail,
         "pose_spikes": part_spikes, "objects": part_objects, "duplicate": part_pair, "landmark_missing": part_missing,
         "filled_gap": part_filled,
         "arena_shape": part_geometry, "arm_width": part_geometry, "mouse_size": part_size, "off_maze": part_outside,
         "oft_dropped": part_outside}


# ---- sheets -------------------------------------------------------------------------------------------------------

def safe(name):
    return re.sub(r"[^\w.-]+", "_", str(name))


def split_flags(text):
    """Flags separated by ';', without a note in brackets after one."""
    out = [re.sub(r"\s*\(.*\)\s*$", "", f).strip() for f in re.split(r";", text or "")]
    return list(dict.fromkeys(f for f in out if f))


def make_sheet(rec, flags, draw, pairs, manifest, opts):
    """The sheet of one recording, and the rows each of its parts shows."""
    r = Recording(rec, opts)
    s = Sheet()
    s.text(r.id, 0.8)
    s.text(f"flags: {'; '.join(flags) or 'none'}      drawn: {'; '.join(draw) or 'none'}", 0.55)
    exp = core.number(rec, "expected_s")
    s.text(f"video {os.path.basename(r.video.path)} ({r.video.width} x {r.video.height}, {r.video_frames} frames at "
           f"{r.fps:.3f}/s = {r.video_frames / r.fps:.1f} s); tracking {os.path.basename(rec['tracking'])} ({r.rows} "
           f"rows)" + (f"; expected {exp:g} s" if exp else "") + "; " +
           (f"scale {r.cm:.4f} cm/px" if r.scaled else f"no scale: {opts.px_per_cm:g} px per cm assumed"),
           0.45, (200, 200, 200))
    items = [("dot", NOSE, r.nose), ("dot", BODY, r.body), ("dot", TAIL, r.tail), ("dot", OTHER, "other body points"),
             ("line", ARENA, "arena outline (tracked in the frame)")]
    items += [("dot", OBJECT_COLOURS[j % len(OBJECT_COLOURS)], o) for j, o in enumerate(r.objects)]
    items += [("cross", MEDIAN, "median over the recording"), ("ring", MARK, "the landmark a frame is about"),
              ("arrow", MOTION, "body motion")]
    s.legend(items)
    shown = {}
    overview = sorted(set(np.linspace(0, r.n - 1, PER_ROW).round().astype(int).tolist()))
    s.text("the clip: first and last frame and four between", 0.55)
    s.tiles([tile(r, i, "first frame" if i == 0 else "last frame" if i == r.n - 1 else "") for i in overview])
    shown["clip"] = overview
    done = set()
    for f in draw:
        m = re.fullmatch(r"video:([^:]+):(.+)", f)
        if m and m.group(1) in STILL_CHECKS:
            s.text(f"{f}: the {m.group(1)} check draws its own still (stills/{m.group(1)}/{safe(r.id)}.jpg)", 0.5)
            continue
        part = PARTS.get(m.group(2) if m else f)
        if part is None:
            s.text(f"{f}: nothing to draw beyond the clip (a readout or summary)", 0.5)
            continue
        if part in done:
            continue
        done.add(part)
        shown[f] = part(r, s, {"flag": f, "pairs": pairs, "manifest": manifest, "opts": opts})
    return s.image(), shown


def read_pairs(path):
    """id -> {other id: offset}, from a table with columns pair ("a = b") and offset_frames."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            ids = [x.strip() for x in (row.get("pair") or "").split(" = ")]
            if len(ids) != 2:
                continue
            off = int(round(float(row.get("offset_frames") or 0)))
            out.setdefault(ids[0], {})[ids[1]] = -off      # a(t + off) matches b(t)
            out.setdefault(ids[1], {})[ids[0]] = off
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("list")
    ap.add_argument("out_dir")
    ap.add_argument("--draw", default="draw", help="list column naming the flags to draw (default draw; all flags "
                                                   "when the list has no such column)")
    ap.add_argument("--id-columns", help="list columns that, joined with spaces, make the manifest id")
    ap.add_argument("--pairs", help="table of duplicate pairs: pair ('a = b'), offset_frames")
    ap.add_argument("--only", help="regular expression on the id")
    ap.add_argument("--px-per-cm", type=float, default=16.0, help="for recordings without a scale")
    ap.add_argument("--ffmpeg"); ap.add_argument("--ffprobe")
    opts = ap.parse_args(argv)
    opts.ffmpeg = core.find_tool("ffmpeg", opts.ffmpeg)
    opts.ffprobe = core.find_tool("ffprobe", opts.ffprobe)

    manifest = {r["id"]: r for r in core.read_manifest(opts.manifest)}
    with open(opts.list, newline="", encoding="utf-8-sig") as fh:
        listed = list(csv.DictReader(fh))
    cols = opts.id_columns.split() if opts.id_columns else None
    pairs = read_pairs(opts.pairs) if opts.pairs else {}
    os.makedirs(opts.out_dir, exist_ok=True)
    index_path = os.path.join(opts.out_dir, "review_sheets.csv")
    before = core.read_table(index_path) if os.path.exists(index_path) else []

    todo = []
    for row in listed:
        rid = " ".join((row.get(c) or "").strip() for c in cols) if cols else (row.get("id") or "").strip()
        if not rid or (opts.only and not re.search(opts.only, rid)):
            continue
        flags = split_flags(row.get("flags"))
        draw = split_flags(row.get(opts.draw)) if opts.draw in row else flags
        if opts.draw in row and not draw:
            continue
        both = dict(pairs.get(rid, {}))
        offsets = [o.strip() for o in (row.get("pair_offset") or "").split(";")]
        for k, other in enumerate(o.strip() for o in (row.get("pair") or "").split(";")):
            if other:
                both[other] = int(float(offsets[k])) if k < len(offsets) and offsets[k] else 0
        todo.append((rid, flags, draw, list(both.items())))
    print(f"{len(todo)} recordings to draw", flush=True)

    index, t0 = [], time.time()
    for k, (rid, flags, draw, both) in enumerate(todo):
        row = {"id": rid, "sheet": "", "drawn": draw, "frames": "", "error": ""}
        if rid not in manifest:
            row["error"] = "not in the manifest"
        else:
            try:
                img, shown = make_sheet(manifest[rid], flags, draw, both, manifest, opts)
                row["sheet"] = safe(rid) + ".jpg"
                cv2.imwrite(os.path.join(opts.out_dir, row["sheet"]), img, [cv2.IMWRITE_JPEG_QUALITY, 85])
                row["frames"] = "; ".join(f"{p}: {' '.join(str(i) for i in rows)}" for p, rows in shown.items())
            except Exception as e:  # one recording's problem stops its sheet, not the run
                row["error"] = f"{type(e).__name__}: {e}"
                if os.environ.get("VIDEO_QC_DEBUG"):
                    traceback.print_exc()
        index.append(row)
        print(f"{k + 1} of {len(todo)}  {rid}: {row['sheet'] or row['error']}  {time.time() - t0:.0f} s", flush=True)

    ids = {r["id"] for r in index}
    if opts.only:
        index += [r for r in before if r["id"] not in ids]
    else:
        written = {r["sheet"] for r in index}
        for r in before:
            path = os.path.join(opts.out_dir, r.get("sheet") or "")
            if r.get("sheet") and r["sheet"] not in written and os.path.isfile(path):
                os.remove(path)
    core.write_csv(index, index_path)
    n_err = sum(1 for r in index if r["error"])
    print(f"{len(index) - n_err} sheets in {opts.out_dir}" + (f"; {n_err} not drawn (see error)" if n_err else ""))


if __name__ == "__main__":
    main()
