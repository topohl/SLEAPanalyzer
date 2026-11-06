"""sync: is the tracked animal where the video shows it, frame by frame?

A tracker that reads frames out of step writes the pose of another frame on a row: the tracked body point then sits
where the animal is a few frames earlier or later, not where the video shows it in that frame. This check finds the
animal in the video itself and compares it with the tracked body point ('body', default bodycentre) row by row.

Where: windows of the clip decoded frame after frame at reduced resolution. --sync-windows gives where they start as
fractions of the video (0 = from the first frame, 1 = ending with the last; default 0,0.5,1) or "all" for the whole
video; --sync-window-s their length. --sync-frames FILE (columns id, frame; a frame may be a range a-b) checks only the
listed rows, 2 s around each, e.g. candidates from another screen; recordings not listed are skipped. Frames are
tracking rows counted from 0.

How: the background is the median of 25 frames spread over the video, each with the tracked animal masked out. The
animal is a foreground blob inside the arena outline (plus a margin for rearing), dark or bright as the manifest's
'animal' says (auto: whichever the tracked point shows; when nothing stands out at the tracked points, whatever differs
most from the background in those frames, so that a tracking made from another video is off the animal rather than an
unreadable video). In each frame it is the blob under the tracked point, else the largest one near it, so that a
reflection in the wall or a second animal is not taken for it; the other side of this is that a tracking which itself
follows a reflection or a second animal is not seen. 'body' should be a landmark near the middle of the body, which
is what the blob's centre is compared with. A blob outside 0.4-1.6 times the animal's usual area (merged with
something, or mostly hidden) is not used. Blob centroid and body point differ by a few mm that depend on posture and
place (tail, shadow, camera angle); that offset is the median over the surrounding 3 s and is removed before comparing
(off_animal uses the recording's median instead, which a lasting offset cannot hide in).

A row is out of step when its tracked point is farther than 0.3 animal sizes from its own frame's animal but within 0.3
sizes of the animal of another frame within 8 frames, at least 0.2 sizes closer there, and the tracking jumps where such
a run of rows begins or ends (by 0.15 animal sizes more than the animal moves): a pose from another frame jumps in and
out, whereas a blob that merges with a reflection or a second animal jumps itself while the tracking moves on smoothly.
Limits are in animal sizes (the equivalent diameter of the blob, measured in each recording; about 5 cm for an adult
mouse, so 1.5 cm, 1 cm and 0.75 cm), which holds for any scale and species. Distances are reported in cm when the
recording has a scale, else in px. A pose read a few frames off shows only while the animal moves: k frames off at
speed v puts it k * v away, so slow walking can hide one.

Writes video_qc_sync_frames.csv: every out-of-step row, its distance, the frame whose animal its pose fits and the
distance there. Stills: the worst rows of a flagged recording, red cross = tracked point, green ring = where the video
shows the body point.

Flags:
  out_of_step  more than 2 rows out of step
  off_animal   the tracked point is farther than 0.6 animal sizes from the animal in more than 10% of the frames, or
               the whole tracking is that far beside it (calib_dx, calib_dy: made from another crop or scale)
  shifted      a whole window fits the video best a constant number of frames away: tracking rows and video frames do
               not line up there. In a window the decoder reached by seeking this can also be the decoder; rerun with
               --sync-windows all to decode from the first frame
  no_animal    the animal does not stand out from the floor, or is clear in fewer than half of the frames: the method
               failed on this video (contrast, lighting, outline, an animal that never moves), which says nothing about
               the tracking; no other flag is set then
"""

from __future__ import annotations

import csv
import os
import re
import warnings

import cv2
import numpy as np

import core

ARENA_PX = 360        # decode so that the arena spans about this many pixels
BACKGROUND_N = 25
LAG = 8               # frames either side a pose is compared with
CALIBRATION_S = 1.5   # the blob-to-body offset is the median over this many seconds either side
CONTEXT_S = 2         # --sync-frames: seconds checked around each listed row
OWN = 0.3             # out of step: this far from its own frame's animal, within this of another frame's (sizes) ...
BETTER = 0.2          # ... and this much closer to another frame's
JUMP = 0.15           # ... and the tracking jumps this much more than the animal where such a run begins or ends
OUT_OF_STEP_MAX = 2   # rows out of step allowed before the flag
FAR = 0.6             # off the animal: farther than this (animal sizes) ...
FAR_SHARE = 0.10      # ... in more than this share of frames
FOUND_MIN = 0.5       # below this share of frames with a clear animal the method failed
VISIBLE = 6           # the animal must differ from the floor by this many times the floor's noise
AREA = (0.4, 1.6)     # a blob is taken as the animal when its area is within this range of the animal's usual area
STILLS_N = 6

_frames_out = []      # out-of-step rows of all recordings, written by finish()


def run(rec, ctx):
    v, opts, track = ctx["video"], ctx["opts"], ctx["track"]
    body = core.field(rec, "body", "bodycentre")
    if body not in track:
        raise ValueError(f"landmark '{body}' not in the tracking")
    bx, by = track[body]
    fps = v.fps
    # containers without a frame count or stream duration (e.g. mkv): the tracking rows; decoding stops at the end
    frames = v.nb_frames or (int(v.duration * fps) if np.isfinite(v.duration * fps) else ctx["rows"])
    total = min(ctx["rows"], frames)
    scale = ctx["scale"]
    unit, conv = ("cm", scale) if scale else ("px", 1.0)

    windows = plan(rec, opts, total, fps)
    if not windows:
        return {"note": "no frames listed" if windows is None else "no listed frame within the video", "flags": []}
    geo = geometry(v, ctx["poly"])
    cal = calibrate(v, geo, bx, by, total, core.field(rec, "animal", "auto"))
    bg, sign, thr, size_px = cal["bg"], cal["sign"], cal["thr"], cal["size"]
    row = {"windows": ";".join(f"{a}-{a + n - 1}" for a, n in windows), "animal": "dark" if sign < 0 else "bright",
           "contrast": core.rnd(cal["contrast"], 1), "floor_noise": core.rnd(cal["noise"], 1)}
    if not cal["visible"]:
        row.update(frames_checked=0, note="the animal does not stand out from the floor", flags=["no_animal"])
        return row
    if cal["elsewhere"]:
        row["note"] = "the animal stands out from the floor, but not at the tracked point"

    tracked, blob, ratio, frame_no, win_id = [], [], [], [], []
    near, area = 1.5 * size_px, np.pi * (size_px / 2) ** 2
    for w, (start, n) in enumerate(windows):
        for k, img in enumerate(v.frames(start, n, width=geo["width"], crop=geo["crop"], gray=True)):
            i = start + k
            p = np.array([bx[i], by[i]])
            tracked.append(p)
            c, r = find_animal(img, bg, geo, sign, thr, p, near, area)
            blob.append(c)
            ratio.append(r)
            frame_no.append(i)
            win_id.append(w)
    P, C = np.array(tracked, float).reshape(-1, 2), np.array(blob, float).reshape(-1, 2)
    ratio, frame_no, win_id = np.array(ratio, float), np.array(frame_no), np.array(win_id)

    has_point = np.isfinite(P).all(axis=1)
    found = has_point & np.isfinite(C).all(axis=1)
    clear = found & (ratio >= AREA[0]) & (ratio <= AREA[1])
    C = np.where(clear[:, None], C, np.nan)
    row["frames_checked"] = int(has_point.sum())
    row["animal_found"] = core.rnd(found.sum() / max(has_point.sum(), 1))
    row["animal_clear"] = core.rnd(clear.sum() / max(has_point.sum(), 1))
    row[f"animal_size_{unit}"] = core.rnd(size_px * conv, 2)
    if clear.sum() < FOUND_MIN * max(has_point.sum(), 1) or clear.sum() < 10:
        row["flags"] = ["no_animal"]
        return row
    found = clear

    shift = np.nanmedian((C - P)[found], axis=0)
    row[f"calib_dx_{unit}"], row[f"calib_dy_{unit}"] = core.rnd(shift[0] * conv, 2), core.rnd(shift[1] * conv, 2)
    dist = np.hypot(*(P + shift - C).T)
    row[f"median_dist_{unit}"] = core.rnd(np.nanmedian(dist[found]) * conv, 2)
    far = found & (dist > FAR * size_px)
    row[f"far_limit_{unit}"] = core.rnd(FAR * size_px * conv, 2)
    row["share_far"] = core.rnd(far.sum() / found.sum())

    own = np.full(len(P), np.nan)
    best = np.full(len(P), np.nan)
    lag = np.zeros(len(P), int)
    lags = []
    for w in range(len(windows)):
        ix = np.where(win_id == w)[0]
        own[ix], best[ix], lag[ix] = out_of_step(P[ix], C[ix], fps)
        lags.append(window_lag(P[ix], C[ix], shift, size_px))
    cand = found & (own > OWN * size_px) & (best < OWN * size_px) & (own - best > BETTER * size_px)
    oos = np.zeros(len(P), bool)
    for w in range(len(windows)):
        ix = np.where(win_id == w)[0]
        oos[ix] = jumps_in_tracking(P[ix], C[ix], cand[ix], JUMP * size_px)
    row[f"out_of_step_limit_{unit}"] = core.rnd(OWN * size_px * conv, 2)
    row["out_of_step"] = int(oos.sum())
    row["out_of_step_longest"] = longest_run(frame_no[oos])
    row["out_of_step_offset"] = int(np.bincount(lag[oos] + LAG).argmax() - LAG) if oos.any() else ""
    row["out_of_step_frames"] = ranges(frame_no[oos])
    row["window_lag"] = ";".join("" if g is None else str(g) for g in lags)

    for j in np.where(oos)[0]:
        _frames_out.append({"id": rec["id"], "frame": int(frame_no[j]), "unit": unit,
                            "dist": core.rnd(own[j] * conv, 2), "fits_frame": int(frame_no[j] + lag[j]),
                            "dist_there": core.rnd(best[j] * conv, 2)})

    # the whole tracking sits beside the animal: made from a differently cropped or scaled video, or another video
    offset = np.hypot(*shift) > FAR * size_px
    flags = []
    if oos.sum() > OUT_OF_STEP_MAX:
        flags.append("out_of_step")
    if far.sum() > FAR_SHARE * found.sum() or offset:
        flags.append("off_animal")
    if any(g not in (None, 0) for g in lags):
        flags.append("shifted")
    row["flags"] = flags

    if opts.stills and ("out_of_step" in flags or "off_animal" in flags):
        mark = np.zeros(2) if offset else shift       # from the blob centre to where the video shows the body point
        gap = np.hypot(*(P + mark - C).T)
        if "out_of_step" in flags:
            pick = np.where(oos)[0][np.argsort(-(own - best)[oos])][:STILLS_N]
        else:
            pick = np.where(found)[0][np.argsort(-gap[found])][:STILLS_N]
        labels = []
        for j in sorted(pick, key=lambda j: frame_no[j]):
            if oos[j]:
                txt = f"{own[j] * conv:.1f} {unit} off, fits {lag[j]:+d}"
            else:
                txt = f"{gap[j] * conv:.1f} {unit} off"
            labels.append((int(frame_no[j]), P[j], C[j] - mark, txt))
        write_strip(v, labels, size_px, os.path.join(opts.stills_dir, "sync", safe(rec["id"]) + ".jpg"))
    return row


def finish(rows, opts):
    keys = ["id", "frame", "unit", "dist", "fits_frame", "dist_there"]
    path = os.path.join(opts.out_dir, "video_qc_sync_frames.csv")
    out = list(_frames_out)
    if getattr(opts, "only", None) and os.path.exists(path):   # keep the rows of recordings this run left out
        ids = {r["id"] for r in rows}
        out += [r for r in core.read_table(path) if r["id"] not in ids]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    _frames_out.clear()


# ---- where to look ------------------------------------------------------------------------------------------------

def plan(rec, opts, total, fps):
    """(start, n) windows in frames, merged where they overlap; None when --sync-frames lists nothing for rec."""
    if getattr(opts, "sync_frames", None):
        listed = listed_frames(opts).get(rec["id"])
        if not listed:
            return None
        pad = int(round(CONTEXT_S * fps))
        wins = [(max(0, a - pad), min(total, b + pad + 1)) for a, b in listed if a < total]
    elif str(opts.sync_windows).strip().lower() == "all":
        wins = [(0, total)]
    else:
        n = min(total, int(round(opts.sync_window_s * fps)))
        wins = []
        for f in str(opts.sync_windows).split(","):
            a = int(round(min(max(float(f), 0), 1) * (total - n)))
            wins.append((a, a + n))
    merged = []
    for a, b in sorted(wins):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(a, b - a) for a, b in merged if b > a]


def listed_frames(opts):
    if not hasattr(opts, "_sync_listed"):
        out = {}
        with open(opts.sync_frames, newline="", encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                m = re.fullmatch(r"\s*(\d+)\s*(?:-\s*(\d+))?\s*", r.get("frame") or "")
                if m:
                    a = int(m.group(1))
                    out.setdefault(r["id"], []).append((a, int(m.group(2) or a)))
        opts._sync_listed = out
    return opts._sync_listed


# ---- finding the animal -------------------------------------------------------------------------------------------

def geometry(v, poly):
    """Crop to the arena (with a margin), the decode width, and the arena mask in decoded pixels."""
    if poly is None:
        poly = np.array([[0, 0], [v.width, 0], [v.width, v.height], [0, v.height]], float)
    lo, hi = poly.min(axis=0), poly.max(axis=0)
    extent = float(max(hi - lo))
    margin = 0.08 * extent
    x0, y0 = (max(0, int(c - margin)) // 2 * 2 for c in lo)
    x1, y1 = min(v.width, int(hi[0] + margin)), min(v.height, int(hi[1] + margin))
    cw, ch = (x1 - x0) // 2 * 2, (y1 - y0) // 2 * 2
    width = max(2, min(cw, int(round(cw * ARENA_PX / extent))) // 2 * 2)
    height = int(round(ch * width / cw / 2) * 2)
    fx, fy = width / cw, height / ch
    inside = np.zeros((height, width), np.uint8)
    cv2.fillPoly(inside, [np.round((poly - [x0, y0]) * [fx, fy]).astype(np.int32)], 1)
    grow = max(1, int(round(0.06 * extent * fx)))
    inside = cv2.dilate(inside, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * grow + 1, 2 * grow + 1)))
    return {"crop": (x0, y0, cw, ch), "width": width, "f": np.array([fx, fy]), "origin": np.array([x0, y0]),
            "inside": inside.astype(bool), "extent": extent}


def to_decoded(geo, p):
    return (p - geo["origin"]) * geo["f"]


def to_source(geo, q):
    return q / geo["f"] + geo["origin"]


def calibrate(v, geo, bx, by, total, animal):
    """Background, polarity (-1 dark, +1 bright), foreground threshold, animal size (source px), whether the animal
    stands out from the floor at all and whether only away from the tracked point, from frames spread over the
    recording."""
    picks = np.linspace(min(30, total // 10), max(total - 30, total * 9 // 10), BACKGROUND_N).astype(int)
    imgs, pts = [], []
    for i in picks:
        img = next(v.frames(int(i), 1, width=geo["width"], crop=geo["crop"], gray=True), None)
        if img is not None:
            imgs.append(img.astype(np.float32))
            pts.append(to_decoded(geo, np.array([bx[i], by[i]])))
    if len(imgs) < 3:
        raise RuntimeError("could not decode frames for the background")
    yy, xx = np.mgrid[:imgs[0].shape[0], :imgs[0].shape[1]]
    hide = 0.12 * geo["extent"] * geo["f"][0]
    stack = np.stack(imgs)
    masked = stack.copy()
    for k, p in enumerate(pts):
        if np.isfinite(p).all():
            masked[k][(xx - p[0]) ** 2 + (yy - p[1]) ** 2 < hide ** 2] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bg = np.nanmedian(masked, axis=0)
    bg = np.where(np.isfinite(bg), bg, np.median(stack, axis=0)).astype(np.float32)

    # the animal against the floor near the tracked point, and the floor's own variation
    if not any(np.isfinite(p).all() for p in pts):
        raise RuntimeError("no tracked point in the background frames")
    diffs = [img - bg for img in imgs]
    noise = float(np.median([1.4826 * np.median(np.abs(d[geo["inside"]] - np.median(d[geo["inside"]])))
                             for d in diffs]))
    limit = max(VISIBLE * noise, 15.0)
    r = 0.4 * hide

    def contrast_at(sign, refs):
        vals = [np.percentile(sign * d[(xx - q[0]) ** 2 + (yy - q[1]) ** 2 < r ** 2], 95)
                for d, q in zip(diffs, refs) if np.isfinite(q).all()]
        return float(np.median(vals)) if vals else 0.0

    signs = [-1 if animal == "dark" else 1] if animal in ("dark", "bright") else [-1, 1]
    cands = [(contrast_at(s, pts), s, pts) for s in signs]
    if max(c[0] for c in cands) < limit:
        # nothing at the tracked point: look where each frame differs most from the background, so that a tracking
        # that does not belong to this video is reported as off the animal, not as a video the method cannot read
        k = 2 * int(r) + 1
        for s in signs:
            spots = []
            for d in diffs:
                m = cv2.blur(s * d, (k, k))
                m[~geo["inside"]] = -np.inf
                y, x = np.unravel_index(np.argmax(m), m.shape)
                spots.append(np.array([x, y], float))
            cands.append((contrast_at(s, spots), s, spots))
    best = max(range(len(cands)), key=lambda j: cands[j][0])
    contrast, sign, refs = cands[best]
    elsewhere = best >= len(signs) and contrast >= limit
    thr = max(5 * noise, 0.4 * contrast, 8.0)

    sizes = []
    for img, q in zip(imgs, refs):
        if not np.isfinite(q).all():
            continue
        areas, cents, _ = blobs(img, bg, geo, sign, thr)
        dd = np.hypot(*(cents - q).T) if len(areas) else np.array([])
        if (dd < hide).any():
            sizes.append(areas[dd < hide].max())
    size = 2 * np.sqrt(np.median(sizes) / np.pi) / geo["f"][0] if sizes else hide / geo["f"][0]
    visible = contrast >= limit and len(sizes) >= BACKGROUND_N // 3
    return {"bg": bg, "sign": sign, "thr": thr, "size": float(size), "contrast": contrast, "noise": noise,
            "visible": visible, "elsewhere": elsewhere}


def blobs(img, bg, geo, sign, thr):
    """Areas and centroids (decoded px) of the foreground blobs inside the arena, and the label image."""
    d = (img.astype(np.float32) - bg) * sign
    fg = ((d > thr) & geo["inside"]).astype(np.uint8)
    fg = cv2.morphologyEx(fg, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    n, lab, st, cen = cv2.connectedComponentsWithStats(fg, connectivity=8)
    return st[1:, cv2.CC_STAT_AREA].astype(float), cen[1:], lab


def find_animal(img, bg, geo, sign, thr, p, near, area):
    """Centroid (source px) and area (in animal areas) of the animal: the blob under the tracked point p, else the
    largest blob within `near` of it, else (p far or missing) the largest blob anywhere; NaN when none is big enough.
    Preferring the blob under the point keeps a reflection in the wall or a second animal from being taken for it."""
    areas, cents, lab = blobs(img, bg, geo, sign, thr)
    areas = areas / (area * geo["f"][0] * geo["f"][1])
    keep = areas >= 0.25
    none = np.array([np.nan, np.nan]), np.nan
    if not keep.any():
        return none
    src = to_source(geo, cents)
    if np.isfinite(p).all():
        q = np.round(to_decoded(geo, p)).astype(int)
        if 0 <= q[1] < lab.shape[0] and 0 <= q[0] < lab.shape[1]:
            j = lab[q[1], q[0]] - 1
            if j >= 0 and keep[j]:
                return src[j], areas[j]
        close = keep & (np.hypot(*(src - p).T) < near)
        if close.any():
            j = np.where(close)[0][np.argmax(areas[close])]
            return src[j], areas[j]
    j = np.where(keep)[0][np.argmax(areas[keep])]
    return src[j], areas[j]


# ---- comparing ------------------------------------------------------------------------------------------------------

def out_of_step(P, C, fps):
    """Per row: distance of the calibrated tracked point to its own frame's animal, the smallest distance to the animal
    of frames within LAG, and that frame's offset."""
    n = len(P)
    half = max(5, int(round(CALIBRATION_S * fps)))
    diff = C - P
    pad = np.full((n + 2 * half, 2), np.nan)
    pad[half:half + n] = diff
    win = np.lib.stride_tricks.sliding_window_view(pad, 2 * half + 1, axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        local = np.nanmedian(win, axis=2)
        enough = np.isfinite(win[:, 0, :]).sum(axis=1) >= 10
        local[~enough] = np.nanmedian(diff, axis=0)
    E = P + local
    own = np.hypot(*(E - C).T)
    best = np.full(n, np.inf)
    lag = np.zeros(n, int)
    for k in range(-LAG, LAG + 1):
        if k == 0 or abs(k) >= n:
            continue
        Ck = np.full((n, 2), np.nan)
        if k > 0:
            Ck[:n - k] = C[k:]
        else:
            Ck[-k:] = C[:n + k]
        dk = np.hypot(*(E - Ck).T)
        better = np.isfinite(dk) & (dk < best)
        best[better], lag[better] = dk[better], k
    best[~np.isfinite(best)] = np.nan
    return own, best, lag


def jumps_in_tracking(P, C, cand, limit):
    """Keep the runs of candidate rows where the tracking, not the video, jumps: entering or leaving the run the tracked
    point moves more than `limit` farther than the animal does. A pose read from another frame jumps in and out; a
    blob that merges with a reflection or a second animal for a while jumps itself while the tracking moves on."""
    keep = np.zeros(len(P), bool)
    idx = np.where(cand)[0]
    if len(idx) == 0:
        return keep
    for run in np.split(idx, np.where(np.diff(idx) != 1)[0] + 1):
        a, b = run[0], run[-1]
        for u, w in ((a - 1, a), (b, b + 1)):
            if u < 0 or w >= len(P) or not (np.isfinite(P[[u, w]]).all() and np.isfinite(C[[u, w]]).all()):
                continue
            if np.hypot(*(P[w] - P[u])) - np.hypot(*(C[w] - C[u])) > limit:
                keep[run] = True
                break
    return keep


def window_lag(P, C, shift, size_px):
    """The constant frame offset at which the whole window fits the video best (0 when it fits as it is), or None
    when the animal moved too little to tell."""
    n = len(C)
    if n < 2 * LAG + 30:
        return None
    moving = np.zeros(n, bool)
    moving[2:-2] = np.hypot(*(C[4:] - C[:-4]).T) > 0.1 * size_px
    if moving.sum() < 30:
        return None
    E = P + shift
    med = {}
    for k in range(-LAG, LAG + 1):
        Ck = np.full((n, 2), np.nan)
        if k >= 0:
            Ck[:n - k] = C[k:]
        else:
            Ck[-k:] = C[:n + k]
        d = np.hypot(*(E - Ck).T)[moving]
        med[k] = np.nanmedian(d) if np.isfinite(d).sum() >= 30 else np.inf
    k = min(med, key=med.get)
    return k if k != 0 and med[k] < 0.7 * med[0] else 0


def longest_run(frames):
    if len(frames) == 0:
        return 0
    runs = np.split(frames, np.where(np.diff(frames) != 1)[0] + 1)
    return int(max(len(r) for r in runs))


def ranges(frames):
    if len(frames) == 0:
        return ""
    runs = np.split(frames, np.where(np.diff(frames) != 1)[0] + 1)
    return " ".join(f"{r[0]}" if len(r) == 1 else f"{r[0]}-{r[-1]}" for r in runs)


# ---- stills ---------------------------------------------------------------------------------------------------------

def write_strip(v, labels, size_px, path):
    """One tile per frame around the animal: red cross = tracked point, green ring = where the video shows it."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    half = int(max(60, 2.0 * size_px))
    tiles = []
    for i, p, q, txt in labels:
        img = v.frame(i, gray=False)
        if img is None:
            continue
        cx, cy = (int(round(c)) for c in (q if np.isfinite(q).all() else p))
        x0, y0 = max(0, cx - half), max(0, cy - half)
        tile = img[y0:y0 + 2 * half, x0:x0 + 2 * half].copy()
        tile = cv2.copyMakeBorder(tile, 0, 2 * half - tile.shape[0], 0, 2 * half - tile.shape[1],
                                  cv2.BORDER_CONSTANT, value=(0, 0, 0))
        r = max(4, half // 12)
        if np.isfinite(q).all():
            cv2.circle(tile, (int(q[0] - x0), int(q[1] - y0)), r, (0, 220, 0), 2)
        if np.isfinite(p).all():
            cv2.drawMarker(tile, (int(p[0] - x0), int(p[1] - y0)), (0, 0, 255), cv2.MARKER_CROSS, 2 * r, 2)
        tile = cv2.resize(tile, (240, 240), interpolation=cv2.INTER_AREA)
        cv2.rectangle(tile, (0, 0), (240, 34), (0, 0, 0), -1)
        cv2.putText(tile, f"frame {i}", (4, 14), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (255, 255, 255), 1, cv2.LINE_AA)
        cv2.putText(tile, txt, (4, 29), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (255, 255, 255), 1, cv2.LINE_AA)
        tiles.append(tile)
    if tiles:
        cv2.imwrite(path, np.hstack(tiles), [cv2.IMWRITE_JPEG_QUALITY, 90])


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)
