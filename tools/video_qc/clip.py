"""clip: does the clip cover the whole trial in the raw recording it was cut from?

The clip is a spatial crop of a raw recording, at the same pixel scale unless the manifest gives clip_scale (clip px
per raw px, e.g. 0.5 for a clip that was halved). The raw is raw_video (one file, or the chapter files of one
recording in order, ';'-separated, read as one continuous video) or is searched for below raw_dir: every video file,
skipping folders named cropped, converted, old or temp, and the clip's own folder unless that is raw_dir itself, with
camera chapters joined (GoPro GX01nnnn, GX02nnnn, ... and GOPRnnnn, GP01nnnn, ...). Recordings whose file name or
number appears in the clip name are tried first (GX010786 for a clip named ..._GX010786_... or ..._GX0786_...);
recordings no larger than the clip (other cuts of the same trial, perhaps) only after all larger ones.

How:
  1. Key frames of a candidate raw are decoded at 192 px wide (the first minute first, the rest only when needed) and
     kept for the run, so the clips of one folder share them.
  2. The crop position: the clip background (median of 15 frames) matched in the median of the raw key frames.
  3. The clip's place in the raw: the first 40 s of the clip and every raw key frame, both minus the clip background
     and blurred, so mostly the animal is left, compared at every offset. Score = summed squared difference over
     summed energy, relative to what a raw showing nothing there would score (0 = the same animal in the same place
     at every key frame, about 1 = unrelated or empty). A raw with a score up to 0.4 is taken. The clip's first and
     last frames are then found to the frame by template matching at 1/4 scale (finer for clips under 1080 px) within
     +-20 frames, which also refines the crop box; when one end matches several frames alike (nothing moves), the
     pair that fits the clip's frame count is taken.
  4. The trial, from the raw inside the crop box: foreground (more than 25 grey levels from the raw's median over the
     clip, brought to the frame's brightness, so that a change of exposure does not count) that reaches from outside
     the arena (beyond the tracked outline widened by about half a body length, or across the edge of the crop) onto
     its floor (inside the outline; the whole crop when there is no outline) and is larger than 4 animals in all is
     an intrusion (a hand, a cup), unless the raw and its background correlate there (r >= 0.8: a shadow or a change
     of light, the scene shows through). The animal is the median largest foreground blob over the clip. Movement
     only outside the floor, or only on it, does not count. An intrusion seen on a single key frame only counts when
     at least 0.5 s of frames around it show it; one shorter than the key-frame spacing (1 s on GoPro) can be
     missed. The trial is the intrusion-free stretch that overlaps the clip most: the release is its start (when an
     intrusion ends there within 60 s before the clip start) and the trial ends at the next intrusion or at the end
     of the raw. An intrusion inside the clip with the animal in the arena before and after it does not end the
     trial (a hand putting something right). Times are refined to the frame.

Columns: raw (file(s)), found_by (raw_video / name / search), align_score, tried (each candidate with its score and
the seconds of it decoded, in the order tried; '-' = could not be placed), raw_start and raw_end (the clip's first
and last frame in the raw, counted over the chapters; start_at / end_at name the file and its own frame),
match_start / match_end (normalised squared difference of those matches; near 0 = the same picture),
raw_span_minus_clip_frames (clip frames dropped or added between raw and clip), clip_frames, raw_frames, raw_fps,
crop box in raw px, release (raw frame and s), start_minus_release_s, expected_end_s (release + expected_s),
end_minus_expected_s (clip end minus that), trial_end (next intrusion or raw end; trial_end_by), next_intrusion_raw_s
and next_intrusion_after_clip_s, alone_to_expected (the raw shows the animal alone up to release + expected_s),
recoverable_s (raw with the animal alone after the clip end, within the protocol time), raw_after_clip_s,
animal_cm2 (or animal_px2 without a scale: the size taken for the animal), intrusion_in_clip_s (s from the clip
start, to the key frame). Other times in s on the raw's clock.

Flags:
  raw_not_found    no candidate matched (align_score is the best try; error says why when none could be scored)
  raw_clip_sized   the raw taken is no larger than the clip although larger recordings were searched: it may be
                   another cut of the same trial, not the camera recording
  no_expected      no expected_s: positions and times only, no end flags
  release_not_seen no intrusion ends in the 60 s before the clip start: the clip start stands in for the release
  starts_early     the clip starts more than 0.5 s before the release (the hand is in it)
  starts_late      the clip starts more than --length-tolerance-s after the release
  cut_short        the clip ends more than the tolerance before the trial does (release + expected_s, or the next
                   intrusion when that comes first) while the raw shows the animal alone: the rest can be recovered
  ended_early      an intrusion or the end of the raw comes more than the tolerance before release + expected_s:
                   the trial itself was short
  ends_late        the clip runs on more than the tolerance past release + expected_s or into the next intrusion
  intrusion_in_clip  a hand (anything larger than 4 animals reaching in from outside) in the arena during the clip,
                   the animal there before and after it
"""

from __future__ import annotations

import collections
import os
import re
import subprocess
import tempfile

import cv2
import numpy as np

import core

THUMB_W = 192          # raw key-frame thumbnails, px wide (1920 -> 1/10)
HEAD_S = 60            # first pass over the candidates: the first minute of each
ALIGN_S = 40           # clip seconds matched against the raw
ACCEPT = 0.4           # alignment score up to which a raw is taken
NULL_MIN = 0.05        # the clip's first seconds must stand out this much from the noise floor to be placed at all
SCENE_MIN = 0.5        # candidates whose view correlates less with the clip background are skipped
FINE = 0.25            # scale of the exact frame search, relative to the raw ...
FINE_W = 270           # ... or larger, so that the clip is at least this many px wide there
FINE_HALF = 20         # +- frames searched around the coarse start and end
FG_DIFF = 25           # grey levels from the background: foreground
INTRUSION_ANIMALS = 4  # foreground of more than this many animals is an intrusion ...
SHADOW_R = 0.8         # ... unless raw and background correlate this well there (a shadow, a change of light)
RELEASE_S = 60         # how far before the clip start a release is looked for
AFTER_S = 60           # how far past the clip end (or the expected end) the next intrusion is looked for
EARLY_S = 0.5          # starts_early: the hand is in the clip for more than this
CONFIRM_S = 0.5        # an intrusion on one key frame counts when this much of the frames around it show it
CACHE_MB = 800         # key-frame thumbnails kept in memory over the run
VIDEO_EXT = {".mp4", ".mov", ".avi", ".mkv", ".m4v", ".mts", ".mpg", ".mpeg", ".wmv"}
SKIP_DIRS = {"cropped", "converted", "old", "temp"}

_RECS = collections.OrderedDict()   # paths -> Recording (or None when unreadable), least recently used first
_DIRS = {}                          # raw_dir -> recordings found below it


# ---- raw recordings ---------------------------------------------------------------------------------------------------

class Recording:
    """One raw recording: its chapter files in order, read as one continuous video, and its key-frame thumbnails."""

    def __init__(self, paths, opts):
        self.paths = list(paths)
        self.videos = [core.Video(p, opts.ffmpeg, opts.ffprobe) for p in self.paths]
        v = self.videos[0]
        self.fps, self.width, self.height = v.fps, v.width, v.height
        self.counts = [x.nb_frames or int(round(x.duration * x.fps)) for x in self.videos]
        self.starts = [int(s) for s in np.cumsum([0] + self.counts[:-1])]
        self.n = int(sum(self.counts))
        self.tw = min(THUMB_W, self.width)
        self.th = int(round(self.height * self.tw / self.width / 2) * 2)
        self.fx, self.fy = self.tw / self.width, self.th / self.height
        self.keys = np.zeros(0, int)
        self.thumbs = np.zeros((0, self.th, self.tw), np.uint8)
        self.upto = 0                      # key frames are known for raw frames [0, upto)
        self._scene = None

    def name(self):
        return "+".join(os.path.basename(p) for p in self.paths)

    def at(self, frame):
        """'file:frame' of a frame counted over the chapters."""
        c = max(i for i, s in enumerate(self.starts) if s <= frame)
        return f"{os.path.basename(self.paths[c])}:{frame - self.starts[c]}"

    def _keys(self, c, a, b):
        """Key frames of chapter c with own frame numbers in [a, b): numbers and thumbnails."""
        v = self.videos[c]
        # -skip_frame, not -discard: dropping the other packets before the decoder shifts the key frames' times by a
        # frame or two on video with B-frames
        cmd = [v.ffmpeg, "-hide_banner", "-nostats", "-loglevel", "info", "-threads", "2", "-skip_frame", "nokey"]
        if a > 0:
            cmd += ["-ss", f"{a / v.fps:.6f}"]
        cmd += ["-t", f"{(b - a) / v.fps:.6f}", "-i", v.path, "-map", "0:v:0", "-vsync", "0",
                "-vf", f"showinfo,scale={self.tw}:{self.th}:flags=area,format=gray",
                "-f", "rawvideo", "-pix_fmt", "gray", "-"]
        size = self.tw * self.th
        frames = []
        with tempfile.TemporaryFile() as err:     # showinfo writes each frame's time to stderr
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err)
            while True:
                buf = p.stdout.read(size)
                if len(buf) < size:
                    break
                frames.append(np.frombuffer(buf, np.uint8).reshape(self.th, self.tw))
            p.stdout.close()
            p.wait()
            err.seek(0)
            times = [float(t) for t in re.findall(rb"pts_time:\s*(-?[\d.]+)", err.read())]
        n = min(len(frames), len(times))
        idx = a + np.round(np.array(times[:n]) * v.fps).astype(int)
        keep = (idx >= a) & (idx < b)
        return idx[keep], [f for f, k in zip(frames[:n], keep) if k]

    def extend(self, upto):
        """Decode the key frames up to raw frame `upto` (not before)."""
        upto = min(int(upto), self.n)
        if upto <= self.upto:
            return
        print(f"  clip: key frames of {self.name()} to {upto / self.fps:.0f} s", flush=True)
        keys, thumbs = list(self.keys), list(self.thumbs)
        gap = max(1, int(round(0.4 * self.fps)))    # all-intra video: keep about 2.5 key frames per second
        for c, (s, cnt) in enumerate(zip(self.starts, self.counts)):
            a, b = max(self.upto, s), min(upto, s + cnt)
            if a >= b:
                continue
            idx, frames = self._keys(c, a - s, b - s)
            for i, f in zip(idx + s, frames):
                if not keys or i - keys[-1] >= gap:
                    keys.append(int(i))
                    thumbs.append(f)
        self.keys = np.array(keys, int)
        self.thumbs = np.stack(thumbs) if thumbs else np.zeros((0, self.th, self.tw), np.uint8)
        self.upto = upto

    def frames(self, start, n, **kw):
        """(raw frame, image) for frames start .. start+n-1, across chapters (core.Video.frames arguments)."""
        start, end = max(0, int(start)), min(self.n, int(start) + int(n))
        for c, (s, cnt) in enumerate(zip(self.starts, self.counts)):
            a, b = max(start, s), min(end, s + cnt)
            if a < b:
                for i, f in enumerate(self.videos[c].frames(a - s, b - a, **kw)):
                    yield a + i, f

    def background(self):
        """Median of up to 25 key frames spread over what is decoded."""
        pick = np.unique(np.linspace(0, len(self.keys) - 1, min(25, len(self.keys))).astype(int))
        return np.median(self.thumbs[pick], axis=0).astype(np.uint8)

    def scene(self):
        """One thumbnail from the middle of the first file, for the quick 'same view?' test."""
        if self._scene is None:
            self._scene = self.videos[0].frame(self.counts[0] // 2, width=self.tw, gray=True)
        return self._scene


def recording(paths, opts):
    """The run's Recording for these files (None when unreadable); the least recently used ones are dropped."""
    key = tuple(os.path.normcase(os.path.normpath(p)) for p in paths)     # the same file however it is written
    if key in _RECS:
        _RECS.move_to_end(key)
        return _RECS[key]
    try:
        r = Recording(paths, opts)
    except Exception as e:  # an unreadable file in a searched folder is skipped
        print(f"  clip: cannot read {'+'.join(map(os.path.basename, paths))}: {e}", flush=True)
        r = None
    _RECS[key] = r
    trim_cache()
    return r


def trim_cache():
    total = sum(x.thumbs.nbytes for x in _RECS.values() if x is not None)
    while total > CACHE_MB * 1e6 and len(_RECS) > 1:
        _, old = _RECS.popitem(last=False)
        total -= old.thumbs.nbytes if old is not None else 0


def recordings_below(raw_dir):
    """Recordings below raw_dir, each a list of files in chapter order."""
    if raw_dir in _DIRS:
        return _DIRS[raw_dir]
    found = []
    for root, dirs, files in os.walk(raw_dir):
        dirs[:] = sorted(d for d in dirs if d.lower() not in SKIP_DIRS)
        groups = {}
        for f in sorted(files):
            stem, ext = os.path.splitext(f)
            if ext.lower() not in VIDEO_EXT:
                continue
            m = re.match(r"^G([HX])(\d\d)(\d{4})$", stem, re.I)
            m_old = re.match(r"^GOPR(\d{4})$", stem, re.I) or re.match(r"^GP(\d\d)(\d{4})$", stem, re.I)
            if m:
                key, chapter = ("G" + m.group(1).upper(), m.group(3)), int(m.group(2))
            elif m_old:
                key = ("GOPR", m_old.groups()[-1])
                chapter = 0 if len(m_old.groups()) == 1 else int(m_old.group(1))
            else:
                key, chapter = (stem,), 0
            groups.setdefault(key, []).append((chapter, os.path.join(root, f).replace(os.sep, "/")))
        for key in sorted(groups):
            found.append([p for _, p in sorted(groups[key])])
    _DIRS[raw_dir] = found
    return found


def name_hint(clip_path, paths):
    """2 when a raw file name is in the clip name, 1 when its number is (GX010786 ~ 'GX0786'), else 0."""
    clip = os.path.splitext(os.path.basename(clip_path))[0].lower()
    groups = re.findall(r"\d{3,}", clip)
    best = 0
    for p in paths:
        stem = os.path.splitext(os.path.basename(p))[0].lower()
        if stem in clip:
            return 2
        m = re.search(r"(\d{3,})$", stem)
        if m and any(g.endswith(m.group(1)[-4:]) for g in groups):
            best = 1
    return best


# ---- the clip ---------------------------------------------------------------------------------------------------------

class Clip:
    def __init__(self, rec, video):
        self.v, self.path = video, rec["video"]
        self.scale = core.number(rec, "clip_scale") or 1.0       # clip px per raw px
        self.n = video.nb_frames or video.count_frames()
        self.fps, self.w, self.h = video.fps, video.width, video.height
        raw_w = self.w / self.scale
        self.fine_w = int(round(raw_w * min(1.0, max(FINE, FINE_W / raw_w))))
        self.fine_s = self.fine_w / raw_w                           # fine px per raw px
        self.bg = video.background(15, width=self.fine_w, gray=True)
        self._seg = {}

    def thumb_size(self, fx, fy):
        return int(round(self.w * fx / self.scale)), int(round(self.h * fy / self.scale))

    def bg_thumb(self, fx, fy):
        return cv2.resize(self.bg, self.thumb_size(fx, fy), interpolation=cv2.INTER_AREA)

    def segment(self, fx, fy, n):
        """The first n frames at the thumbnail scale of a raw."""
        key = (round(fx, 6), round(fy, 6))
        if key not in self._seg or len(self._seg[key]) < n:
            size = self.thumb_size(fx, fy)
            self._seg[key] = np.stack([cv2.resize(f, size, interpolation=cv2.INTER_AREA)
                                       for f in self.v.frames(0, n, width=self.fine_w, gray=True)])
        return self._seg[key][:n]

    def fine(self, i, gray=True, width=None):
        for f in self.v.frames(max(0, i), 1, width=width or self.fine_w, gray=gray):
            return f
        return None


# ---- finding the clip in a raw ----------------------------------------------------------------------------------------

def place(tmpl, img):
    """Sub-pixel (x, y) of tmpl in img and the peak normalised correlation."""
    tmpl = tmpl[:img.shape[0], :img.shape[1]].astype(np.float32)
    res = cv2.matchTemplate(img.astype(np.float32), tmpl, cv2.TM_CCOEFF_NORMED)
    _, peak, _, (x, y) = cv2.minMaxLoc(res)

    def sub(a, i):  # vertex of the parabola through the peak and its neighbours
        if 0 < i < len(a) - 1:
            den = a[i - 1] - 2 * a[i] + a[i + 1]
            if den < 0:
                return i + 0.5 * (a[i - 1] - a[i + 1]) / den
        return float(i)
    return sub(res[y, :], x), sub(res[:, x], y), float(peak)


def _feature(imgs, ox, oy, size):
    """Thumbnails shifted by the sub-pixel offset onto the raw grid, halved and blurred (float32, one row each)."""
    wt, ht = size
    m = np.float32([[1, 0, -ox], [0, 1, -oy]])
    out = []
    for im in imgs:
        if ox or oy:
            im = cv2.warpAffine(im, m, (wt, ht), flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_REPLICATE)
        else:
            im = im[:ht, :wt]
        im = cv2.resize(im.astype(np.float32), (wt // 2, ht // 2), interpolation=cv2.INTER_AREA)
        out.append(cv2.GaussianBlur(im, (0, 0), 1.0).ravel())
    return np.array(out, np.float32)


def align(r, clip, rate):
    """Best place of the clip's first seconds in what is decoded of r: dict(d, score, ...) or None."""
    if len(r.keys) < 3:
        return None
    x, y, peak = place(clip.bg_thumb(r.fx, r.fy), r.background())
    if peak < SCENE_MIN:
        return {"score": float("nan"), "peak": peak, "r": r}
    w0, h0 = clip.thumb_size(r.fx, r.fy)
    xa, ya = int(np.ceil(x - 1e-6)), int(np.ceil(y - 1e-6))
    ox, oy = xa - x, ya - y
    size = (min(w0 - 1, r.tw - xa) // 2 * 2, min(h0 - 1, r.th - ya) // 2 * 2)
    gap = int(np.median(np.diff(r.keys))) if len(r.keys) > 1 else 1
    step = next((s for s in (5, 4, 6, 3, 2) if gap % s == 0), 1)
    seg = clip.segment(r.fx, r.fy, min(clip.n, int(round(max(ALIGN_S, 12 * gap / r.fps) * clip.fps))))
    j = np.arange(0, len(seg), step)
    jr = np.round(j * rate).astype(int)                     # raw frames after the clip start
    seg = seg[j]
    bg = _feature([clip.bg_thumb(r.fx, r.fy)], ox, oy, size)[0]
    cd = _feature(seg, ox, oy, size) - bg
    rd = _feature(r.thumbs[:, ya:ya + size[1], xa:xa + size[0]], 0, 0, size) - bg
    cd -= np.median(cd, axis=1, keepdims=True)                # brightness drift
    rd -= np.median(rd, axis=1, keepdims=True)
    ec, er = (cd ** 2).sum(1), (rd ** 2).sum(1)
    cross = rd @ cd.T
    off = int(jr[-1])
    num = np.zeros(r.upto + off + 1)
    den = np.zeros_like(num)
    sig = np.zeros_like(num)                                # the clip's own signal (the animal) at each offset
    cnt = np.zeros(len(num), int)
    for k, rk in enumerate(r.keys):
        idx = rk - jr + off
        num[idx] += er[k] + ec - 2 * cross[k]
        den[idx] += er[k] + ec
        sig[idx] += ec
        cnt[idx] += 1
    d = np.arange(len(num)) - off
    expect = np.searchsorted(r.keys, d + off, "right") - np.searchsorted(r.keys, d, "left")
    last = int(round((clip.n - 1) * rate))
    ok = (d >= 0) & (d + off < r.upto) & (d + last < r.n + 2 * r.fps) & (cnt >= np.maximum(4, 0.5 * expect))
    if not ok.any():
        return {"score": float("nan"), "peak": peak, "r": r}
    eps = cd.shape[1] * 4.0                                 # noise floor: 2 grey levels on each side
    score = np.where(ok, num / (den + cnt * eps + 1e-9), np.inf)
    i = int(np.argmin(score))
    # A raw showing nothing there (an empty arena, a still scene) scores sig / (sig + floor), which is small for a faint
    # or small animal; the score is taken relative to that, so such a raw scores 1.
    null = sig[i] / (sig[i] + cnt[i] * eps)
    if null < NULL_MIN:
        return {"score": float("nan"), "peak": peak, "r": r, "still": True}
    return {"r": r, "d": int(d[i]), "score": float(score[i] / null), "peak": peak,
            "box": (x / r.fx, y / r.fy, clip.w / clip.scale, clip.h / clip.scale)}


def find(cands, clip, opts):
    """The first candidate recording the clip aligns with: (match, best try, 'name score' per try). Candidates no
    larger than the clip (another cut of the same trial, perhaps) come after all larger ones."""
    recs, log = [], []
    for paths in cands:
        r = recording(paths, opts)
        if r is None:
            log.append(f"{'+'.join(map(os.path.basename, paths))} unreadable")
            continue
        rt = r.fps / clip.fps if abs(r.fps / clip.fps - 1) > 0.01 else 1.0
        if r.n < (clip.n - 1) * rt + 1 - 2 * r.fps:            # too short to hold the clip
            continue
        small = r.width <= clip.w / clip.scale + 2 and r.height <= clip.h / clip.scale + 2
        recs.append((r, rt, name_hint(clip.path, paths), small))
    if len(recs) > 1:   # quick view test, then: larger than the clip, named ones, ones decoded already, the most alike
        scored = []
        for r, rt, h, small in recs:
            sc = r.scene()
            peak = place(clip.bg_thumb(r.fx, r.fy), sc)[2] if sc is not None else -1
            if peak >= SCENE_MIN or h:
                scored.append(((small, -h, -min(r.upto, 1), -peak), (r, rt, h, small)))
        recs = [c for _, c in sorted(scored, key=lambda t: t[0])]
    best, tried = None, set()
    larger = any(not small for *_, small in recs)
    for tier in (False, True):
        for whole in (False, True):
            for r, rate, _, small in recs:
                if small != tier or (whole and r.upto >= r.n and id(r) in tried):
                    continue
                r.extend(r.n if whole else max(r.upto, HEAD_S * r.fps))
                trim_cache()
                m = align(r, clip, rate)
                tried.add(id(r))
                ok = m is not None and np.isfinite(m["score"])
                log.append(f"{r.name()} {core.rnd(m['score']) if ok else '-'} ({r.upto / r.fps:.0f} s)")
                if m is None:
                    continue
                m.update(rate=rate, clip_sized=small and larger)
                if best is None or (np.isfinite(m["score"]) and not best["score"] <= m["score"]):
                    best = m
                if m["score"] <= ACCEPT:
                    return m, best, log
    return None, best, log


def locate(r, tmpl, box, lo, n, fine_s, extend=True):
    """How well tmpl (a clip frame at the fine scale) matches the crop region of raw frames lo .. lo+n-1:
    {raw frame: (normalised squared difference, refined box)}. When the best frame is on the edge of the window,
    the window next to it is searched as well."""
    x0, y0, w, h = box
    X, Y, W, H = even_crop(r, x0 - 24, y0 - 24, w + 48, h + 48)    # 24 raw px of play around the crop box
    t = tmpl[2:-2, 2:-2].astype(np.float32)
    out = {}
    lo = max(0, int(lo))
    for _ in range(2 if extend else 1):
        for i, f in r.frames(lo, n, width=int(round(W * fine_s)), crop=(X, Y, W, H), gray=True):
            if i in out or f.shape[0] < t.shape[0] or f.shape[1] < t.shape[1]:
                continue
            res = cv2.matchTemplate(f.astype(np.float32), t, cv2.TM_SQDIFF_NORMED)
            mn, _, loc, _ = cv2.minMaxLoc(res)
            sx, sy = f.shape[1] / W, f.shape[0] / H
            out[i] = (mn, (X + loc[0] / sx - 2 / fine_s, Y + loc[1] / sy - 2 / fine_s, w, h))
        if not out:
            break
        best = min(out, key=lambda i: out[i][0])
        if best == min(out) and best > 0:
            lo = max(0, best - n + 1)
        elif best == max(out) and best < r.n - 1:
            lo = best
        else:
            break
    return out


def even_crop(r, x, y, w, h):
    """A crop inside the raw frame with even position and size: ffmpeg rounds odd ones down on 4:2:0 video."""
    x, y = max(0, int(x)) // 2 * 2, max(0, int(y)) // 2 * 2
    w = min(int(np.ceil(w)), r.width - x) // 2 * 2
    h = min(int(np.ceil(h)), r.height - y) // 2 * 2
    return x, y, w, h


def place_clip(r, clip, m, rate):
    """The clip's first and last frame in the raw, to the frame: (start, end, start score, end score, box).
    Each end is matched on its own; when the two disagree with the clip's frame count and one end is ambiguous (no
    motion: several frames match alike), the pair that fits the count is taken if it matches nearly as well."""
    first, last = clip.fine(0), clip.fine(clip.n - 1)
    span = int(round((clip.n - 1) * rate))
    S = locate(r, first, m["box"], m["d"] - FINE_HALF, 2 * FINE_HALF + 1, clip.fine_s)
    if not S:
        return m["d"], m["d"] + span, float("nan"), float("nan"), m["box"]
    s = min(S, key=lambda i: S[i][0])
    box = S[s][1]
    half = FINE_HALF + int(0.003 * clip.n)
    E = locate(r, last, box, s + span - half, 2 * half + 1, clip.fine_s)
    if not E:
        return s, min(s + span, r.n - 1), S[s][0], float("nan"), box
    e = min(E, key=lambda i: E[i][0])
    if e - s != span:
        if e - span not in S and e - span >= 0:
            S.update(locate(r, first, box, e - span, 1, clip.fine_s, extend=False))
        if s + span not in E and s + span < r.n:
            E.update(locate(r, last, box, s + span, 1, clip.fine_s, extend=False))
        pairs = [(S[i][0] + E[i + span][0], i) for i in S if i + span in E]
        if pairs:
            total, i = min(pairs)
            if total <= 1.5 * (S[s][0] + E[e][0]) + 0.002:
                s, e = i, i + span
    return s, e, S[s][0], E[e][0], S[s][1]


# ---- the trial around the clip ----------------------------------------------------------------------------------------

class Scene:
    """Foreground in the crop box of the raw at thumbnail scale: animal size, intrusions."""

    def __init__(self, r, box, s, e, poly, clip_scale):
        self.r = r
        x0, y0, w, h = box
        self.xa, self.ya = max(0, int(round(x0 * r.fx))), max(0, int(round(y0 * r.fy)))
        self.wt, self.ht = min(int(w * r.fx), r.tw - self.xa), min(int(h * r.fy), r.th - self.ya)
        crop = r.thumbs[:, self.ya:self.ya + self.ht, self.xa:self.xa + self.wt]
        inside = np.flatnonzero((r.keys >= s) & (r.keys <= e))
        if len(inside) < 3:
            inside = np.arange(len(r.keys))
        pick = inside[np.unique(np.linspace(0, len(inside) - 1, min(60, len(inside))).astype(int))]
        self.bg = np.median(crop[pick], axis=0).astype(np.uint8)
        masks = [self.mask(c) for c in crop]
        largest = np.array([self.largest(m) for m in masks])
        self.animal = max(float(np.median(largest[inside])), 0.001 * self.wt * self.ht, 4.0)
        self.floor = self.arena = np.ones((self.ht, self.wt), bool)
        if poly is not None:
            pts = np.column_stack(((poly[:, 0] / clip_scale + x0) * r.fx - self.xa,
                                   (poly[:, 1] / clip_scale + y0) * r.fy - self.ya))
            a = np.zeros((self.ht, self.wt), np.uint8)
            cv2.fillPoly(a, [np.round(pts).astype(np.int32)], 1)
            rad = max(1, int(round(1.5 * np.sqrt(self.animal))))
            self.floor = a > 0
            self.arena = cv2.dilate(a, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * rad + 1,) * 2)) > 0
        # a hand comes from beyond the widened outline, or across the edge of the crop, onto the floor; between the
        # two (walls, their reflections) is neither
        self.outside = ~self.arena
        self.outside[[0, -1], :] = self.outside[:, [0, -1]] = True
        self.key_intrusion = np.array([self.intrusion(c, m) for c, m in zip(crop, masks)])
        self.key_animal = largest >= 0.3 * self.animal

    def mask(self, img):
        # the background at this frame's brightness: a change of exposure (an arm in view) is not foreground
        lit = self.bg > 30
        gain = float(np.median(img[lit] / self.bg[lit].astype(np.float32))) if lit.any() else 1.0
        ref = np.clip(self.bg * gain, 0, 255).astype(np.uint8)
        m = (cv2.absdiff(img, ref) > FG_DIFF).astype(np.uint8)
        return cv2.morphologyEx(m, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))

    @staticmethod
    def largest(m):
        n, _, stats, _ = cv2.connectedComponentsWithStats(m)
        return int(stats[1:, cv2.CC_STAT_AREA].max()) if n > 1 else 0

    def intrusion(self, img, m):
        """Foreground blobs that reach from outside the arena onto its floor, together larger than INTRUSION_ANIMALS
        animals, and that are not a shadow. Blobs only on the floor (the animal, a second animal in a cup) or only
        outside it (people passing, reflections on the walls) do not count."""
        limit = INTRUSION_ANIMALS * self.animal
        if m.sum() <= limit:
            return False
        _, lab, stats, _ = cv2.connectedComponentsWithStats(m)
        touching = np.intersect1d(lab[self.floor & (m > 0)], lab[self.outside & (m > 0)])
        touching = touching[touching > 0]
        if stats[touching, cv2.CC_STAT_AREA].sum() <= limit:
            return False
        px = np.isin(lab, touching)
        r = np.corrcoef(img[px].astype(float), self.bg[px].astype(float))[0, 1]
        return not r >= SHADOW_R

    def frames_intrusion(self, a, b):
        """Intrusion per raw frame in [a, b)."""
        out = []
        for _, f in self.r.frames(a, b - a, width=self.r.tw, gray=True):
            c = f[self.ya:self.ya + self.ht, self.xa:self.xa + self.wt]
            out.append(self.intrusion(c, self.mask(c)))
        return np.array(out, bool)


def trial(sc, s, e, lo, hi):
    """The intrusion-free stretch of the raw that overlaps the clip [s, e] most, within raw frames [lo, hi):
    (release or None, trial end or None when nothing ended it before hi, ended by 'intrusion' / 'raw_end' / '',
    intrusions inside the clip as (first, last) raw key frames). An intrusion inside the clip with the animal in the
    arena before and after it (a hand putting something right) does not end the trial; it is listed instead."""
    r = sc.r
    sel = np.flatnonzero((r.keys >= lo) & (r.keys < hi))
    keys, intr = r.keys[sel], sc.key_intrusion[sel].copy()
    need = int(round(CONFIRM_S * r.fps))
    checked = 0
    for i in np.flatnonzero(intr):   # an intrusion on one key frame only: look at the frames around it
        if (i > 0 and intr[i - 1]) or (i + 1 < len(intr) and intr[i + 1]):
            continue
        if checked >= 12:
            break
        a = keys[i - 1] + 1 if i > 0 else keys[i]
        b = keys[i + 1] if i + 1 < len(keys) else keys[i] + 1
        intr[i] = sc.frames_intrusion(a, b).sum() >= need
        checked += 1
    # stretches between intrusions, as key indices [first, last]
    stretches, start = [], None
    for i, v in enumerate(intr):
        if not v and start is None:
            start = i
        if v and start is not None:
            stretches.append((start, i - 1))
            start = None
    if start is not None:
        stretches.append((start, len(intr) - 1))
    if not stretches:
        return None, None, "", []
    animal, merged, inside = sc.key_animal[sel], [list(stretches[0])], []
    for a, b in stretches[1:]:
        pa, pb = merged[-1]
        if keys[pb + 1] > s and keys[a - 1] < e and animal[pa:pb + 1].mean() >= 0.5 and animal[a:b + 1].mean() >= 0.5:
            inside.append((int(keys[pb + 1]), int(keys[a - 1])))
            merged[-1][1] = b
        else:
            merged.append([a, b])
    stretches = merged
    overlap = [min(e, keys[b]) - max(s, keys[a]) for a, b in stretches]
    a, b = stretches[int(np.argmax(overlap))]
    release = None
    if a > 0:                       # an intrusion ends here: find its last frame
        win = sc.frames_intrusion(keys[a - 1] + 1, keys[a])
        hit = np.flatnonzero(win)
        release = int(keys[a - 1] + 1 + hit[-1] + 1) if len(hit) else int(keys[a - 1] + 1)
    end, by = None, ""
    if b + 1 < len(keys):           # the next intrusion: find its first frame
        win = sc.frames_intrusion(keys[b] + 1, keys[b + 1] + 1)
        hit = np.flatnonzero(win)
        end, by = (int(keys[b] + 1 + hit[0]) if len(hit) else int(keys[b + 1])), "intrusion"
    elif hi >= r.n:
        end, by = r.n, "raw_end"
    return release, end, by, [(p, q) for p, q in inside if keys[a] <= p and q <= keys[b]]


# ---- the check --------------------------------------------------------------------------------------------------------

def run(rec, ctx):
    opts = ctx["opts"]
    clip = Clip(rec, ctx["video"])
    raw_video, raw_dir = core.field(rec, "raw_video"), core.field(rec, "raw_dir")
    if raw_video:
        cands, how = [[p.strip() for p in raw_video.split(";") if p.strip()]], "raw_video"
    elif raw_dir:
        # not the clip itself, nor the other clips next to it (unless the clips sit in raw_dir itself)
        norm = lambda p: os.path.normcase(os.path.abspath(p))  # noqa: E731
        here, clip_dir = norm(rec["video"]), os.path.dirname(norm(rec["video"]))
        own = clip_dir != norm(raw_dir)
        cands = [c for c in recordings_below(raw_dir)
                 if not any(norm(p) == here or (own and os.path.dirname(norm(p)) == clip_dir) for p in c)]
        how = "search"
    else:
        return {"flags": ["raw_not_found"], "error": "no raw_video or raw_dir"}
    exp = core.number(rec, "expected_s")
    m, best, log = find(cands, clip, opts)
    row = {"raw": "", "found_by": "", "align_score": core.rnd(best["score"] if best else None), "tried": log,
           "clip_frames": clip.n}
    if m is None:
        row["flags"] = ["raw_not_found"] + ([] if exp else ["no_expected"])
        if best:
            row["best_try"] = best["r"].name()
        unread = [t[:-len(" unreadable")] for t in log if t.endswith(" unreadable")]
        if not cands:
            row["error"] = "no video files below raw_dir" + ("" if os.path.isdir(raw_dir) else " (no such folder)")
        elif len(unread) == len(cands):
            row["error"] = "cannot read " + ", ".join(unread)
        elif best is not None and best.get("still"):
            row["error"] = "the clip's first seconds show too little movement to place them"
        elif len(unread) == len(log):
            row["error"] = "no candidate is long enough and shows the clip's view"
        elif best is None or not np.isfinite(best["score"]):
            row["error"] = "not placed: the view differs, or the clip holds fewer than 4 of the raw's key frames"
        return row
    r, rate, fps = m["r"], m["rate"], m["r"].fps
    row.update(raw=r.paths, found_by="name" if how == "search" and name_hint(clip.path, r.paths) else how,
               raw_fps=core.rnd(fps, 4), raw_frames=r.n)

    s, e, sq_s, sq_e, box = place_clip(r, clip, m, rate)
    e1 = e + max(1, int(round(rate)))            # the first raw frame after the clip
    row.update(raw_start=s, raw_end=e, start_at=r.at(s), end_at=r.at(e),
               match_start=core.rnd(sq_s, 4), match_end=core.rnd(sq_e, 4),
               raw_span_minus_clip_frames=int(round((e - s) / rate)) + 1 - clip.n,
               crop_x_px=int(round(box[0])), crop_y_px=int(round(box[1])),
               crop_w_px=int(round(box[2])), crop_h_px=int(round(box[3])))

    exp_f = int(round(exp * fps)) if exp else 0
    lo = max(0, s - int(RELEASE_S * fps))
    hi = min(r.n, max(e, s + exp_f) + int(AFTER_S * fps))
    r.extend(hi)
    sc = Scene(r, box, s, e, ctx.get("poly"), clip.scale)
    release, end, by, inside = trial(sc, s, e, lo, hi)
    flags = ["raw_clip_sized"] if m.get("clip_sized") else []
    if release is None:
        flags.append("release_not_seen")
    if inside:
        flags.append("intrusion_in_clip")
        row["intrusion_in_clip_s"] = [f"{(p - s) / fps:.1f}-{(q - s) / fps:.1f}" for p, q in inside]
    ref = release if release is not None else s
    tol = opts.length_tolerance_s * fps
    row.update(release=release if release is not None else "",
               release_s=core.rnd(release / fps, 2) if release is not None else "",
               start_minus_release_s=core.rnd((s - release) / fps, 2) if release is not None else "")
    area = sc.animal / (r.fx * r.fy) * clip.scale ** 2                 # clip px^2
    if ctx.get("scale"):
        row["animal_cm2"] = core.rnd(area * ctx["scale"] ** 2, 1)
    else:
        row["animal_px2"] = core.rnd(area, 0)
    if release is not None:
        if release - s > EARLY_S * fps:
            flags.append("starts_early")
        elif s - release > tol:
            flags.append("starts_late")
    row.update(trial_end=end if end is not None else "", trial_end_by=by,
               next_intrusion_raw_s=core.rnd(end / fps, 2) if by == "intrusion" else "",
               next_intrusion_after_clip_s=core.rnd((end - e1) / fps, 2) if by == "intrusion" else "",
               raw_after_clip_s=core.rnd((r.n - e1) / fps, 2))
    stop = end if end is not None else hi        # the trial lasted at least to here
    if exp:
        x = ref + exp_f
        tail = (r.keys >= e1) & (r.keys < min(x, stop))
        seen = sc.key_animal[tail].mean() if tail.any() else 1.0
        alone = stop >= x and seen >= 0.5
        upto = min(x, stop)
        row.update(expected_end_s=core.rnd(x / fps, 2), end_minus_expected_s=core.rnd((e1 - x) / fps, 2),
                   alone_to_expected="yes" if alone else "no",
                   recoverable_s=core.rnd(max(0, upto - e1) / fps, 2))
        if stop < x - tol:
            flags.append("ended_early")
        if upto - e1 > tol and seen >= 0.5:
            flags.append("cut_short")
        if e1 - upto > tol:
            flags.append("ends_late")
    else:
        flags.append("no_expected")
        if end is not None and e1 - end > tol:
            flags.append("ends_late")
    row["flags"] = flags
    if opts.stills and flags:
        still(rec, opts, clip, r, box, s, e, release, ref + exp_f if exp else None, flags, inside, rate)
    return row


def still(rec, opts, clip, r, box, s, e, release, x, flags, inside, rate):
    """Clip last frame next to the raw at release + expected_s; clip first frame next to the raw at the release when
    the start is flagged; the clip at the first intrusion inside it."""
    W = 360

    def raw_tile(i):
        i = int(min(max(i, 0), r.n - 1))
        crop = even_crop(r, round(box[0]), round(box[1]), box[2], box[3])
        for _, f in r.frames(i, 1, width=W, crop=crop, gray=False):
            return f, i
        return None, i

    def label(img, text):
        if img is None:
            img = np.zeros((W, W, 3), np.uint8)
        img = cv2.resize(img, (W, int(img.shape[0] * W / img.shape[1])))
        cv2.rectangle(img, (0, 0), (W, 22), (0, 0, 0), -1)
        cv2.putText(img, text, (4, 16), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 255), 1)
        return img

    rows = []
    if release is not None and ("starts_early" in flags or "starts_late" in flags):
        f, i = raw_tile(release)
        rows.append([label(clip.fine(0, gray=False, width=W), f"clip first frame = raw {s}"),
                     label(f, f"raw {i}: release ({i / r.fps:.1f} s)")])
    if inside:
        k = int(round((inside[0][0] - s) / rate))
        f, i = raw_tile(inside[0][0])
        rows.append([label(clip.fine(k, gray=False, width=W), f"clip frame {k} = raw {i}"),
                     label(f, f"raw {i}: intrusion ({i / r.fps:.1f} s)")])
    target = x if x is not None else e + 1
    f, i = raw_tile(target)
    if target >= r.n:
        text = f"raw {i} = last; expected end {target / r.fps:.1f} s" if x is not None else f"raw {i} = last"
    else:
        text = f"raw {i} ({i / r.fps:.1f} s): " + ("release + expected" if x is not None else "after the clip")
    rows.append([label(clip.fine(clip.n - 1, gray=False, width=W), f"clip last frame = raw {e}"), label(f, text)])
    h = max(t.shape[0] for row in rows for t in row)
    img = np.vstack([np.hstack([cv2.copyMakeBorder(t, 0, h - t.shape[0], 0, 0, cv2.BORDER_CONSTANT) for t in row])
                     for row in rows])
    bar = np.zeros((26, img.shape[1], 3), np.uint8)
    cv2.putText(bar, f"{rec['id']}: {', '.join(flags)}", (4, 18), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1)
    out = os.path.join(opts.stills_dir, "clip")
    os.makedirs(out, exist_ok=True)
    cv2.imwrite(os.path.join(out, re.sub(r"[^\w.-]+", "_", rec["id"]) + ".jpg"), np.vstack([bar, img]),
                [cv2.IMWRITE_JPEG_QUALITY, 85])
