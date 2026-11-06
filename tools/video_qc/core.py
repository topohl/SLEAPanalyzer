"""Shared pieces of the video QC: the manifest, tracking files, video probing and decoding, scale, output.

Nothing here knows about a particular experiment. A recording is one manifest row:

  id            any unique name (required)
  tracking      tracking csv whose rows are the video frames, DLC/SLEAP style: three header lines (scorer,
                bodyparts, coords), then frame, x, y, likelihood per landmark (required)
  video         the clip the tracking was made from (a row without one is reported unreadable)
  raw_video     the recording the clip was cut from; several files of one recording (camera chapters) in
                order, separated by ';' (optional; the clip check)
  raw_dir       a folder to search for the raw recording when raw_video is not given (optional)
  clip_scale    clip pixels per raw pixel, when the clip was resized after cropping (default 1; the clip check)
  expected_s    protocol duration of the trial in seconds (optional; frames and clip checks)
  group         recordings expected to look alike, e.g. a batch or a camera setup (optional; group checks)
  corners       arena corner landmarks in order around the arena, default "tl tr br bl"
  outline       landmarks around the arena floor in order, when they are not just the corners (e.g. the twelve
                corners of a plus maze); masks use it, the scale uses corners
  arena_cm      real arena size: "W x H" (area of the four corners) or one length (the first two corners)
  scale_cm_per_px  a known scale instead of arena_cm
  objects       landmarks of fixed objects to check against the video, e.g. "objL objR" (optional)
  object_set    recordings with the same objects, e.g. one test and trial: the objects check compares groups only
                within it (optional)
  body          the landmark that follows the animal, default "bodycentre"
  animal        "dark", "bright" or "auto" (default): the animal against the floor

Any other column is carried through to the outputs.
"""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
from fractions import Fraction

import numpy as np


# ---- tools ------------------------------------------------------------------------------------------------------------

def find_tool(name: str, explicit: str | None = None) -> str:
    """ffmpeg or ffprobe: an explicit path, the FFMPEG/FFPROBE variable, PATH, imageio-ffmpeg, or next to Python."""
    cands = [explicit, os.environ.get(name.upper())]
    cands.append(shutil.which(name))
    if name == "ffmpeg":
        try:
            import imageio_ffmpeg
            cands.append(imageio_ffmpeg.get_ffmpeg_exe())
        except ImportError:
            pass
    exe = name + (".exe" if os.name == "nt" else "")
    for base in (sys.prefix, os.path.dirname(sys.executable)):
        cands += [os.path.join(base, "Library", "bin", exe), os.path.join(base, "bin", exe), os.path.join(base, exe)]
    for c in cands:
        if c and os.path.isfile(c):
            return c
    raise SystemExit(f"{name} not found: pass --{name}, set {name.upper()}, or put it on PATH")


# ---- manifest ---------------------------------------------------------------------------------------------------------

def read_manifest(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    for k, r in enumerate(rows):
        for need in ("id", "tracking"):
            if not (r.get(need) or "").strip():
                raise SystemExit(f"manifest row {k + 2}: '{need}' is empty")
    ids = [r["id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise SystemExit("manifest ids are not unique")
    return rows


def field(rec: dict, name: str, default=None):
    v = (rec.get(name) or "").strip()
    return v if v else default


def names(rec: dict, name: str, default: str = "") -> list[str]:
    return [s for s in (field(rec, name, default) or "").replace(",", " ").split() if s]


def number(rec: dict, name: str):
    v = field(rec, name)
    try:
        return float(v) if v is not None else None
    except ValueError:
        return None


# ---- tracking ---------------------------------------------------------------------------------------------------------

def read_tracking(path: str) -> tuple[dict, int]:
    """landmark -> (x, y) arrays in pixels, and the number of rows (= video frames)."""
    with open(path, newline="") as fh:
        fh.readline()
        parts = fh.readline().rstrip("\r\n").split(",")
        coords = fh.readline().rstrip("\r\n").split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=3)
    if data.ndim == 1:
        data = data[None, :]
    out = {}
    for name in dict.fromkeys(parts[1:]):
        ix = [i for i, p in enumerate(parts) if p == name]
        xi = [i for i in ix if coords[i] == "x"]
        yi = [i for i in ix if coords[i] == "y"]
        if xi and yi:
            out[name] = (data[:, xi[0]].astype(float), data[:, yi[0]].astype(float))
    return out, data.shape[0]


def median_point(track: dict, name: str):
    if name not in track:
        return None
    x, y = track[name]
    if not np.isfinite(x).any():
        return None
    return np.array([np.nanmedian(x), np.nanmedian(y)])


def arena_polygon(track: dict, corners: list[str]):
    """Median positions of the outline landmarks, in order (pixels), or None."""
    pts = [median_point(track, c) for c in corners]
    if len(pts) < 3 or any(p is None for p in pts):
        return None
    return np.array(pts)


def polygon_area(poly) -> float:
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def cm_per_px(rec: dict, track: dict):
    """The scale of this recording, or None (results then stay in pixels)."""
    s = number(rec, "scale_cm_per_px")
    if s:
        return s
    size = field(rec, "arena_cm")
    corners = names(rec, "corners", "tl tr br bl")
    if not size:
        return None
    size = size.lower().replace("×", "x").replace(" ", "")
    if "x" in size:
        w, h = (float(v) for v in size.split("x"))
        poly = arena_polygon(track, corners)
        return math.sqrt(w * h / polygon_area(poly)) if poly is not None and len(poly) >= 3 else None
    a, b = median_point(track, corners[0]), median_point(track, corners[1]) if len(corners) > 1 else None
    if a is None or b is None:
        return None
    return float(size) / float(np.hypot(*(a - b)))


# ---- video ------------------------------------------------------------------------------------------------------------

class Video:
    """Probe and decode one video with ffprobe/ffmpeg. Frames are numbered from 0 in decoding order."""

    def __init__(self, path: str, ffmpeg: str, ffprobe: str):
        if not path:
            raise ValueError("no video in the manifest")
        self.path, self.ffmpeg, self.ffprobe = path, ffmpeg, ffprobe
        out = subprocess.run([ffprobe, "-v", "error", "-select_streams", "v:0", "-show_entries",
                              "stream=width,height,r_frame_rate,avg_frame_rate,nb_frames,duration:format=duration",
                              "-of", "json", path], capture_output=True, text=True)
        if out.returncode != 0:
            raise RuntimeError(f"ffprobe failed on {path}: {out.stderr.strip()[:200]}")
        probe = json.loads(out.stdout or "{}")
        if not probe.get("streams"):
            raise RuntimeError(f"no video stream in {path}")
        info = {k: str(v) for k, v in probe["streams"][0].items()}
        self.width, self.height = int(info.get("width", 0)), int(info.get("height", 0))
        rate = info.get("avg_frame_rate") or "0/0"
        if rate in ("0/0", "0/1"):
            rate = info.get("r_frame_rate") or "0/0"
        self.fps = float(Fraction(rate)) if rate not in ("0/0", "0/1") else float("nan")
        nominal = info.get("r_frame_rate") or "0/0"
        self.nominal_fps = float(Fraction(nominal)) if nominal not in ("0/0", "0/1") else float("nan")
        nb = info.get("nb_frames", "N/A")
        self.nb_frames = int(nb) if nb.isdigit() and int(nb) > 0 else None
        # some containers (mkv, webm) give no stream duration or frame count, only the container's duration
        self.duration = float("nan")
        for d in (info.get("duration"), (probe.get("format") or {}).get("duration")):
            try:
                self.duration = float(d)
                break
            except (TypeError, ValueError):
                pass

    def count_frames(self) -> int:
        """Exact frame count by decoding (slow); for containers without a frame count."""
        out = subprocess.run([self.ffprobe, "-v", "error", "-select_streams", "v:0", "-count_frames", "-show_entries",
                              "stream=nb_read_frames", "-of", "csv=p=0", self.path], capture_output=True, text=True)
        return int(out.stdout.strip().split(",")[0])

    def frames(self, start: int = 0, n: int | None = None, width: int | None = None, crop=None, gray: bool = True):
        """Yield decoded frames start, start+1, ... (n of them, or to the end), optionally cropped (x, y, w, h in source
        pixels) and scaled to `width`. Starting at 0 decodes from the first frame; later starts seek by time, which is
        frame-exact for constant-rate video."""
        vf = []
        if crop is not None:
            x, y, w, h = (int(round(v)) for v in crop)
            vf.append(f"crop={w}:{h}:{x}:{y}:exact=1")
            src_w, src_h = w, h
        else:
            src_w, src_h = self.width, self.height
        if width:
            height = int(round(src_h * width / src_w / 2) * 2)
            vf.append(f"scale={width}:{height}:flags=area")
            out_w, out_h = width, height
        else:
            out_w, out_h = src_w, src_h
        pix = "gray" if gray else "bgr24"
        vf.append(f"format={pix}")
        cmd = [self.ffmpeg, "-v", "error"]
        if start > 0:
            cmd += ["-ss", f"{start / self.fps:.6f}"]
        cmd += ["-i", self.path]
        if n is not None:
            cmd += ["-frames:v", str(int(n))]
        cmd += ["-vf", ",".join(vf), "-f", "rawvideo", "-pix_fmt", pix, "-"]
        size = out_w * out_h * (1 if gray else 3)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        try:
            while True:
                buf = p.stdout.read(size)
                if len(buf) < size:
                    break
                a = np.frombuffer(buf, np.uint8)
                yield a.reshape(out_h, out_w) if gray else a.reshape(out_h, out_w, 3)
        finally:
            p.stdout.close()
            p.wait()

    def frame(self, index: int, width: int | None = None, gray: bool = False):
        for f in self.frames(index, 1, width=width, gray=gray):
            return f
        return None

    def background(self, n: int = 15, width: int | None = None, gray: bool = False, margin_frames: int = 60):
        """Median of n frames spread over the video: the static scene without the moving animal."""
        total = self.nb_frames
        if not total:
            total = int(self.duration * self.fps) if math.isfinite(self.duration * self.fps) else self.count_frames()
        lo, hi = min(margin_frames, total // 4), max(total - margin_frames, total // 4 * 3)
        picks = np.linspace(lo, hi, n).astype(int)
        frames = [f for f in (self.frame(int(i), width=width, gray=gray) for i in picks) if f is not None]
        if not frames:
            return None
        return np.median(np.stack(frames), axis=0).astype(np.uint8)


# ---- output -----------------------------------------------------------------------------------------------------------

def read_table(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(rows: list[dict], path: str):
    """Rows to csv, the union of their keys as columns; keys starting with '_' (a check's private data) are left out."""
    keys = []
    for r in rows:
        for k in r:
            if k not in keys and not k.startswith("_"):
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (";".join(v) if isinstance(v, (list, tuple)) else v) for k, v in r.items() if k in keys})


def rnd(v, d=3):
    try:
        return round(float(v), d) if v is not None and np.isfinite(v) else ""
    except (TypeError, ValueError):
        return v
