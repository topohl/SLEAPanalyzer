"""Render review material for QC episodes: one still strip and one short clip per episode, plus index.html.

Works for any assay whose tracking file rows are the video frames (the staged SLEAPanalyzer inputs).

usage:
  python -I tools/qc_review.py <episodes.csv> <videos.csv> <out_dir> [--pad-s 1] [--max-clip-s 45] [--width 540]

episodes.csv  one row per episode. Required columns: episode_id, tracking_file, start_frame, end_frame (1-based
              positions in the tracking file, as written by SLEAPanalyzer), label (free text shown on the clip).
              Optional: target (a landmark drawn as a cross, e.g. socr), group (index section), note.
videos.csv    tracking_file, video_path: the video each tracking file was made from (rows = video frames).

Markers: red = nose, green = body centre, blue = tail base, yellow cross = target landmark. Frames inside the
episode carry an orange border and the label; the padding before and after has no border.
"""
import argparse, csv, html, os, subprocess, sys
import numpy as np
import cv2

try:
    import imageio_ffmpeg
    FFMPEG = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:  # fall back to ffmpeg on PATH
    FFMPEG = "ffmpeg"


def read_tracking(path):
    with open(path) as fh:
        fh.readline(); parts = fh.readline().strip().split(","); coords = fh.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=3)
    def get(name):
        ix = [i for i, p in enumerate(parts) if p == name]
        if not ix:
            return None
        xi = [i for i in ix if coords[i] == "x"][0]; yi = [i for i in ix if coords[i] == "y"][0]
        return data[:, xi], data[:, yi]
    return get, data.shape[0]


def annotate(img, get_pts, frame, target, inside, text):
    for name, colour in (("nose", (0, 0, 255)), ("bodycentre", (0, 255, 0)), ("tailbase", (255, 0, 0))):
        pts = get_pts.get(name)
        if pts is None or frame >= len(pts[0]):
            continue
        x, y = pts[0][frame], pts[1][frame]
        if np.isfinite(x) and np.isfinite(y):
            cv2.circle(img, (int(x), int(y)), 7, colour, 2)
    if target is not None:
        cv2.drawMarker(img, target, (0, 255, 255), cv2.MARKER_CROSS, 20, 2)
    h, w = img.shape[:2]
    if inside:
        cv2.rectangle(img, (0, 0), (w - 1, h - 1), (0, 140, 255), 12)
    cv2.rectangle(img, (0, 0), (w, 44), (0, 0, 0), -1)
    cv2.putText(img, text, (12, 31), cv2.FONT_HERSHEY_SIMPLEX, 0.9, (255, 255, 255), 2)
    return img


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("episodes"); ap.add_argument("videos"); ap.add_argument("out")
    ap.add_argument("--pad-s", type=float, default=1.0)
    ap.add_argument("--max-clip-s", type=float, default=45.0)
    ap.add_argument("--width", type=int, default=540)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "stills"), exist_ok=True)
    os.makedirs(os.path.join(a.out, "clips"), exist_ok=True)
    videos = {r["tracking_file"]: r["video_path"] for r in csv.DictReader(open(a.videos))}
    episodes = list(csv.DictReader(open(a.episodes)))
    cache = {}
    rows = []
    for ep in episodes:
        tf = ep["tracking_file"]
        if tf not in videos:
            sys.exit(f"no video for {tf}")
        if tf not in cache:
            get, n = read_tracking(tf)
            pts = {k: get(k) for k in ("nose", "bodycentre", "tailbase")}
            cache[tf] = (get, n, pts)
        get, n, pts = cache[tf]
        target = None
        if ep.get("target"):
            t = get(ep["target"])
            if t is not None:
                target = (int(np.nanmedian(t[0])), int(np.nanmedian(t[1])))
        cap = cv2.VideoCapture(videos[tf])
        fps = cap.get(cv2.CAP_PROP_FPS) or 30.0
        first, last = int(ep["start_frame"]) - 1, int(ep["end_frame"]) - 1   # 0-based video frames
        pad = int(round(a.pad_s * fps))
        c0, c1 = max(0, first - pad), min(n - 1, last + pad)
        if (c1 - c0 + 1) / fps > a.max_clip_s:      # long episode: keep the start and the end
            keep = int(a.max_clip_s * fps / 2)
            spans = [(c0, c0 + keep), (c1 - keep, c1)]
        else:
            spans = [(c0, c1)]
        eid = ep["episode_id"]
        clip = os.path.join("clips", f"{eid}.mp4")
        out_w = a.width; out_h = None
        proc = None
        stills = []
        want = sorted({first, first + (last - first) // 3, first + 2 * (last - first) // 3, last})
        for s0, s1 in spans:
            cap.set(cv2.CAP_PROP_POS_FRAMES, s0)
            for f in range(s0, s1 + 1):
                ok, img = cap.read()
                if not ok:
                    break
                inside = first <= f <= last
                text = f"{eid}  t={f / fps:6.1f}s  {ep['label'] if inside else ''}"
                img = annotate(img, pts, f, target, inside, text)
                if f in want:
                    stills.append(img.copy())
                if proc is None:
                    out_h = int(round(img.shape[0] * out_w / img.shape[1] / 2) * 2)
                    proc = subprocess.Popen([FFMPEG, "-loglevel", "error", "-y", "-f", "rawvideo", "-pix_fmt", "bgr24",
                                             "-s", f"{out_w}x{out_h}", "-r", f"{fps}", "-i", "-", "-c:v", "libx264",
                                             "-pix_fmt", "yuv420p", "-crf", "23", "-preset", "veryfast",
                                             os.path.join(a.out, clip)], stdin=subprocess.PIPE)
                proc.stdin.write(cv2.resize(img, (out_w, out_h)).tobytes())
        cap.release()
        if proc is not None:
            proc.stdin.close(); proc.wait()
        if stills:
            tiles = [cv2.resize(s, (out_w, out_h)) for s in stills]
            cv2.imwrite(os.path.join(a.out, "stills", f"{eid}.jpg"), np.hstack(tiles), [cv2.IMWRITE_JPEG_QUALITY, 85])
        rows.append(dict(ep, still=os.path.join("stills", f"{eid}.jpg"), clip=clip,
                         clip_note="start and end only" if len(spans) > 1 else ""))
        print("rendered", eid, flush=True)
    with open(os.path.join(a.out, "review_index.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    groups = {}
    for r in rows:
        groups.setdefault(r.get("group", "") or "episodes", []).append(r)
    parts = ["<!doctype html><meta charset='utf-8'><title>QC review</title>",
             "<style>body{font-family:sans-serif;margin:16px} img{width:100%;max-width:1400px} "
             "td{vertical-align:top;padding:6px;border-bottom:1px solid #ccc}</style>",
             "<p>Red = nose, green = body centre, blue = tail base, yellow cross = cup or object point. "
             "Orange border = frames inside the episode. Stills: first, one-third, two-thirds and last episode frame.</p>"]
    for g, rs in groups.items():
        parts.append(f"<h2>{html.escape(g)} ({len(rs)})</h2><table>")
        for r in rs:
            meta = html.escape(" | ".join(f"{k}: {r[k]}" for k in r if k not in ("still", "clip", "tracking_file", "group") and r[k]))
            parts.append(f"<tr><td><b>{html.escape(r['episode_id'])}</b><br>{meta}<br>"
                         f"<a href='{r['clip']}'>clip</a></td><td><img src='{r['still']}'></td></tr>")
        parts.append("</table>")
    with open(os.path.join(a.out, "index.html"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts))


if __name__ == "__main__":
    main()
