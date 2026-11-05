"""tools/qc_review.py renders a still, a clip and an index for each episode of a synthetic recording."""
import csv, os, subprocess, sys, tempfile
import numpy as np
import cv2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "qc_review.py")


def make_recording(d, n=60, fps=30):
    video = os.path.join(d, "rec.mp4")
    w = cv2.VideoWriter(video, cv2.VideoWriter_fourcc(*"mp4v"), fps, (200, 160))
    rows = []
    for i in range(n):
        img = np.full((160, 200, 3), 90, np.uint8)
        x, y = 20 + 2 * i, 80
        cv2.circle(img, (x, y), 8, (20, 20, 20), -1)
        w.write(img)
        rows.append([i, x + 6, y, 1, x, y, 1, x - 8, y, 1, 150, 40, 1])
    w.release()
    track = os.path.join(d, "rec.csv")
    with open(track, "w", newline="") as fh:
        c = csv.writer(fh)
        parts = ["nose", "bodycentre", "tailbase", "socr"]
        c.writerow(["scorer"] * 13)
        c.writerow(["bodyparts"] + [p for p in parts for _ in range(3)])
        c.writerow(["coords"] + ["x", "y", "likelihood"] * 4)
        c.writerows(rows)
    return video, track


def test_renders_still_clip_and_index():
    with tempfile.TemporaryDirectory() as d:
        video, track = make_recording(d)
        eps = os.path.join(d, "episodes.csv"); vids = os.path.join(d, "videos.csv"); out = os.path.join(d, "review")
        with open(eps, "w", newline="") as fh:
            c = csv.writer(fh); c.writerow(["episode_id", "tracking_file", "start_frame", "end_frame", "label", "target", "group"])
            c.writerow(["ep01", track, 21, 40, "test episode", "socr", "synthetic"])
        with open(vids, "w", newline="") as fh:
            c = csv.writer(fh); c.writerow(["tracking_file", "video_path"]); c.writerow([track, video])
        subprocess.run([sys.executable, "-I", TOOL, eps, vids, out, "--pad-s", "0.2", "--width", "200"], check=True)
        assert os.path.getsize(os.path.join(out, "stills", "ep01.jpg")) > 0
        cap = cv2.VideoCapture(os.path.join(out, "clips", "ep01.mp4"))
        frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT)); cap.release()
        assert frames == 20 + 2 * 6          # the episode plus 0.2 s on each side
        index = open(os.path.join(out, "review_index.csv")).read()
        assert "ep01" in index and "clips" in index
        assert "ep01" in open(os.path.join(out, "index.html"), encoding="utf-8").read()


def test_repaired_file_is_drawn():
    with tempfile.TemporaryDirectory() as d:
        video, track = make_recording(d)
        fixed = os.path.join(d, "fixed.csv")
        lines = open(track).read().splitlines()
        rows = [l.split(",") for l in lines[3:]]
        for r in rows:
            r[2] = str(float(r[2]) + 40)                      # the repaired nose 40 px below the original
            r[4] = str(float(r[4]) - 40)                      # the repaired body centre 40 px to the left
        open(fixed, "w").write("\n".join(lines[:3] + [",".join(r) for r in rows]) + "\n")
        eps = os.path.join(d, "episodes.csv"); vids = os.path.join(d, "videos.csv")
        with open(vids, "w", newline="") as fh:
            c = csv.writer(fh); c.writerow(["tracking_file", "video_path"]); c.writerow([track, video])
        stills = {}
        for name, rep in (("plain", ""), ("repaired", fixed)):
            with open(eps, "w", newline="") as fh:
                c = csv.writer(fh); c.writerow(["episode_id", "tracking_file", "start_frame", "end_frame", "label", "repaired_file"])
                c.writerow(["ep01", track, 21, 40, "test", rep])
            out = os.path.join(d, name)
            subprocess.run([sys.executable, "-I", TOOL, eps, vids, out, "--pad-s", "0", "--width", "200"], check=True)
            stills[name] = cv2.imread(os.path.join(out, "stills", "ep01.jpg"))
        diff = np.abs(stills["repaired"].astype(int) - stills["plain"].astype(int)).sum(axis=2)
        rows_changed = np.where(diff.max(axis=1) > 60)[0]
        # the repaired nose ring (y = 120) and body centre ring (y = 80, 40 px left of the original), radius 12
        assert any(66 <= r <= 94 for r in rows_changed) and any(106 <= r <= 136 for r in rows_changed)
        assert 60 <= rows_changed.min() and rows_changed.max() <= 137


if __name__ == "__main__":
    test_renders_still_clip_and_index(); test_repaired_file_is_drawn(); print("ok")
