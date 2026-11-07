"""tools/video_qc review sheets: files written, rows skipped and cleared, and planted faults (a corner jump, a pose
spike, nose and tail swapped while moving) in the frames shown."""
import csv, os, subprocess, sys, tempfile
import numpy as np
import cv2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "video_qc", "review.py")
sys.path.insert(0, os.path.join(ROOT, "tools", "video_qc"))
import core  # noqa: E402

N, W, H = 300, 200, 160
CORNERS = {"tl": (20, 20), "tr": (180, 20), "br": (180, 140), "bl": (20, 140)}   # 160 x 120 px = 32 x 24 cm
JUMP, SPIKE, SWAP = range(200, 205), 100, range(250, 260)


def body_path():
    """A mouse running in a circle (about 19 cm/s), nose ahead and tail base behind."""
    t = 2 * np.pi * np.arange(N) / 100
    body = np.column_stack([100 + 50 * np.sin(t), 80 + 30 * np.cos(t)])
    v = np.column_stack([np.cos(t), -0.6 * np.sin(t)])
    v /= np.hypot(*v.T)[:, None]
    return body, body + 8 * v, body - 8 * v


def make_files(d):
    body, nose, tail = body_path()
    video = os.path.join(d, "clip.mp4")
    w = cv2.VideoWriter(video, cv2.VideoWriter_fourcc(*"mp4v"), 30, (W, H))
    for i in range(N):
        img = np.full((H, W, 3), 170, np.uint8)
        cv2.rectangle(img, CORNERS["tl"], CORNERS["br"], (255, 255, 255), 2)
        cv2.circle(img, (60, 80), 6, (40, 120, 200), -1)
        cv2.ellipse(img, tuple(int(v) for v in body[i]), (9, 5), 0, 0, 360, (30, 30, 30), -1)
        w.write(img)
    w.release()
    s = slice(SWAP.start, SWAP.stop)
    nose[s], tail[s] = tail[s].copy(), nose[s].copy()                # nose and tail swapped while it runs
    pts = {k: np.tile(np.array(p, float), (N, 1)) for k, p in CORNERS.items()}
    pts["tr"][JUMP.start:JUMP.stop] = (140, 60)                     # the corner jumps onto the floor
    pts["objL"] = np.tile([60.0, 80.0], (N, 1))
    for m in (nose, body, tail):                                     # a pose from elsewhere, across the motion
        m[SPIKE] += (0, 15)
    pts.update(nose=nose, bodycentre=body, tailbase=tail)
    tracking = os.path.join(d, "track.csv")
    with open(tracking, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(["scorer"] * (1 + 3 * len(pts)))
        c.writerow(["bodyparts"] + [k for k in pts for _ in range(3)])
        c.writerow(["coords"] + ["x", "y", "likelihood"] * len(pts))
        for i in range(N):
            c.writerow([i] + [v for k in pts for v in (pts[k][i, 0], pts[k][i, 1], 1)])
    man = os.path.join(d, "manifest.csv")
    with open(man, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(["id", "tracking", "video", "expected_s", "arena_cm", "objects"])
        for rid in ("rec one", "rec two", "rec three"):
            c.writerow([rid, tracking, video, 10, "32x24", "objL"])
    return man


def run(d, man, rows, *args):
    lst = os.path.join(d, "list.csv")
    with open(lst, "w", newline="") as fh:
        c = csv.writer(fh); c.writerow(["id", "flags", "draw", "pair"]); c.writerows(rows)
    out = os.path.join(d, "out")
    subprocess.run([sys.executable, "-I", TOOL, man, lst, out] + list(args), check=True, capture_output=True)
    return out, {r["id"]: r for r in core.read_table(os.path.join(out, "review_sheets.csv"))}


def shown(row):
    return {p.split(":")[0]: [int(v) for v in p.split(":")[1].split()] for p in row["frames"].split("; ")}


def test_sheets():
    with tempfile.TemporaryDirectory() as d:
        man = make_files(d)
        out, idx = run(d, man, [
            ["rec one", "corner; length; readout_outlier",
             "corner (tracking changed since its decision); pose_spikes; nose_tail_order; duplicate", "rec two"],
            ["rec two", "readout_outlier", "readout_outlier", ""],
            ["rec three", "corner", "", ""]])                       # nothing to draw: skipped
        assert sorted(idx) == ["rec one", "rec two"] and not idx["rec one"]["error"], idx
        img = cv2.imread(os.path.join(out, "rec_one.jpg"))
        assert img is not None and img.shape[1] > 1000 and img.shape[0] > 600, None if img is None else img.shape
        assert not os.path.exists(os.path.join(out, "rec_three.jpg"))
        assert idx["rec one"]["drawn"] == "corner;pose_spikes;nose_tail_order;duplicate", idx["rec one"]
        f = shown(idx["rec one"])
        assert f["clip"][0] == 0 and f["clip"][-1] == N - 1, f
        assert len(f["corner"]) == 1 and f["corner"][0] in JUMP, f        # one jump, in the frames it happens
        assert f["pose_spikes"] == [SPIKE - 1, SPIKE, SPIKE + 1], f
        assert f["nose_tail_order"] and all(i in SWAP for i in f["nose_tail_order"]), f
        assert len(f["duplicate"]) == 3, f
        assert list(shown(idx["rec two"])) == ["clip"], idx["rec two"]

        # a later run without rec one removes its sheet; --only keeps the other rows
        out, idx = run(d, man, [["rec two", "corner", "corner", ""]])
        assert list(idx) == ["rec two"] and not os.path.exists(os.path.join(out, "rec_one.jpg"))
        out, idx = run(d, man, [["rec one", "length", "length", ""]], "--only", "one")
        assert sorted(idx) == ["rec one", "rec two"] and os.path.exists(os.path.join(out, "rec_two.jpg")), idx
        assert shown(idx["rec one"])["length"][0] == 0 and shown(idx["rec one"])["length"][-1] == N - 1


if __name__ == "__main__":
    test_sheets(); print("ok")
