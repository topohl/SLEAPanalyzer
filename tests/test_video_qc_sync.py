"""tools/video_qc sync check: on synthetic recordings it finds rows whose pose belongs to another frame, a tracked
point off the animal, a tracking shifted against its video, a tracking that belongs to another video and a video
without a visible animal, and passes clean recordings; also a small arena in a larger frame at 25 fps, an mkv without
a frame count, and recordings that cannot be checked (missing video, missing landmark, listed frames past the end)."""
import csv, os, subprocess, sys, tempfile
import numpy as np
import cv2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "video_qc")
sys.path.insert(0, TOOL)
import core  # noqa: E402
PX_PER_CM = 9.0                    # the arena corners are 360 px = 40 cm apart
SPEED = 5.0                        # px per frame: 0.56 cm, a walking mouse


def truth(n, r=120):
    """True body centre and heading of a mouse walking circles around the arena centre (n + 20 frames, so that
    planted poses may come from later frames)."""
    th = np.arange(n + 20) * SPEED / r
    return np.column_stack([200 + r * np.cos(th), 200 + r * np.sin(th)]), th


def make_recording(d, name, n=300, fps=30, animal="dark", visible=True, steps=None, shift=0, off=None, second=False,
                   seed=1, turn=0.0, canvas=None, corners=True, move=(0, 0)):
    """A 400 x 400 clip of one mouse-like ellipse with a tail, static marks on the floor and sensor noise, and its
    tracking. steps: {row: k} gives row the pose of frame row + k; shift: every row gets the pose of frame row + shift;
    off: (first, last) rows whose tracked point lies 4 cm beside the animal; second: a second animal fidgeting just
    outside the path, which the walking one touches and merges with on every lap (a stimulus mouse, a reflection);
    turn: the walk starts this many radians further round; canvas: (width, height, x, y, f) puts the scene, scaled by f,
    at x, y of a larger grey frame; corners: False leaves the arena corners out of the tracking; move: (dx, dy) px added
    to every tracked landmark (a tracking made from a differently cropped clip)."""
    rng = np.random.default_rng(seed)
    pos, th = truth(n)
    if turn:
        rot = np.array([[np.cos(turn), -np.sin(turn)], [np.sin(turn), np.cos(turn)]])
        pos, th = 200 + (pos - 200) @ rot.T, th + turn
    floor, body = (175, 40) if animal == "dark" else (60, 200)
    video = os.path.join(d, name + ".mp4")
    cw, ch, cx, cy, f = canvas or (400, 400, 0, 0, 1.0)
    w = cv2.VideoWriter(video, cv2.VideoWriter_fourcc(*"mp4v"), fps, (cw, ch))
    for i in range(n):
        img = np.full((400, 400, 3), floor, np.uint8)
        cv2.rectangle(img, (60, 300), (90, 330), (110, 110, 110), -1)        # static marks the check must ignore
        cv2.circle(img, (300, 70), 7, (30, 30, 30), -1)
        if second:
            s = np.array([200 + 136, 200 + 35 * np.sin(i / 7.5)])
            cv2.ellipse(img, (int(s[0]), int(s[1])), (38, 15), 90, 0, 360, (body,) * 3, -1)
        if visible:
            c = pos[i]
            heading = th[i] + np.pi / 2
            back = c - 34 * np.array([np.cos(heading), np.sin(heading)])
            tip = c - 75 * np.array([np.cos(heading), np.sin(heading)])
            cv2.line(img, tuple(int(v) for v in back), tuple(int(v) for v in tip), (body,) * 3, 2)
            cv2.ellipse(img, (int(c[0]), int(c[1])), (36, 14), np.degrees(heading), 0, 360, (body,) * 3, -1)
        if canvas:
            big = np.full((ch, cw, 3), 110, np.uint8)
            s = int(round(400 * f))
            big[cy:cy + s, cx:cx + s] = cv2.resize(img, (s, s), interpolation=cv2.INTER_AREA)
            img = big
        img = np.clip(img + rng.normal(0, 3, img.shape), 0, 255).astype(np.uint8)
        w.write(img)
    w.release()

    tracked = pos[:n] + rng.normal(0, 0.8, (n, 2))
    if shift:
        tracked = pos[shift:n + shift].copy()
    for row, k in (steps or {}).items():
        tracked[row] = pos[row + k]
    if off:
        a, b = off
        out = np.column_stack([np.cos(th[a:b + 1]), np.sin(th[a:b + 1])])
        tracked[a:b + 1] = pos[a:b + 1] - 4 * PX_PER_CM * out             # 4 cm towards the arena centre
    tracked = tracked * f + [cx + move[0], cy + move[1]]
    track = os.path.join(d, name + ".csv")
    parts = ["tl", "tr", "br", "bl", "bodycentre"] if corners else ["bodycentre"]
    corners = [(x * f + cx + move[0], y * f + cy + move[1]) for x, y in ((20, 20), (380, 20), (380, 380), (20, 380))
               ] if corners else []
    with open(track, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(["scorer"] * (1 + 3 * len(parts)))
        c.writerow(["bodyparts"] + [p for p in parts for _ in range(3)])
        c.writerow(["coords"] + ["x", "y", "likelihood"] * len(parts))
        for i in range(n):
            row = [i]
            for x, y in corners:
                row += [x, y, 1]
            row += [round(tracked[i, 0], 2), round(tracked[i, 1], 2), 1]
            c.writerow(row)
    return video, track


def run_tool(d, recs, *args):
    """recs: (id, video, tracking, arena_cm, animal[, body]). Returns the sync rows by id and the out-of-step frame
    table."""
    man = os.path.join(d, "manifest.csv")
    with open(man, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(["id", "tracking", "video", "arena_cm", "animal", "body"])
        for rid, video, track, arena, animal, *body in recs:
            c.writerow([rid, track, video, arena, animal, *body])
    out = os.path.join(d, "out")
    subprocess.run([sys.executable, "-I", TOOL, man, out, "--checks", "sync", *args], check=True,
                   stdout=subprocess.DEVNULL)
    with open(os.path.join(out, "video_qc_sync.csv"), newline="", encoding="utf-8") as fh:
        rows = {r["id"]: r for r in csv.DictReader(fh)}
    with open(os.path.join(out, "video_qc_sync_frames.csv"), newline="", encoding="utf-8") as fh:
        frames = list(csv.DictReader(fh))
    return rows, frames, out


def flags(row):
    return [f for f in row["flags"].split(";") if f]


def test_finds_rows_out_of_step_and_passes_clean():
    with tempfile.TemporaryDirectory() as d:
        steps = {40: 5, 100: 4, 101: 4, 102: 4, 103: 4, 200: -5, 201: -5}
        clean = make_recording(d, "clean")
        crowded = make_recording(d, "crowded", second=True)
        bad = make_recording(d, "bad", steps=steps)
        rows, frames, out = run_tool(d, [("clean", *clean, "40x40", "auto"), ("crowded", *crowded, "40x40", "auto"),
                                         ("bad", *bad, "40x40", "dark")], "--sync-windows", "all", "--stills")
        r = rows["clean"]
        assert flags(r) == [] and r["out_of_step"] == "0" and r["window_lag"] == "0", r
        assert float(r["median_dist_cm"]) < 0.3 and r["animal"] == "dark", r
        assert 4.0 < float(r["animal_size_cm"]) < 6.0, r     # the ellipse is 8 x 3.1 cm: 5.0 cm across
        r = rows["crowded"]
        assert flags(r) == [] and r["out_of_step"] == "0", r
        r = rows["bad"]
        assert flags(r) == ["out_of_step"], r
        assert r["out_of_step_frames"] == "40 100-103 200-201", r
        assert r["out_of_step"] == "7" and r["out_of_step_longest"] == "4" and r["out_of_step_offset"] == "4", r
        fits = {int(f["frame"]): int(f["fits_frame"]) for f in frames if f["id"] == "bad"}
        assert fits == {row: row + k for row, k in steps.items()}, fits
        assert all(f["unit"] == "cm" and float(f["dist"]) > 1.5 for f in frames)
        assert os.path.getsize(os.path.join(out, "stills", "sync", "bad.jpg")) > 0
        assert not os.path.exists(os.path.join(out, "stills", "sync", "clean.jpg"))


def test_windows_seek_to_the_right_frames():
    with tempfile.TemporaryDirectory() as d:
        # 30 s; windows of 4 s at the start, the middle and the end: 0-119, 390-509, 780-899
        rec = make_recording(d, "long", n=900, steps={60: 5, 300: 5, 420: 5, 421: 5, 850: -4})
        rows, frames, _ = run_tool(d, [("long", *rec, "40x40", "auto")], "--sync-window-s", "4")
        r = rows["long"]
        assert r["windows"] == "0-119;390-509;780-899", r
        assert r["window_lag"] == "0;0;0", r          # a seek one frame off would show here
        assert r["out_of_step_frames"] == "60 420-421 850" and flags(r) == ["out_of_step"], r


def test_listed_frames_only():
    with tempfile.TemporaryDirectory() as d:
        rec = make_recording(d, "rec", steps={40: 5, 100: 4, 101: 4, 200: -5, 201: -5})
        other = make_recording(d, "other", seed=2)
        listed = os.path.join(d, "candidates.csv")
        with open(listed, "w", newline="") as fh:
            c = csv.writer(fh)
            c.writerow(["id", "frame"]); c.writerow(["rec", "40"]); c.writerow(["rec", "100-101"])
        rows, _, _ = run_tool(d, [("rec", *rec, "40x40", "auto"), ("other", *other, "40x40", "auto")],
                              "--sync-frames", listed)
        r = rows["rec"]
        assert r["windows"] == "0-161", r             # 2 s (60 frames) around each listed row, merged
        assert r["out_of_step_frames"] == "40 100-101", r
        assert rows["other"]["note"] == "no frames listed" and flags(rows["other"]) == [], rows["other"]


def test_off_animal_shifted_and_no_animal():
    with tempfile.TemporaryDirectory() as d:
        off = make_recording(d, "off", off=(100, 199))
        shifted = make_recording(d, "shifted", shift=3)
        empty = make_recording(d, "empty", visible=False)
        rows, _, _ = run_tool(d, [("off", *off, "40x40", "auto"), ("shifted", *shifted, "40x40", "auto"),
                                  ("empty", *empty, "40x40", "auto")], "--sync-windows", "all")
        r = rows["off"]
        assert flags(r) == ["off_animal"] and abs(float(r["share_far"]) - 1 / 3) < 0.03, r
        r = rows["shifted"]
        assert flags(r) == ["shifted"] and r["window_lag"] == "3", r
        r = rows["empty"]
        assert flags(r) == ["no_animal"], r           # the method failed: nothing is said about the tracking


def test_bright_animal_without_scale_in_px():
    with tempfile.TemporaryDirectory() as d:
        rec = make_recording(d, "white", animal="bright", steps={150: 5, 151: 5, 152: 5})
        rows, frames, _ = run_tool(d, [("white", *rec, "", "auto")], "--sync-windows", "all")
        r = rows["white"]
        assert r["animal"] == "bright" and "median_dist_px" in r and "median_dist_cm" not in r, r
        assert r["out_of_step_frames"] == "150-152" and flags(r) == ["out_of_step"], r
        assert all(f["unit"] == "px" and float(f["dist"]) > 15 for f in frames)


def test_tracking_of_another_video_or_crop_is_off_the_animal():
    with tempfile.TemporaryDirectory() as d:
        video, _ = make_recording(d, "a")
        _, other = make_recording(d, "b", turn=np.pi, seed=2)          # the same walk, half a lap away
        _, moved = make_recording(d, "c", move=(45, 0))                 # the same clip cropped 5 cm further left
        rows, _, out = run_tool(d, [("paired", video, other, "40x40", "auto"),
                                    ("moved", video, moved, "40x40", "auto")], "--sync-windows", "all", "--stills")
        r = rows["paired"]
        assert flags(r) == ["off_animal"] and float(r["share_far"]) > 0.9, r     # not no_animal: the video is fine
        assert r["out_of_step"] == "0", r                                        # nothing to be out of step with
        assert r["animal"] == "dark" and "not at the tracked point" in r["note"], r
        r = rows["moved"]
        assert flags(r) == ["off_animal"] and abs(float(r["calib_dx_cm"]) + 5) < 0.5, r
        assert os.path.getsize(os.path.join(out, "stills", "sync", "moved.jpg")) > 0


def test_small_arena_in_a_larger_frame_at_25_fps():
    with tempfile.TemporaryDirectory() as d:
        rec = make_recording(d, "small", fps=25, animal="bright", steps={60: 5, 61: 5, 62: 5, 200: -4},
                             canvas=(640, 480, 330, 180, 0.5), corners=False)
        rows, _, _ = run_tool(d, [("small", *rec, "", "auto")], "--sync-windows", "all")
        r = rows["small"]
        assert r["animal"] == "bright" and "median_dist_px" in r, r
        assert r["out_of_step_frames"] == "60-62 200" and flags(r) == ["out_of_step"], r


def test_container_without_frame_count():
    """mkv gives neither a frame count nor a stream duration; seeking into it must still reach the right frames."""
    with tempfile.TemporaryDirectory() as d:
        video, track = make_recording(d, "long", n=900, steps={60: 5, 420: 5, 421: 5, 850: -4})
        mkv = os.path.join(d, "long.mkv")
        subprocess.run([core.find_tool("ffmpeg"), "-v", "error", "-i", video, "-c", "copy", mkv], check=True)
        rows, _, _ = run_tool(d, [("long", mkv, track, "40x40", "auto")], "--sync-window-s", "4")
        r = rows["long"]
        assert r["windows"] == "0-119;390-509;780-899" and r["window_lag"] == "0;0;0", r
        assert r["out_of_step_frames"] == "60 420-421 850" and flags(r) == ["out_of_step"], r


def test_recordings_that_cannot_be_checked():
    with tempfile.TemporaryDirectory() as d:
        video, track = make_recording(d, "rec")
        listed = os.path.join(d, "late.csv")
        with open(listed, "w", newline="") as fh:
            csv.writer(fh).writerows([["id", "frame"], ["rec", "5000"], ["nose", "40"]])
        rows, _, _ = run_tool(d, [("rec", video, track, "40x40", "auto"),
                                  ("missing", os.path.join(d, "gone.mp4"), track, "40x40", "auto"),
                                  ("nose", video, track, "40x40", "auto", "nose")], "--sync-frames", listed)
        r = rows["rec"]
        assert r["note"] == "no listed frame within the video" and flags(r) == [], r
        assert flags(rows["missing"]) == ["unreadable"] and rows["missing"]["error"], rows["missing"]
        assert flags(rows["nose"]) == ["error"] and "nose" in rows["nose"]["error"], rows["nose"]


if __name__ == "__main__":
    test_finds_rows_out_of_step_and_passes_clean(); test_windows_seek_to_the_right_frames()
    test_listed_frames_only(); test_off_animal_shifted_and_no_animal(); test_bright_animal_without_scale_in_px()
    test_tracking_of_another_video_or_crop_is_off_the_animal(); test_small_arena_in_a_larger_frame_at_25_fps()
    test_container_without_frame_count(); test_recordings_that_cannot_be_checked()
    print("ok")
