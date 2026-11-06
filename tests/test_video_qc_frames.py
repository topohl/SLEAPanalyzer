"""tools/video_qc frames check and command line: rows against frames, length, rate, containers, and how runs into one
folder combine."""
import csv, os, subprocess, sys, tempfile
import numpy as np
import cv2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "video_qc")
sys.path.insert(0, TOOL)
import core  # noqa: E402


def make_clip(d, name, n, fps=30):
    path = os.path.join(d, name + ".mp4")
    w = cv2.VideoWriter(path, cv2.VideoWriter_fourcc(*"mp4v"), fps, (160, 120))
    for i in range(n):
        img = np.full((120, 160, 3), 150, np.uint8)
        cv2.circle(img, (20 + i % 120, 60), 6, (20, 20, 20), -1)
        w.write(img)
    w.release()
    return path


def make_tracking(d, name, n):
    path = os.path.join(d, name + ".csv")
    parts = ["bodycentre", "tl", "tr", "br", "bl"]
    with open(path, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(["scorer"] * (1 + 3 * len(parts)))
        c.writerow(["bodyparts"] + [p for p in parts for _ in range(3)])
        c.writerow(["coords"] + ["x", "y", "likelihood"] * len(parts))
        for i in range(n):
            c.writerow([i, 20 + i % 120, 60, 1, 5, 5, 1, 155, 5, 1, 155, 115, 1, 5, 115, 1])
    return path


def run(d, rows, checks="frames", table="video_qc_frames.csv"):
    man = os.path.join(d, "manifest.csv")
    with open(man, "w", newline="") as fh:
        c = csv.writer(fh); c.writerow(["id", "tracking", "video", "expected_s", "group"]); c.writerows(rows)
    out = os.path.join(d, "out")
    subprocess.run([sys.executable, "-I", TOOL, man, out, "--checks", checks], check=True, capture_output=True)
    with open(os.path.join(out, table), newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}


def test_rows_and_length():
    with tempfile.TemporaryDirectory() as d:
        ok_v, ok_t = make_clip(d, "ok", 300), make_tracking(d, "ok", 300)
        short_v, short_t = make_clip(d, "short", 200), make_tracking(d, "short", 200)
        lost_t = make_tracking(d, "lost", 290)                      # tracking that lost 10 frames of the ok clip
        res = run(d, [["ok", ok_t, ok_v, 10, "a"], ["short", short_t, short_v, 10, "a"], ["lost", lost_t, ok_v, 10, "a"],
                      ["missing", ok_t, os.path.join(d, "none.mp4"), 10, "a"]])
        assert res["ok"]["flags"] == "", res["ok"]
        assert res["ok"]["video_frames"] == "300" and res["ok"]["rows"] == "300"
        assert res["short"]["flags"] == "length" and float(res["short"]["length_minus_expected_s"]) < -3
        assert res["lost"]["flags"] == "rows_differ" and res["lost"]["rows_minus_frames"] == "-10"
        assert res["missing"]["flags"] == "unreadable" and res["missing"]["error"]


def test_later_run_keeps_earlier_checks():
    """A run of another check into the same folder keeps the frames flags in video_qc.csv and names both checks; a
    row without a video is unreadable, not a refused manifest."""
    with tempfile.TemporaryDirectory() as d:
        short_v, short_t = make_clip(d, "short", 200), make_tracking(d, "short", 200)
        rows = [["short", short_t, short_v, 10, "a"], ["novideo", short_t, "", 10, "a"]]
        run(d, rows)
        res = run(d, rows, checks="objects", table="video_qc.csv")
        assert res["short"]["checks"] == "frames;objects" and res["short"]["flags"] == "frames:length", res["short"]
        assert res["novideo"]["flags"] == "frames:unreadable;objects:unreadable", res["novideo"]


def test_only_keeps_the_other_rows():
    """A rerun of a few recordings with --only into the same folder replaces their rows and keeps the others."""
    with tempfile.TemporaryDirectory() as d:
        ok_v, ok_t = make_clip(d, "ok", 300), make_tracking(d, "ok", 300)
        lost_t = make_tracking(d, "lost", 290)
        rows = [["ok", ok_t, ok_v, 10, "a"], ["lost", lost_t, ok_v, 10, "a"]]
        run(d, rows)
        rows[1][1] = ok_t                                            # the tracking of 'lost' corrected
        man = os.path.join(d, "manifest.csv")
        with open(man, "w", newline="") as fh:
            c = csv.writer(fh); c.writerow(["id", "tracking", "video", "expected_s", "group"]); c.writerows(rows)
        out = os.path.join(d, "out")
        subprocess.run([sys.executable, "-I", TOOL, man, out, "--checks", "frames", "--only", "^lost$"], check=True,
                       capture_output=True)
        tab = {r["id"]: r for r in core.read_table(os.path.join(out, "video_qc_frames.csv"))}
        both = {r["id"]: r for r in core.read_table(os.path.join(out, "video_qc.csv"))}
        assert list(tab) == ["ok", "lost"] and tab["lost"]["flags"] == "" and tab["ok"]["video_frames"] == "300"
        assert list(both) == ["ok", "lost"] and both["lost"]["flags"] == "", both
        r = subprocess.run([sys.executable, "-I", TOOL, man, out, "--checks", "frames", "--only", "none"],
                           check=True, capture_output=True, text=True)
        assert "nothing written" in r.stdout and len(core.read_table(os.path.join(out, "video_qc.csv"))) == 2


def test_containers_and_rates():
    """An mkv without a stream duration or frame count is read through the container; a variable frame rate is
    flagged; a check's private keys stay out of its table."""
    ffmpeg = core.find_tool("ffmpeg")
    with tempfile.TemporaryDirectory() as d:
        ok_v, ok_t = make_clip(d, "ok", 300), make_tracking(d, "ok", 300)
        mkv, vfr = os.path.join(d, "ok.mkv"), os.path.join(d, "vfr.mp4")
        subprocess.run([ffmpeg, "-v", "error", "-y", "-i", ok_v, "-c", "copy", mkv], check=True)
        subprocess.run([ffmpeg, "-v", "error", "-y", "-i", ok_v, "-vf", "setpts='if(lt(N,150),N,N*2-150)/30/TB'",
                        "-vsync", "passthrough", "-c:v", "mpeg4", vfr], check=True)
        v = core.Video(mkv, ffmpeg, core.find_tool("ffprobe"))
        assert v.nb_frames is None and abs(v.duration - 10) < 0.1 and v.background(5, width=80) is not None
        res = run(d, [["mkv", ok_t, mkv, 10, "a"], ["vfr", ok_t, vfr, 10, "a"]])
        assert res["mkv"]["flags"] == "" and res["mkv"]["video_frames"] == "300", res["mkv"]
        assert "rate" in res["vfr"]["flags"].split(";") and res["vfr"]["nominal_fps"] == "30.0", res["vfr"]
        path = os.path.join(d, "t.csv")
        core.write_csv([{"id": "a", "_still": b"\xff", "flags": []}], path)
        assert core.read_table(path) == [{"id": "a", "flags": ""}]


if __name__ == "__main__":
    test_rows_and_length(); test_later_run_keeps_earlier_checks(); test_only_keeps_the_other_rows()
    test_containers_and_rates(); print("ok")
