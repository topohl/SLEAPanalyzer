"""tools/video_qc runs: a second run reuses every row and writes the same tables; a changed tracking file is checked
again alone (the sync frames of the others stay, stills follow the flags); --fresh and --only; rows of an earlier
version taken over with --adopt; a run killed part way is continued by the next; the objects group comparison after a
partial re-check is the one of a fresh run, stills included."""
import csv, os, re, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(os.path.dirname(HERE), "tools", "video_qc")
sys.path.insert(0, HERE)
import test_video_qc_frames as fra   # noqa: E402  synthetic recordings of the other tests
import test_video_qc_objects as obj  # noqa: E402
import test_video_qc_sync as syn     # noqa: E402
import core                          # noqa: E402  (on the path through the tests above)


def write_manifest(path, header, rows):
    with open(path, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(header)
        c.writerows(rows)


def run(man, out, *args):
    r = subprocess.run([sys.executable, "-I", TOOL, man, out, *args], check=True, capture_output=True, text=True)
    return r.stdout


def counts(stdout, check):
    """(reused, to check) as the run printed them."""
    m = re.search(rf"^{check}: (\d+) reused, (\d+) to check$", stdout, re.M)
    return int(m.group(1)), int(m.group(2))


def tables(out):
    """Every csv in out_dir, as bytes."""
    res = {}
    for f in sorted(os.listdir(out)):
        if f.endswith(".csv"):
            with open(os.path.join(out, f), "rb") as fh:
                res[f] = fh.read()
    return res


def by_id(path):
    return {r["id"]: r for r in core.read_table(path)}


def stills(out, check):
    d = os.path.join(out, "stills", check)
    return set(os.listdir(d)) if os.path.isdir(d) else set()


def plant(track, rows, k):
    """Rows of a tracking csv get the body point (last landmark) of row + k: poses read out of step."""
    with open(track, newline="") as fh:
        lines = list(csv.reader(fh))
    body = {i: lines[3 + i + k][-3:-1] for i in rows}
    for i, xy in body.items():
        lines[3 + i][-3:-1] = xy
    with open(track, "w", newline="") as fh:
        csv.writer(fh).writerows(lines)


def test_reuse_change_fresh_only_adopt():
    with tempfile.TemporaryDirectory() as d:
        recs = {"clean": syn.make_recording(d, "clean"),
                "bad": syn.make_recording(d, "bad", steps={40: 5, 100: 4, 101: 4}),
                "other": syn.make_recording(d, "other", seed=2, steps={150: 4, 151: 4, 152: 4})}
        man, out = os.path.join(d, "manifest.csv"), os.path.join(d, "out")
        write_manifest(man, ["id", "tracking", "video", "arena_cm", "animal", "group"],
                       [[rid, track, video, "40x40", "auto", "a"] for rid, (video, track) in recs.items()])
        args = ("--checks", "frames,objects,sync", "--sync-windows", "all", "--stills")
        so = run(man, out, *args)
        assert all(counts(so, c) == (0, 3) for c in ("frames", "objects", "sync")), so
        first = tables(out)
        assert set(first) == {"video_qc.csv", "video_qc_frames.csv", "video_qc_objects.csv", "video_qc_sync.csv",
                              "video_qc_sync_frames.csv"}, set(first)
        assert all(len(r["inputs"]) == 16 for r in core.read_table(os.path.join(out, "video_qc_sync.csv")))
        assert stills(out, "sync") == {"bad.jpg", "other.jpg"}
        t_other = os.path.getmtime(os.path.join(out, "stills", "sync", "other.jpg"))

        # nothing changed: every row reused, the same tables
        so = run(man, out, *args)
        assert all(counts(so, c) == (3, 0) for c in ("frames", "objects", "sync")), so
        assert tables(out) == first

        # a table of a version without inputs: taken over with --adopt, else checked again
        frames_csv = os.path.join(out, "video_qc_frames.csv")
        def drop_inputs():
            rows = core.read_table(frames_csv)
            core.write_csv([{k: v for k, v in r.items() if k != "inputs"} for r in rows], frames_csv)
        drop_inputs()
        assert counts(run(man, out, "--checks", "frames", "--adopt"), "frames") == (3, 0)
        assert tables(out) == first
        drop_inputs()
        assert counts(run(man, out, "--checks", "frames"), "frames") == (0, 3)
        assert tables(out) == first

        # clean gets rows out of step, bad is tracked again without them: only these two are checked
        plant(recs["clean"][1], [60, 61, 62], 5)
        syn.make_recording(d, "bad")
        so = run(man, out, *args)
        assert all(counts(so, c) == (1, 2) for c in ("frames", "objects", "sync")), so
        sync_new = by_id(os.path.join(out, "video_qc_sync.csv"))
        sync_was = {r["id"]: r for r in csv.DictReader(first["video_qc_sync.csv"].decode().splitlines())}
        assert sync_new["other"] == sync_was["other"]
        assert sync_new["clean"]["out_of_step_frames"] == "60-62" and sync_new["clean"]["flags"] == "out_of_step"
        assert sync_new["bad"]["flags"] == "" and sync_new["bad"]["out_of_step"] == "0", sync_new["bad"]
        frames_now = core.read_table(os.path.join(out, "video_qc_sync_frames.csv"))
        frames_was = list(csv.DictReader(first["video_qc_sync_frames.csv"].decode().splitlines()))
        assert [r for r in frames_now if r["id"] == "other"] == [r for r in frames_was if r["id"] == "other"]
        assert sorted(int(r["frame"]) for r in frames_now if r["id"] == "clean") == [60, 61, 62]
        assert not any(r["id"] == "bad" for r in frames_now)
        assert stills(out, "sync") == {"clean.jpg", "other.jpg"}                  # bad's went with its flag
        assert os.path.getmtime(os.path.join(out, "stills", "sync", "other.jpg")) == t_other
        changed = tables(out)

        # --fresh checks all again, with the same results; --only keeps the rows of the others
        so = run(man, out, *args, "--fresh")
        assert all(counts(so, c) == (0, 3) for c in ("frames", "objects", "sync")), so
        assert tables(out) == changed
        so = run(man, out, *args, "--fresh", "--only", "^other$")
        assert all(counts(so, c) == (0, 1) for c in ("frames", "objects", "sync")), so
        assert tables(out) == changed and stills(out, "sync") == {"clean.jpg", "other.jpg"}


def test_killed_run_is_continued():
    with tempfile.TemporaryDirectory() as d:
        n = 80
        video = fra.make_clip(d, "clip", 60)
        rows = [[f"r{k:02d}", fra.make_tracking(d, f"r{k:02d}", 60 if k % 7 else 59), video, 2, "a"] for k in range(n)]
        man = os.path.join(d, "manifest.csv")
        write_manifest(man, ["id", "tracking", "video", "expected_s", "group"], rows)
        out = os.path.join(d, "out")
        table = os.path.join(out, "video_qc_frames.csv")
        p = subprocess.Popen([sys.executable, "-I", TOOL, man, out, "--checks", "frames"],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        t0 = time.time()
        while not os.path.exists(table) and p.poll() is None and time.time() - t0 < 120:
            time.sleep(0.02)
        p.kill()
        p.wait()
        part = core.read_table(table)
        assert 20 <= len(part) < n, f"{len(part)} rows: the run ended before it was killed"
        assert not os.path.exists(os.path.join(out, "video_qc.csv"))
        so = run(man, out, "--checks", "frames")
        assert counts(so, "frames") == (len(part), n - len(part)), so
        fresh = os.path.join(d, "fresh")
        run(man, fresh, "--checks", "frames")
        assert tables(out) == tables(fresh)
        flagged = [r["id"] for r in core.read_table(table) if r["flags"]]
        assert flagged == [f"r{k:02d}" for k in range(0, n, 7)], flagged      # 59 rows: rows_differ


def test_group_comparison_after_a_partial_recheck():
    with tempfile.TemporaryDirectory() as d:
        video = os.path.join(d, "arena.mp4")
        obj.make_video(video)
        shifts = {"g1": (0, 0), "g2": (3, 0), "g3": (0, 3), "g4": (12, 9)}       # px; g4 is 1.5 cm off

        def tracking(rid, dx, dy):
            path = os.path.join(d, rid + ".csv")
            obj.write_tracking(path, obj.base_points(objL=obj.fixed(obj.DISK, dx, dy), objR=obj.fixed(obj.BLOCK, -dx, dy),
                                                     flat=obj.fixed(obj.FLAT, dx, dy)))
            return path
        rows = [[f"{g}_{k}", tracking(f"{g}_{k}", *s), video, g, "32x32", "objL objR flat lid"]
                for g, s in shifts.items() for k in range(2)]
        man = os.path.join(d, "manifest.csv")
        write_manifest(man, ["id", "tracking", "video", "group", "arena_cm", "objects"], rows)
        out, fresh = os.path.join(d, "out"), os.path.join(d, "fresh")
        run(man, out, "--checks", "objects", "--stills")
        before = by_id(os.path.join(out, "video_qc_objects.csv"))
        assert {k for k, r in before.items() if r["flags"]} == {"g4_0", "g4_1"}
        assert stills(out, "objects") == {"g4_0.jpg", "g4_1.jpg"}

        tracking("g4_0", 0, 0)                     # g4 now half on its place: no longer off as a group
        tracking("g3_0", 30, 3)                    # g3 half 3 cm off: off as a group
        so = run(man, out, "--checks", "objects", "--stills")
        assert counts(so, "objects") == (6, 2), so
        after = by_id(os.path.join(out, "video_qc_objects.csv"))
        run(man, fresh, "--checks", "objects", "--stills")
        want = by_id(os.path.join(fresh, "video_qc_objects.csv"))
        assert after == want, [(k, after[k], want[k]) for k in want if after[k] != want[k]]
        # kept rows whose group flag changed: g3_1 gained it (its still was made), g4_1 lost it (its still went)
        assert "group_offset" in after["g3_1"]["flags"] and after["g4_1"]["flags"] == "", (after["g3_1"], after["g4_1"])
        assert "1 flagged recordings without a still checked again" in so, so
        assert stills(out, "objects") == stills(fresh, "objects") and "g3_1.jpg" in stills(out, "objects")


if __name__ == "__main__":
    test_reuse_change_fresh_only_adopt(); test_killed_run_is_continued(); test_group_comparison_after_a_partial_recheck()
    print("ok")
