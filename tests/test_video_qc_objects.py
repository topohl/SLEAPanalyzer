"""tools/video_qc objects check: points on, off and between synthetic objects, one group placed differently, groups
compared within object sets, and other framings (wide frame, light animal, no outline, lines on the floor)."""
import csv, os, subprocess, sys, tempfile
from types import SimpleNamespace
import numpy as np
import cv2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "video_qc")
sys.path.insert(0, TOOL)
import objects  # noqa: E402

# 400 x 400 px video, arena 32 x 32 cm from (40, 40) to (360, 360): 10 px per cm
DISK = (130, 250)          # blue disk, radius 25 px
BLOCK = (270, 250)         # orange block, 50 x 36 px
FLAT = (200, 120)          # floor-coloured square with a thin dark outline: found only by its edges
LID = (335, 90)            # a lid like FLAT reaching past the arena edge: closed by the outline
CORNERS = {"tl": (40, 40), "tr": (360, 40), "br": (360, 360), "bl": (40, 360)}


def make_video(path, n=160, fps=30):
    w = cv2.VideoWriter(path, cv2.VideoWriter_fourcc(*"mp4v"), fps, (400, 400))
    for i in range(n):
        img = np.full((400, 400, 3), 60, np.uint8)
        cv2.rectangle(img, (40, 40), (360, 360), (200, 200, 205), -1)
        cv2.circle(img, DISK, 25, (200, 60, 30), -1)
        cv2.rectangle(img, (BLOCK[0] - 25, BLOCK[1] - 18), (BLOCK[0] + 25, BLOCK[1] + 18), (30, 120, 240), -1)
        cv2.rectangle(img, (FLAT[0] - 30, FLAT[1] - 30), (FLAT[0] + 30, FLAT[1] + 30), (150, 150, 150), 1)
        cv2.rectangle(img, (320, 60), (380, 120), (150, 150, 150), 1)
        a = 2 * np.pi * i / 50                                  # the animal runs circles in the open middle
        cv2.ellipse(img, (int(200 + 70 * np.cos(a)), int(220 + 30 * np.sin(a))), (12, 8), 0, 0, 360, (40, 40, 40), -1)
        w.write(img)
    w.release()


def write_tracking(path, points, n=160):
    """points: landmark -> function(frame) -> (x, y)."""
    names = list(points)
    with open(path, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(["scorer"] * (1 + 3 * len(names)))
        c.writerow(["bodyparts"] + [p for p in names for _ in range(3)])
        c.writerow(["coords"] + ["x", "y", "likelihood"] * len(names))
        for i in range(n):
            row = [i]
            for p in names:
                x, y = points[p](i)
                row += [round(x, 2), round(y, 2), 1]
            c.writerow(row)


def fixed(p, dx=0, dy=0):
    return lambda i: (p[0] + dx, p[1] + dy)


def share(p_main, p_other, frac_other):
    """On p_main, but on p_other in the given share of frames."""
    return lambda i: p_other if (i % 20) < 20 * frac_other else p_main


def base_points(**over):
    pts = {k: fixed(v) for k, v in CORNERS.items()}
    pts.update(objL=fixed(DISK), objR=fixed(BLOCK), flat=fixed(FLAT), lid=fixed(LID), bodycentre=lambda i: (200, 220))
    pts.update(over)
    return {k: v for k, v in pts.items() if v is not None}


def run_tool(d, recs, extra=()):
    video = os.path.join(d, "arena.mp4")
    if not os.path.exists(video):
        make_video(video)
    man = os.path.join(d, "manifest.csv")
    with open(man, "w", newline="") as fh:
        c = csv.writer(fh)
        c.writerow(["id", "tracking", "video", "group", "arena_cm", "objects"])
        for rid, pts, group, arena in recs:
            track = os.path.join(d, rid + ".csv")
            write_tracking(track, pts)
            c.writerow([rid, track, video, group, arena, "objL objR flat lid"])
    out = os.path.join(d, "out")
    subprocess.run([sys.executable, "-I", TOOL, man, out, "--checks", "objects", *extra], check=True,
                   stdout=subprocess.DEVNULL)
    with open(os.path.join(out, "video_qc_objects.csv"), newline="") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}, out


def test_points_on_off_jumping_and_merged():
    with tempfile.TemporaryDirectory() as d:
        rows, out = run_tool(d, [
            ("clean", base_points(), "", "32x32"),
            ("off", base_points(objR=fixed(BLOCK, 0, 70)), "", "32x32"),            # 7 cm below the block
            ("jumps", base_points(objL=share(DISK, (200, 320), 0.3)), "", "32x32"),  # on the floor 30% of frames
            ("merged", base_points(objL=share(DISK, BLOCK, 0.75)), "", "32x32"),     # on the other object 75%
            ("missing", base_points(flat=None), "", "32x32"),
            ("noscale", base_points(), "", ""),
        ], extra=("--stills",))
        flags = {k: set(filter(None, r["flags"].split(";"))) for k, r in rows.items()}
        assert flags["clean"] == set(), flags["clean"]
        c = rows["clean"]
        for o in ("objL", "objR", "flat", "lid"):
            assert c[f"{o}_on_region"] == "True" and float(c[f"{o}_to_region_cm"]) == 0, o
            assert float(c[f"{o}_near_median"]) == 1
        for o in ("objL", "objR", "flat"):
            assert abs(float(c[f"{o}_dx_cm"])) < 0.3 and abs(float(c[f"{o}_dy_cm"])) < 0.3, o
        assert float(c["lid_region_area_cm2"]) > 12                # the lid's floor inside the outline, not only its edges
        assert 18 < float(c["objL_region_area_cm2"]) < 25       # disk radius 2.5 cm: 19.6 cm2, + the 0.2 cm rim
        assert len({c["objL_region"], c["objR_region"], c["flat_region"], c["lid_region"]}) == 4
        assert abs(float(c["objL_objR_sep_cm"]) - 14) < 0.1 and abs(float(c["objL_objR_region_sep_cm"]) - 14) < 0.3
        assert flags["off"] == {"off_object"}, flags["off"]
        assert rows["off"]["objR_on_region"] == "False" and 4.5 < float(rows["off"]["objR_to_region_cm"]) < 6
        assert flags["jumps"] == {"object_jumps"}, flags["jumps"]
        assert abs(float(rows["jumps"]["objL_near_median"]) - 0.7) < 0.01
        assert flags["merged"] == {"objects_merged", "object_jumps"}, flags["merged"]
        assert rows["merged"]["objL_region"] == rows["merged"]["objR_region"]
        assert flags["missing"] == {"missing_landmark"}, flags["missing"]
        assert flags["noscale"] == set() and "objL_dx_px" in rows["noscale"] and rows["noscale"]["objL_dx_cm"] == ""
        stills = set(os.listdir(os.path.join(out, "stills", "objects")))
        assert stills == {"off.jpg", "jumps.jpg", "merged.jpg", "missing.jpg"}, stills
        img = cv2.imread(os.path.join(out, "stills", "objects", "off.jpg"))
        assert img.shape[0] == 430 and img.shape[1] > 400


def test_one_group_placed_differently():
    with tempfile.TemporaryDirectory() as d:
        recs = []
        shifts = {"g1": (0, 0), "g2": (3, 0), "g3": (0, 3), "g4": (12, 9)}       # px; g4 is 1.5 cm off
        for g, (dx, dy) in shifts.items():
            for k in range(2):
                recs.append((f"{g}_{k}", base_points(objL=fixed(DISK, dx, dy), objR=fixed(BLOCK, -dx, dy),
                                                     flat=fixed(FLAT, dx, dy)), g, "32x32"))
        # a recording whose objL sits on the block must not move g1's median
        recs.append(("g1_merged", base_points(objL=fixed(BLOCK, 20, 0)), "g1", "32x32"))
        rows, _ = run_tool(d, recs)
        for rid, r in rows.items():
            f = set(filter(None, r["flags"].split(";")))
            want = {"group_offset"} if rid.startswith("g4") else set()
            if rid == "g1_merged":
                want = {"objects_merged"}
            assert f == want, (rid, f)
        assert float(rows["g4_0"]["objL_group_shift_cm"]) > 1.2
        assert float(rows["g1_0"]["objL_group_shift_cm"]) < 0.6


def wide_arena(i, grid=False):
    """640 x 360 frame, arena 32 x 32 cm (5 px per cm) in its upper right, dark floor, a light animal."""
    img = np.full((360, 640, 3), 200, np.uint8)
    cv2.rectangle(img, (420, 30), (580, 190), (35, 35, 40), -1)
    if grid:                                                     # painted lines every 4 cm
        for k in range(0, 161, 20):
            cv2.line(img, (420 + k, 30), (420 + k, 190), (90, 90, 90), 1)
            cv2.line(img, (420, 30 + k), (580, 30 + k), (90, 90, 90), 1)
    cv2.circle(img, (460, 120), 12, (40, 200, 230), -1)
    cv2.rectangle(img, (530, 100), (552, 140), (230, 230, 230), -1)
    a = 2 * np.pi * i / 40
    cv2.ellipse(img, (int(500 + 25 * np.cos(a)), int(70 + 10 * np.sin(a))), (7, 4), 0, 0, 360, (235, 235, 235), -1)
    return img


def test_other_framings():
    """A light animal on a dark floor in an arena that fills a corner of a wide 25 fps frame passes; tracking without
    the outline landmarks and a floor with lines on it are said so, not reported as merged objects."""
    with tempfile.TemporaryDirectory() as d:
        videos = {}
        for name, grid in (("plain", False), ("grid", True)):
            videos[name] = os.path.join(d, name + ".mp4")
            w = cv2.VideoWriter(videos[name], cv2.VideoWriter_fourcc(*"mp4v"), 25, (640, 360))
            for i in range(120):
                w.write(wide_arena(i, grid))
            w.release()
        corners = {"tl": fixed((420, 30)), "tr": fixed((580, 30)), "br": fixed((580, 190)), "bl": fixed((420, 190))}
        pts = dict(objL=fixed((460, 120)), objR=fixed((541, 120)), bodycentre=fixed((500, 70)))
        with_corners, without = os.path.join(d, "t.csv"), os.path.join(d, "t_nocorners.csv")
        write_tracking(with_corners, dict(corners, **pts), n=120)
        write_tracking(without, pts, n=120)
        man = os.path.join(d, "manifest.csv")
        with open(man, "w", newline="") as fh:
            c = csv.writer(fh)
            c.writerow(["id", "tracking", "video", "arena_cm", "objects"])
            c.writerow(["light_animal", with_corners, videos["plain"], "32x32", "objL objR"])
            c.writerow(["no_corners", without, videos["plain"], "32x32", "objL objR"])
            c.writerow(["lines", with_corners, videos["grid"], "32x32", "objL objR"])
            c.writerow(["no_video", with_corners, os.path.join(d, "missing.mp4"), "32x32", "objL objR"])
        out = os.path.join(d, "out")
        subprocess.run([sys.executable, "-I", TOOL, man, out, "--checks", "objects", "--stills"], check=True,
                       stdout=subprocess.DEVNULL)
        with open(os.path.join(out, "video_qc_objects.csv"), newline="") as fh:
            rows = {r["id"]: r for r in csv.DictReader(fh)}
        flags = {k: set(filter(None, r["flags"].split(";"))) for k, r in rows.items()}
        assert flags["light_animal"] == set(), flags["light_animal"]
        r = rows["light_animal"]
        assert r["objL_on_region"] == "True" and r["objR_on_region"] == "True" and r["objL_region"] != r["objR_region"]
        assert abs(float(r["objL_dx_cm"])) < 0.3 and float(r["regions_share"]) < 0.1
        assert flags["no_corners"] == {"no_outline"}, flags["no_corners"]
        assert rows["no_corners"]["objL_on_region"] == "" and float(rows["no_corners"]["objL_near_median"]) == 1
        assert flags["lines"] == {"uneven_floor"}, flags["lines"]
        assert flags["no_video"] == {"unreadable"}, flags["no_video"]
        assert os.listdir(os.path.join(out, "stills", "objects")) == ["lines.jpg"]


def test_finish_compares_within_object_sets():
    """Two tests share landmark names but not object shapes; group a ran only test A. Pooled, every group looks
    shifted; compared within each object set, none is."""
    row = lambda rid, g, s, dy: {"id": rid, "group": g, "object_set": s, "objL_dx_cm": 0.0, "objL_dy_cm": dy, "flags": []}
    def rows(with_sets):
        out = [row("a1", "a", "A", 0.0), row("a2", "a", "A", 0.1)]
        for g in "bc":
            out += [row(g + "1", g, "A", 0.0), row(g + "2", g, "B", -2.9), row(g + "3", g, "B", -2.9)]
        return [dict(r, object_set=r["object_set"] if with_sets else "") for r in out]
    pooled, by_set = rows(False), rows(True)
    objects.finish(pooled, SimpleNamespace(stills=False))
    objects.finish(by_set, SimpleNamespace(stills=False))
    assert all("group_offset" in r["flags"] for r in pooled if r["group"] == "a")
    assert all(r["flags"] == [] for r in by_set), [r["flags"] for r in by_set]


def test_finish_needs_two_groups():
    opts = SimpleNamespace(stills=False)
    row = lambda g, dx: {"id": g + str(dx), "group": g, "objL_dx_cm": dx, "objL_dy_cm": 0, "flags": []}
    rows = [row("a", 0.0), row("a", 0.1), row("", 3.0)]                      # one group; no group does not count
    objects.finish(rows, opts)
    assert all(r["flags"] == [] for r in rows) and "objL_group_shift_cm" not in rows[2]
    rows = [row("a", 0.0), row("b", 2.0)]                                     # two groups: both are flagged
    objects.finish(rows, opts)
    assert all(r["flags"] == ["group_offset"] for r in rows)
    rows = [{"id": "x", "group": "a", "objL_dx_px": 0.0, "objL_dy_px": 0, "flags": []},
            {"id": "y", "group": "b", "objL_dx_px": 12.0, "objL_dy_px": 0, "flags": []}]
    objects.finish(rows, SimpleNamespace(stills=False, object_offset_cm=1.0, object_px_per_cm=16.0))
    assert all(r["flags"] == [] for r in rows)                                # 12 px is under 1 cm at 16 px per cm
    # cups placed differently in batch b must not flag batch b's recordings of another test
    rows = [{"id": "s1", "group": "a", "socl_dx_cm": 0.0, "socl_dy_cm": 0, "flags": []},
            {"id": "s2", "group": "b", "socl_dx_cm": 3.0, "socl_dy_cm": 0, "flags": []},
            {"id": "n1", "group": "b", "objL_dx_cm": 0.0, "objL_dy_cm": 0, "flags": []}]
    objects.finish(rows, opts)
    assert rows[1]["flags"] == ["group_offset"] and rows[2]["flags"] == [] and "socl_group_shift_cm" not in rows[2]


if __name__ == "__main__":
    test_finish_needs_two_groups(); test_finish_compares_within_object_sets()
    test_points_on_off_jumping_and_merged(); test_one_group_placed_differently(); test_other_framings()
    print("ok")
