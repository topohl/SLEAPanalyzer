"""tools/video_qc clip check: synthetic raw recordings with a hand that releases the animal and comes back, clips cut
from them well and badly, found by raw_video or by searching a folder with decoys and camera chapters; a shadow
passing over the arena must not count as a hand. Then harder cases: a hand reaching in mid-trial, a second animal
inside the arena, an older cut of the same trial next to the camera recording, a faint animal against an empty-arena
decoy, a light animal on a dark floor, a 30 fps clip of a 60 fps raw, a raw with B-frames and key frames 10 s apart,
and a missing file."""
import csv, os, shutil, subprocess, sys, tempfile
import numpy as np
import cv2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "video_qc")
W, H, FPS = 320, 180, 30
X0, Y0, CW, CH = 84, 10, 160, 160          # where the clips are cut from the raw
ARENA = (100, 26, 228, 154)                # arena floor in raw px
RELEASE, EXPECTED_S = 45, 12               # the hand leaves at raw frame 45; the trial lasts 12 s
END = RELEASE + EXPECTED_S * FPS           # first frame after the trial: 405


def scene(floor=200):
    img = np.full((H, W), 70, np.uint8)
    rng = np.random.default_rng(0)
    for _ in range(30):                     # texture around the arena, so the crop has one place
        x, y = int(rng.integers(0, W - 20)), int(rng.integers(0, H - 20))
        img[y:y + int(rng.integers(5, 20)), x:x + int(rng.integers(5, 20))] = int(rng.integers(30, 130))
    cv2.rectangle(img, ARENA[:2], ARENA[2:], floor, -1)
    cv2.rectangle(img, ARENA[:2], ARENA[2:], 140, 2)
    return img


def walk(n, seed):
    """Animal path: held at the centre until the release, then a smooth random walk on the floor."""
    rng = np.random.default_rng(seed)
    p, v, out = np.array([164.0, 90.0]), np.zeros(2), []
    lo, hi = np.array([ARENA[0] + 10, ARENA[1] + 10]), np.array([ARENA[2] - 10, ARENA[3] - 10])
    for i in range(n):
        if i >= RELEASE:
            v = 0.9 * v + rng.normal(0, 0.9, 2)
            v *= min(1.0, 4 / max(np.hypot(*v), 1e-9))
            p = p + v
            bounce = (p < lo) | (p > hi)
            v[bounce] *= -1
            p = np.clip(p, lo, hi)
        out.append(p.copy())
    return np.array(out)


def frames(n, seed, hand_back=None, shadow=None, mid=None, second=None, floor=200, animal=20, hand=35):
    """Raw frames: the hand (arm from the top edge) until RELEASE, in frames mid[0] .. mid[1]-1 and again from
    hand_back; a shadow (the scene darkened to 60 %) over the top left of the arena and its surround in frames
    shadow[0] .. shadow[1]-1; a second, larger animal inside the arena (as in a cup) in frames second[0] ..
    second[1]-1. Grey levels of the floor, the animal and the hand can be set."""
    bg, path = scene(floor), walk(n, seed)
    for i in range(n):
        img = bg.copy()
        x, y = path[i]
        cv2.ellipse(img, (int(x), int(y)), (8, 4), 30, 0, 360, animal, -1)
        if second and second[0] <= i < second[1]:
            cv2.ellipse(img, (140, 62), (20, 10), 0, 0, 360, animal, -1)
        if i < RELEASE or (hand_back is not None and i >= hand_back) or (mid and mid[0] <= i < mid[1]):
            cv2.rectangle(img, (150, 0), (178, 80), hand, -1)
            cv2.ellipse(img, (164, 90), (30, 22), 0, 0, 360, hand, -1)
        if shadow and shadow[0] <= i < shadow[1]:
            img[0:110, 60:200] = (img[0:110, 60:200] * 0.6).astype(np.uint8)
        yield img, path[i]


def write(path, imgs, size, fps=FPS):
    w = cv2.VideoWriter(path, cv2.VideoWriter_fourcc(*"mp4v"), fps, size)
    for im in imgs:
        w.write(cv2.cvtColor(im, cv2.COLOR_GRAY2BGR))
    w.release()


def make_raw(path, n, seed, hand_back=None, shadow=None, fps=FPS, **kw):
    data = list(frames(n, seed, hand_back, shadow, **kw))
    write(path, [im for im, _ in data], (W, H), fps)
    return data


def make_clip(d, name, data, a, b, scale=1.0, step=1, fps=FPS, folder="cropped"):
    """Clip of raw frames a..b (inclusive; every step-th) and its tracking csv (corners and body centre in clip px)."""
    size = (int(CW * scale), int(CH * scale))
    data = data[a:b + 1:step]
    crops = [cv2.resize(im[Y0:Y0 + CH, X0:X0 + CW], size, interpolation=cv2.INTER_AREA) for im, _ in data]
    video = os.path.join(d, folder, name + ".mp4")
    os.makedirs(os.path.dirname(video), exist_ok=True)
    write(video, crops, size, fps)
    corners = [(ARENA[0], ARENA[1]), (ARENA[2], ARENA[1]), (ARENA[2], ARENA[3]), (ARENA[0], ARENA[3])]
    track = os.path.join(d, name + ".csv")
    with open(track, "w", newline="") as fh:
        c = csv.writer(fh)
        parts = ["tl", "tr", "br", "bl", "bodycentre"]
        c.writerow(["scorer"] * 16)
        c.writerow(["bodyparts"] + [p for p in parts for _ in range(3)])
        c.writerow(["coords"] + ["x", "y", "likelihood"] * 5)
        for k, (_, (x, y)) in enumerate(data):
            pts = corners + [(x, y)]
            c.writerow([k] + [v for px, py in pts for v in ((px - X0) * scale, (py - Y0) * scale, 1)])
    return video, track


def run_tool(d, rows):
    manifest = os.path.join(d, "manifest.csv")
    keys = ["id", "tracking", "video", "raw_video", "raw_dir", "expected_s", "clip_scale"]
    with open(manifest, "w", newline="") as fh:
        c = csv.DictWriter(fh, fieldnames=keys)
        c.writeheader()
        for r in rows:
            c.writerow({k: r.get(k, "") for k in keys})
    out = os.path.join(d, "out")
    subprocess.run([sys.executable, "-I", TOOL, manifest, out, "--checks", "clip", "--stills"], check=True)
    with open(os.path.join(out, "video_qc_clip.csv"), newline="", encoding="utf-8") as fh:
        return {r["id"]: r for r in csv.DictReader(fh)}, out


def flags(row):
    return sorted(f for f in row["flags"].split(";") if f)


def test_clip_check():
    with tempfile.TemporaryDirectory() as d:
        raws = os.path.join(d, "raw")
        os.makedirs(raws)
        clean = make_raw(os.path.join(raws, "clean.mp4"), 450, 1, hand_back=420)
        alone = make_raw(os.path.join(raws, "alone.mp4"), 450, 2)                 # nobody comes back
        handed = make_raw(os.path.join(raws, "handed.mp4"), 450, 3, hand_back=255)
        short = make_raw(os.path.join(raws, "short.mp4"), 300, 4)                # the camera stops at 10 s
        shaded = make_raw(os.path.join(raws, "shaded.mp4"), 450, 5, hand_back=420, shadow=(150, 180))
        # a camera folder: the clean recording in two chapters, two other trials of the same view
        cam = os.path.join(d, "session", "C1")
        os.makedirs(cam)
        write(os.path.join(cam, "GX010042.mp4"), [im for im, _ in clean[:240]], (W, H))
        write(os.path.join(cam, "GX020042.mp4"), [im for im, _ in clean[240:]], (W, H))
        for k, seed in ((43, 7), (44, 8)):
            make_raw(os.path.join(cam, f"GX01{k:04d}.mp4"), 450, seed, hand_back=420)
        decoys = os.path.join(d, "decoys")
        os.makedirs(decoys)
        for k in (43, 44):
            shutil.copy(os.path.join(cam, f"GX01{k:04d}.mp4"), decoys)

        v_ok, t_ok = make_clip(d, "ok", clean, RELEASE, END - 1)
        v_cut, t_cut = make_clip(d, "cut", alone, RELEASE, 254)
        v_hand, t_hand = make_clip(d, "hand", handed, RELEASE, 254)
        v_short, t_short = make_clip(d, "short", short, RELEASE, 299)
        v_late, t_late = make_clip(d, "late", clean, RELEASE + 120, END - 1)
        v_early, t_early = make_clip(d, "early", clean, 0, END - 1)
        v_half, t_half = make_clip(d, "half", clean, RELEASE, END - 1, scale=0.5)
        v_shade, t_shade = make_clip(d, "shade", shaded, RELEASE, END - 1)
        rawp = lambda name: os.path.join(raws, name)
        rows = [
            dict(id="ok", tracking=t_ok, video=v_ok, raw_video=rawp("clean.mp4"), expected_s=EXPECTED_S),
            dict(id="cut", tracking=t_cut, video=v_cut, raw_video=rawp("alone.mp4"), expected_s=EXPECTED_S),
            dict(id="hand", tracking=t_hand, video=v_hand, raw_video=rawp("handed.mp4"), expected_s=EXPECTED_S),
            dict(id="short", tracking=t_short, video=v_short, raw_video=rawp("short.mp4"), expected_s=EXPECTED_S),
            dict(id="late", tracking=t_late, video=v_late, raw_video=rawp("clean.mp4"), expected_s=EXPECTED_S),
            dict(id="early", tracking=t_early, video=v_early, raw_video=rawp("clean.mp4"), expected_s=EXPECTED_S),
            dict(id="search", tracking=t_ok, video=v_ok, raw_dir=os.path.join(d, "session"), expected_s=EXPECTED_S),
            dict(id="missing", tracking=t_ok, video=v_ok, raw_dir=decoys, expected_s=EXPECTED_S),
            dict(id="noexp", tracking=t_ok, video=v_ok, raw_video=rawp("clean.mp4")),
            dict(id="half", tracking=t_half, video=v_half, raw_video=rawp("clean.mp4"), expected_s=EXPECTED_S,
                 clip_scale=0.5),
            dict(id="shade", tracking=t_shade, video=v_shade, raw_video=rawp("shaded.mp4"), expected_s=EXPECTED_S),
        ]
        res, out = run_tool(d, rows)

        assert flags(res["ok"]) == [], res["ok"]
        assert abs(int(res["ok"]["raw_start"]) - RELEASE) <= 1 and abs(int(res["ok"]["raw_end"]) - (END - 1)) <= 1
        assert abs(int(res["ok"]["crop_x_px"]) - X0) <= 1 and abs(int(res["ok"]["crop_y_px"]) - Y0) <= 1
        assert abs(int(res["ok"]["release"]) - RELEASE) <= 2 and abs(int(res["ok"]["trial_end"]) - 420) <= 2
        assert res["ok"]["alone_to_expected"] == "yes"

        assert flags(res["cut"]) == ["cut_short"], res["cut"]
        assert abs(float(res["cut"]["recoverable_s"]) - 5) <= 0.2 and res["cut"]["alone_to_expected"] == "yes"
        assert res["cut"]["trial_end_by"] == "raw_end"

        assert flags(res["hand"]) == ["ended_early"], res["hand"]
        assert abs(float(res["hand"]["next_intrusion_after_clip_s"])) <= 0.2
        assert flags(res["short"]) == ["ended_early"], res["short"]
        assert res["short"]["trial_end_by"] == "raw_end" and res["short"]["alone_to_expected"] == "no"

        assert flags(res["late"]) == ["starts_late"], res["late"]
        assert abs(float(res["late"]["start_minus_release_s"]) - 4) <= 0.2
        assert flags(res["early"]) == ["starts_early"], res["early"]
        assert abs(float(res["early"]["start_minus_release_s"]) + 1.5) <= 0.2

        assert flags(res["search"]) == [], res["search"]
        assert res["search"]["found_by"] == "search" and "GX020042" in res["search"]["raw"]
        assert abs(int(res["search"]["raw_start"]) - RELEASE) <= 1 and res["search"]["end_at"].startswith("GX020042")
        assert flags(res["missing"]) == ["raw_not_found"], res["missing"]
        assert float(res["missing"]["align_score"]) > 0.4

        assert flags(res["noexp"]) == ["no_expected"], res["noexp"]
        assert flags(res["half"]) == [], res["half"]
        assert abs(int(res["half"]["raw_start"]) - RELEASE) <= 1 and abs(int(res["half"]["crop_w_px"]) - CW) <= 2
        assert flags(res["shade"]) == [], res["shade"]          # a shadow over the arena is not a hand
        assert abs(int(res["shade"]["trial_end"]) - 420) <= 2

        stills = os.listdir(os.path.join(out, "stills", "clip"))
        assert {"cut.jpg", "hand.jpg", "late.jpg", "early.jpg"} <= set(stills) and "ok.jpg" not in stills


def test_clip_check_harder():
    sys.path.insert(0, TOOL)
    import core
    ffmpeg = core.find_tool("ffmpeg")
    with tempfile.TemporaryDirectory() as d:
        raws = os.path.join(d, "raw")
        os.makedirs(raws)
        rawp = lambda name: os.path.join(raws, name)
        rows = []
        # a hand reaches in for 1 s in mid-trial; the animal is there before and after
        mid = make_raw(rawp("mid.mp4"), 450, 11, hand_back=420, mid=(200, 230))
        v, t = make_clip(d, "midhand", mid, RELEASE, END - 1)
        rows.append(dict(id="midhand", tracking=t, video=v, raw_video=rawp("mid.mp4"), expected_s=EXPECTED_S))
        # a second animal shows up inside the arena (never reaching out of it) in mid-trial: not a hand
        cup = make_raw(rawp("cup.mp4"), 450, 19, hand_back=420, second=(150, 300))
        v, t = make_clip(d, "cup", cup, RELEASE, END - 1)
        rows.append(dict(id="cup", tracking=t, video=v, raw_video=rawp("cup.mp4"), expected_s=EXPECTED_S))
        # an older cut of the same trial (from raw 0, the hand in it) in a folder next to the camera folder; its name
        # is in the clip's name, so it would be tried first
        sess = os.path.join(d, "session")
        os.makedirs(os.path.join(sess, "C1"))
        dup = make_raw(os.path.join(sess, "C1", "GX010077.mp4"), 450, 12, hand_back=420)
        make_clip(sess, "M1K1", dup, 0, END - 1, folder="first_cut")
        v, t = make_clip(sess, "M1K1_recut", dup, RELEASE, END - 1, folder="clips")
        rows.append(dict(id="recut", tracking=t, video=v, raw_dir=sess, expected_s=EXPECTED_S))
        # the same, but the camera recording is gone: the older cut is taken, and flagged
        gone = os.path.join(d, "gone")
        for sub in ("clips", "first_cut"):
            shutil.copytree(os.path.join(sess, sub), os.path.join(gone, sub))
        os.makedirs(os.path.join(gone, "C1"))
        make_raw(os.path.join(gone, "C1", "GX010078.mp4"), 450, 20, hand_back=420)              # another trial
        v = os.path.join(gone, "clips", "M1K1_recut.mp4")
        rows.append(dict(id="recut_gone", tracking=t, video=v, raw_dir=gone, expected_s=EXPECTED_S))
        # the older cut in the clip's own folder: the other clips there are not candidates
        beside = os.path.join(d, "beside")
        shutil.copytree(os.path.join(gone, "C1"), os.path.join(beside, "C1"))
        shutil.copytree(os.path.join(gone, "clips"), os.path.join(beside, "clips"))
        shutil.copy(os.path.join(gone, "first_cut", "M1K1.mp4"), os.path.join(beside, "clips"))
        rows.append(dict(id="recut_beside", tracking=t, video=os.path.join(beside, "clips", "M1K1_recut.mp4"),
                         raw_dir=beside, expected_s=EXPECTED_S))
        # a faint animal (28 grey levels on the floor): found when given; an empty-arena recording and another trial
        # of the same view must not pass for its raw
        faint = make_raw(rawp("faint.mp4"), 450, 13, hand_back=420, animal=172)
        decoys = os.path.join(d, "decoys")
        os.makedirs(decoys)
        make_raw(os.path.join(decoys, "GX010001.mp4"), 450, 14, hand_back=420, animal=200)   # nothing to see
        make_raw(os.path.join(decoys, "GX010002.mp4"), 450, 15, hand_back=420, animal=172)
        v, t = make_clip(d, "faint", faint, RELEASE, END - 1)
        rows.append(dict(id="faint", tracking=t, video=v, raw_video=rawp("faint.mp4"), expected_s=EXPECTED_S))
        rows.append(dict(id="faint_decoys", tracking=t, video=v, raw_dir=decoys, expected_s=EXPECTED_S))
        # a light animal and hand on a dark floor
        light = make_raw(rawp("light.mp4"), 450, 16, hand_back=420, floor=40, animal=230, hand=210)
        v, t = make_clip(d, "light", light, RELEASE, END - 1)
        rows.append(dict(id="light", tracking=t, video=v, raw_video=rawp("light.mp4"), expected_s=EXPECTED_S))
        # a 60 fps raw cut into a 30 fps clip (every second frame): whole, and 5 s short
        r60 = make_raw(rawp("r60.mp4"), 900, 17, hand_back=RELEASE + 2 * EXPECTED_S * FPS + 15, fps=2 * FPS)
        end60 = RELEASE + 2 * EXPECTED_S * FPS
        v, t = make_clip(d, "fps", r60, RELEASE, end60 - 1, step=2)
        rows.append(dict(id="fps", tracking=t, video=v, raw_video=rawp("r60.mp4"), expected_s=EXPECTED_S))
        v, t = make_clip(d, "fps_cut", r60, RELEASE, end60 - 1 - 10 * FPS, step=2)
        rows.append(dict(id="fps_cut", tracking=t, video=v, raw_video=rawp("r60.mp4"), expected_s=EXPECTED_S))
        # B-frames and key frames 10 s apart; a 60 s trial cut 10 s short
        long_ = make_raw(rawp("long.mp4"), 2400, 18, hand_back=RELEASE + 60 * FPS + 120)
        subprocess.run([ffmpeg, "-v", "error", "-y", "-i", rawp("long.mp4"), "-c:v", "mpeg4", "-g", "300", "-bf", "2",
                        "-q:v", "3", rawp("bframes.mp4")], check=True)
        v, t = make_clip(d, "bframes", long_, RELEASE, RELEASE + 50 * FPS - 1)
        rows.append(dict(id="bframes", tracking=t, video=v, raw_video=rawp("bframes.mp4"), expected_s=60))
        rows.append(dict(id="missing", tracking=t, video=v, raw_video=rawp("nothing.mp4"), expected_s=60))
        res, out = run_tool(d, rows)

        assert flags(res["midhand"]) == ["intrusion_in_clip"], res["midhand"]
        assert abs(int(res["midhand"]["release"]) - RELEASE) <= 2 and abs(int(res["midhand"]["trial_end"]) - 420) <= 2
        first = float(res["midhand"]["intrusion_in_clip_s"].split("-")[0])
        assert abs(first - (200 - RELEASE) / FPS) <= 1.0, res["midhand"]["intrusion_in_clip_s"]
        assert flags(res["cup"]) == [] and abs(int(res["cup"]["trial_end"]) - 420) <= 2, res["cup"]
        assert flags(res["recut"]) == [] and "GX010077" in res["recut"]["raw"], res["recut"]
        gone = res["recut_gone"]
        assert flags(gone) == ["raw_clip_sized"] and gone["raw"].endswith("M1K1.mp4"), gone
        assert flags(res["recut_beside"]) == ["raw_not_found"] and "M1K1" not in res["recut_beside"]["tried"], res["recut_beside"]
        assert flags(res["faint"]) == [], res["faint"]
        assert flags(res["faint_decoys"]) == ["raw_not_found"], res["faint_decoys"]
        assert float(res["faint_decoys"]["align_score"]) > 0.4
        assert flags(res["light"]) == [], res["light"]
        assert abs(int(res["light"]["trial_end"]) - 420) <= 2
        assert flags(res["fps"]) == [] and int(res["fps"]["raw_span_minus_clip_frames"]) == 0, res["fps"]
        assert abs(float(res["fps"]["end_minus_expected_s"])) <= 0.05
        assert flags(res["fps_cut"]) == ["cut_short"] and abs(float(res["fps_cut"]["recoverable_s"]) - 5) <= 0.2
        assert flags(res["bframes"]) == ["cut_short"], res["bframes"]
        assert abs(float(res["bframes"]["recoverable_s"]) - 10) <= 0.2
        assert abs(int(res["bframes"]["raw_start"]) - RELEASE) <= 1
        assert flags(res["missing"]) == ["raw_not_found"] and res["missing"]["error"].startswith("cannot read")
        assert "midhand.jpg" in os.listdir(os.path.join(out, "stills", "clip"))


if __name__ == "__main__":
    test_clip_check()
    test_clip_check_harder()
    print("ok")
