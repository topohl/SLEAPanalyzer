"""frames: does the video match its tracking, and is the trial as long as the protocol?

Per recording: the video's frame count (from the container; decoded when the container has none), its frame rate
and length, the tracking rows, and the length the trial should have (expected_s).

Flags:
  rows_differ   tracking rows != video frames: the tracking does not belong to this clip as it is now, or frames
                were dropped or added between the two
  length        frames / frame rate differs from expected_s by more than --length-tolerance-s
  rate          the video's frame rate is unknown, or its average differs from its nominal rate by more than 0.2%:
                a variable frame rate, under which frame numbers and times do not convert
"""

import math

import core


def run(rec, ctx):
    v, opts = ctx["video"], ctx["opts"]
    frames = v.nb_frames if v.nb_frames else v.count_frames()
    row = {"rows": ctx["rows"], "video_frames": frames, "fps": core.rnd(v.fps, 4),
           "video_s": core.rnd(frames / v.fps if v.fps else float("nan"), 2),
           "width": v.width, "height": v.height}
    flags = []
    if frames != ctx["rows"]:
        flags.append("rows_differ")
        row["rows_minus_frames"] = ctx["rows"] - frames
    if not (v.fps and math.isfinite(v.fps)):
        flags.append("rate")
    elif math.isfinite(v.nominal_fps) and abs(v.fps / v.nominal_fps - 1) > 0.002:
        flags.append("rate")
        row["nominal_fps"] = core.rnd(v.nominal_fps, 4)
    exp = core.number(rec, "expected_s")
    if exp and v.fps:
        dev = frames / v.fps - exp
        row["length_minus_expected_s"] = core.rnd(dev, 2)
        if abs(dev) > opts.length_tolerance_s:
            flags.append("length")
    row["flags"] = flags
    return row
