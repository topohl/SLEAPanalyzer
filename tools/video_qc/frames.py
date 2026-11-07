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

VERSION = 1           # results version: raise it when a change can alter the measurements (see core.VERSION)
RATE_LIMIT = 0.002    # rate: average and nominal frame rate differ by more than this share


def run(rec, ctx):
    v, opts = ctx["video"], ctx["opts"]
    frames = v.nb_frames if v.nb_frames else v.count_frames()
    row = {"rows": ctx["rows"], "video_frames": frames, "fps": core.rnd(v.fps, 4),
           "video_s": core.rnd(frames / v.fps if v.fps else float("nan"), 2),
           "width": v.width, "height": v.height}
    if frames != ctx["rows"]:
        row["rows_minus_frames"] = ctx["rows"] - frames
    if v.fps and math.isfinite(v.fps) and math.isfinite(v.nominal_fps):
        row["nominal_fps"] = core.rnd(v.nominal_fps, 4)
    exp = core.number(rec, "expected_s")
    if exp and v.fps:
        row["length_minus_expected_s"] = core.rnd(frames / v.fps - exp, 2)
    row["flags"] = flags(row, rec, opts)
    return row


def flags(row, rec, opts):
    """The flags of a row, from its columns (a row of this run or one kept from an earlier run)."""
    out = []
    rows, frames = core.num(row, "rows"), core.num(row, "video_frames")
    if rows is not None and frames is not None and rows != frames:
        out.append("rows_differ")
    fps, nominal = core.num(row, "fps"), core.num(row, "nominal_fps")
    if not fps:
        out.append("rate")
    elif nominal is not None:
        if abs(fps / nominal - 1) > RATE_LIMIT:
            out.append("rate")
    elif core.had(row, "rate"):          # kept from before nominal_fps was always written
        out.append("rate")
    dev = core.num(row, "length_minus_expected_s")
    if dev is not None and abs(dev) > opts.length_tolerance_s:
        out.append("length")
    return out
