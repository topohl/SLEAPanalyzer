"""Video QC: checks of tracked recordings that need the video itself. Report only; nothing here changes a score.

    python -I tools/video_qc <manifest.csv> <out_dir> [--checks frames,objects,sync,clip] [--only REGEX] [--stills]

Checks (each one file next to this one; see their docstrings):
  frames   the video has as many frames as the tracking has rows, at the expected rate and length
  objects  fixed objects (the tracked points of objects, cups) sit on the objects seen in the video, the same way in
           every group
  sync     the tracked animal is where the video shows it, frame by frame (a pose read from another frame shows up)
  clip     the clip covers the whole trial in the raw recording it was cut from (start, end, cut short or not)

The manifest is described in core.py. Writes video_qc_<check>.csv per check (rewritten every 20 recordings, so a run
that stops keeps what it had; with --only, the rows of the other recordings stay) and video_qc.csv: one row per
recording of the manifest with the flags of every check whose table is in out_dir, also those of an earlier run that
this run did not repeat (e.g. the slower sync and clip checks), and the column checks naming them. Delete a check's
table to leave it out. stills/ with --stills. Group comparisons (objects) are made over the recordings of one run.
Tests: tests/test_video_qc_<check>.py.
"""

from __future__ import annotations

import argparse
import importlib
import os
import re
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core  # noqa: E402

ALL = ["frames", "objects", "sync", "clip"]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("out_dir")
    ap.add_argument("--checks", default="frames,objects", help="comma list of " + ",".join(ALL))
    ap.add_argument("--only", help="regular expression on the manifest id")
    ap.add_argument("--stills", action="store_true", help="write review stills for flagged recordings")
    ap.add_argument("--ffmpeg"); ap.add_argument("--ffprobe")
    ap.add_argument("--sync-windows", default="0,0.5,1", help="where the sync check looks: fractions of the video, or all")
    ap.add_argument("--sync-window-s", type=float, default=20, help="length of each sync window, seconds")
    ap.add_argument("--sync-frames", help="sync check: csv (id, frame) of tracking rows to check instead of the windows")
    ap.add_argument("--length-tolerance-s", type=float, default=2, help="frames and clip checks: allowed deviation")
    ap.add_argument("--object-offset-cm", type=float, default=1.0,
                    help="objects check: a group's median point offset may differ this much from the other groups'")
    ap.add_argument("--object-px-per-cm", type=float, default=16.0,
                    help="objects check: converts its cm limits for recordings without a scale")
    opts = ap.parse_args(argv)

    checks = [c.strip() for c in opts.checks.split(",") if c.strip()]
    bad = [c for c in checks if c not in ALL]
    if bad:
        ap.error("unknown check(s): " + ", ".join(bad))
    mods = {c: importlib.import_module(c) for c in checks}
    opts.ffmpeg = core.find_tool("ffmpeg", opts.ffmpeg)
    opts.ffprobe = core.find_tool("ffprobe", opts.ffprobe)
    os.makedirs(opts.out_dir, exist_ok=True)
    opts.stills_dir = os.path.join(opts.out_dir, "stills")
    if opts.stills:
        os.makedirs(opts.stills_dir, exist_ok=True)

    every = core.read_manifest(opts.manifest)
    recs = [r for r in every if re.search(opts.only, r["id"])] if opts.only else every
    print(f"{len(recs)} recordings, checks: {', '.join(checks)}", flush=True)
    if not recs:
        print("nothing to check; nothing written")
        return
    order = {r["id"]: i for i, r in enumerate(every)}

    def table(c):
        return os.path.join(opts.out_dir, f"video_qc_{c}.csv")

    def save(c):
        """This run's rows; with --only, also the rows an earlier run wrote for the recordings this one left out."""
        rows = results[c]
        if opts.only and os.path.exists(table(c)):
            ids = {r["id"] for r in rows}
            rows = rows + [r for r in core.read_table(table(c)) if r["id"] not in ids]
            rows.sort(key=lambda r: order.get(r["id"], len(order)))
        core.write_csv(rows, table(c))

    results = {c: [] for c in checks}
    t0 = time.time()
    if os.path.exists(os.path.join(opts.out_dir, "video_qc.csv")):
        os.remove(os.path.join(opts.out_dir, "video_qc.csv"))   # stale until this run has finished
    for k, rec in enumerate(recs):
        ctx = {"opts": opts}
        try:
            ctx["track"], ctx["rows"] = core.read_tracking(rec["tracking"])
            ctx["video"] = core.Video(rec["video"], opts.ffmpeg, opts.ffprobe)
            ctx["scale"] = core.cm_per_px(rec, ctx["track"])
            ctx["poly"] = core.arena_polygon(ctx["track"], core.names(rec, "outline") or core.names(rec, "corners", "tl tr br bl"))
            setup_error = None
        except Exception as e:  # a missing file stops this recording, not the run
            setup_error = f"{type(e).__name__}: {e}"
        for c in checks:
            row = {"id": rec["id"]}
            if setup_error:
                row.update(flags=["unreadable"], error=setup_error)
            else:
                try:
                    row.update(mods[c].run(rec, ctx))
                except Exception as e:
                    row.update(flags=["error"], error=f"{type(e).__name__}: {e}")
                    if os.environ.get("VIDEO_QC_DEBUG"):
                        traceback.print_exc()
            row.setdefault("flags", [])
            results[c].append(row)
        done = k + 1
        if done % 10 == 0 or done == len(recs):
            print(f"{done} of {len(recs)}  {time.time() - t0:.0f} s", flush=True)
        if done % 20 == 0 and done < len(recs):   # partial tables, without the group comparisons of finish()
            for c in checks:
                save(c)

    for c in checks:
        if hasattr(mods[c], "finish"):
            mods[c].finish(results[c], opts)   # group comparisons: over the recordings of this run
        save(c)

    # flags per recording, from every check's table in out_dir: this run's and those of earlier runs
    flags_of = {}
    for c in ALL:
        if not os.path.exists(table(c)):
            continue
        flags_of[c] = {r["id"]: [f for f in (r.get("flags") or "").split(";") if f] for r in core.read_table(table(c))}
        if c not in checks:
            n = sum(1 for rec in every if rec["id"] in flags_of[c])
            print(f"kept the {c} results of {time.strftime('%Y-%m-%d %H:%M', time.localtime(os.path.getmtime(table(c))))} "
                  f"({n} of {len(every)} recordings)")
    extra = [k for k in every[0].keys() if k not in ("tracking", "video", "raw_video", "raw_dir")]
    combined = []
    for rec in every:
        if not any(rec["id"] in flags_of[c] for c in flags_of):
            continue
        row = {k: rec.get(k, "") for k in extra}
        ran = [c for c in ALL if c in flags_of and rec["id"] in flags_of[c]]
        row["checks"] = ran
        row["flags"] = [f"{c}:{f}" for c in ran for f in flags_of[c][rec["id"]]]
        combined.append(row)
    core.write_csv(combined, os.path.join(opts.out_dir, "video_qc.csv"))
    n_flag = sum(1 for r in combined if r["flags"])
    print(f"{n_flag} of {len(combined)} recordings flagged; results in {opts.out_dir}")


if __name__ == "__main__":
    main()
