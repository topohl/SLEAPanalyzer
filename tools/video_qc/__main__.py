"""Video QC: checks of tracked recordings that need the video itself. Report only; nothing here changes a score.

    python -I tools/video_qc <manifest.csv> <out_dir> [--checks frames,objects,sync,clip] [--only REGEX] [--stills]
                             [--fresh] [--adopt]

Checks (each one file next to this one; see their docstrings):
  frames   the video has as many frames as the tracking has rows, at the expected rate and length
  objects  fixed objects (the tracked points of objects, cups) sit on the objects seen in the video, the same way in
           every group
  sync     the tracked animal is where the video shows it, frame by frame (a pose read from another frame shows up)
  clip     the clip covers the whole trial in the raw recording it was cut from (start, end, cut short or not)

The manifest is described in core.py. Writes video_qc_<check>.csv per check (with --only, the rows of the other
recordings stay) and video_qc.csv: one row per recording of the manifest with the flags of every check whose table is
in out_dir, also those of an earlier run that this run did not repeat (e.g. the slower sync and clip checks), and the
column checks naming them. Delete a check's table to leave it out. stills/ with --stills: one per flagged recording.

Runs are incremental. Every row carries 'inputs', what its measurements depend on, as data:options:code: the content
of the tracking file, size and time of the video (and of raw_video) and the manifest row; the options that change the
check's measurements; the results versions of the check and of core.py (VERSION in each file, raised when a change can
alter measurements; tests/test_video_qc_versions.py asks for that decision whenever the code changes). A run keeps the
rows in out_dir whose inputs are unchanged and checks the other recordings, saying why ("N reused, M to check (data 2,
new 1)"); rows that ended in an error are always checked again; --fresh checks them all. The tables are rewritten every
20 recordings or minute, so a run that stops is continued by the next one, and a new release re-checks only what
changed. Not part of the inputs: the files below raw_dir that the clip check searches (--fresh or --only after adding
raw recordings there) and the ffmpeg version. --adopt takes rows written by an earlier version of this tool (no inputs,
or inputs in one part) as checked when the table is newer than the recording's tracking and video: once, to keep a
full run of that version.
Flags are not measurements: each check makes them from a row's columns (flags()), for kept rows as for new ones in
every run, so a changed limit (--length-tolerance-s) or flag rule needs no new run. A row that gains a flag that way
has no still until it is checked again (review.py draws review sheets from the video for any flag). Group comparisons
(objects) are made over all rows of the table, kept and new.
Tests: tests/test_video_qc_<check>.py, tests/test_video_qc_runs.py, tests/test_video_qc_versions.py.

A check module has run(rec, ctx) -> row and VERSION, and may have: flags(row, rec, opts), the row's flags from its
columns (numbers in rows of this run, text in kept rows: core.num reads both; rec is None for a row whose recording
left the manifest); start(opts) at the start of a run; inputs(rec, opts), more of what a row depends on; reuse(row,
opts) -> False when an unchanged row must be checked again all the same; save(rows, opts), its other files, written
with its table; finish(rows, opts), group comparisons over the whole table, which may return ids to check again (a
still missing). Rows kept from an earlier run carry '_kept' and their values as text.
"""

from __future__ import annotations

import argparse
import collections
import csv
import importlib
import os
import re
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core  # noqa: E402

ALL = ["frames", "objects", "sync", "clip"]
# options that change a check's measurements; an option named neither here nor in NOT_RESULTS counts for every check
OPTIONS = {"objects": ["object_px_per_cm"], "sync": ["sync_windows", "sync_window_s"]}
NOT_RESULTS = {"manifest", "out_dir", "checks", "only", "stills", "ffmpeg", "ffprobe", "fresh", "adopt",
               "sync_frames",         # its rows for the recording count instead (sync.inputs)
               "length_tolerance_s",  # a flag limit: flags are made again from the columns in every run
               "object_offset_cm"}    # only the group comparison, made again over all rows in every run
FAILED = {"error", "unreadable"}      # flags of a row whose check failed: checked again in every run
WHY = {"data": "data", "options": "options", "code": "code", "earlier version": "earlier version",
       "error": "error before", "reuse": "asked by the check", "new": "new"}
SAVE_EVERY, SAVE_S = 20, 60           # partial tables: every this many recordings or seconds


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("out_dir")
    ap.add_argument("--checks", default="frames,objects", help="comma list of " + ",".join(ALL))
    ap.add_argument("--only", help="regular expression on the manifest id")
    ap.add_argument("--stills", action="store_true", help="write review stills for flagged recordings")
    ap.add_argument("--fresh", action="store_true", help="check every recording again, also those with unchanged inputs")
    ap.add_argument("--adopt", action="store_true",
                    help="take rows without inputs (an earlier version's run) as checked when newer than their files")
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
    given = dict(vars(opts))

    checks = [c.strip() for c in opts.checks.split(",") if c.strip()]
    bad = [c for c in checks if c not in ALL]
    if bad:
        ap.error("unknown check(s): " + ", ".join(bad))
    if opts.sync_frames and "sync" in checks and not os.path.isfile(opts.sync_frames):
        ap.error("--sync-frames: no such file: " + opts.sync_frames)
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
    every_id = {r["id"]: r for r in every}
    by_id = {r["id"]: r for r in recs}
    run_ids = set(by_id)

    def table(c):
        return os.path.join(opts.out_dir, f"video_qc_{c}.csv")

    # what each result depends on, and the rows of earlier runs that stay
    t0 = time.time()
    base = {r["id"]: core.recording_inputs(r) for r in recs}
    named = {k for v in OPTIONS.values() for k in v}
    sig, results, others = {}, {}, {}
    for c in checks:
        m = mods[c]
        if hasattr(m, "start"):
            m.start(opts)
        options = {k: v for k, v in sorted(given.items())
                   if k not in NOT_RESULTS and (k in OPTIONS.get(c, ()) or k not in named)}
        sig[c] = {r["id"]: core.row_inputs([c, base[r["id"]], m.inputs(r, opts) if hasattr(m, "inputs") else None],
                                           options, m.VERSION) for r in recs}
        earlier = read_earlier(table(c))
        made = os.path.getmtime(table(c)) if earlier else 0
        # with --only, the rows of the other recordings stay as they are (their flags made again below)
        others[c] = [dict(r, _kept=True) for r in earlier if r["id"] not in run_ids] if opts.only else []
        results[c] = {}
        why = collections.Counter()
        for r in earlier:
            rid = r["id"]
            if opts.fresh or rid not in run_ids:
                continue
            old = r.get("inputs") or ""
            if opts.adopt and old.count(":") != 2 and r.get("_whole") and newer(made, by_id[rid]):
                old = r["inputs"] = sig[c][rid]
            if old != sig[c][rid]:
                why[core.changed_part(old, sig[c][rid])] += 1
                continue
            if FAILED & set(r["flags"]):
                why["error"] += 1
                continue
            r = dict(r, _kept=True)
            if hasattr(m, "reuse") and not m.reuse(r, opts):
                why["reuse"] += 1
                continue
            results[c][rid] = r
        why["new"] = len(recs) - len(results[c]) - sum(why.values())
        for r in list(results[c].values()) + others[c]:
            reflag(m, r, every_id.get(r["id"]), opts)
        told = "" if opts.fresh else ", ".join(f"{WHY[k]} {n}" for k, n in why.items() if n)
        print(f"{c}: {len(results[c])} reused, {len(recs) - len(results[c])} to check" + (f" ({told})" if told else ""),
              flush=True)
    print(f"  (inputs of {len(recs)} recordings and the earlier tables read in {time.time() - t0:.0f} s)", flush=True)

    def rows_of(c):
        """The check's table: this run's rows (kept and checked) and, with --only, those of the other recordings."""
        rows = [results[c][r["id"]] for r in recs if r["id"] in results[c]] + others[c]
        rows.sort(key=lambda r: order.get(r["id"], len(order)))
        return rows

    def save(c):
        rows = rows_of(c)
        if hasattr(mods[c], "save"):
            mods[c].save(rows, opts)       # first, so that the table never lists a row whose other files are missing
        core.write_csv(rows, table(c), last=("inputs",))

    def check(rec, which):
        """Run the checks `which` on one recording; its old stills go first."""
        for c in which:
            old = os.path.join(opts.stills_dir, c, core.still_name(rec["id"]))
            if os.path.exists(old):
                os.remove(old)
        try:
            ctx = core.open_recording(rec, opts)
            setup_error = None
        except Exception as e:  # a missing file stops this recording, not the run
            ctx, setup_error = None, f"{type(e).__name__}: {e}"
        for c in which:
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
            row["inputs"] = sig[c][rec["id"]]
            results[c][rec["id"]] = row

    t0 = last = time.time()
    if os.path.exists(os.path.join(opts.out_dir, "video_qc.csv")):
        os.remove(os.path.join(opts.out_dir, "video_qc.csv"))   # stale until this run has finished
    todo = [r for r in recs if any(r["id"] not in results[c] for c in checks)]
    for k, rec in enumerate(todo):
        check(rec, [c for c in checks if rec["id"] not in results[c]])
        done = k + 1
        if done % 10 == 0 or done == len(todo):
            print(f"{done} of {len(todo)}  {time.time() - t0:.0f} s", flush=True)
        if done < len(todo) and (done % SAVE_EVERY == 0 or time.time() - last >= SAVE_S):
            try:                            # partial tables, without the group comparisons of finish()
                for c in checks:
                    save(c)
            except OSError as e:
                print(f"partial tables not written ({e}); the run goes on", flush=True)
            last = time.time()

    for c in checks:
        if hasattr(mods[c], "finish"):
            again = set(mods[c].finish(rows_of(c), opts) or ())
            redo = [r for r in recs if r["id"] in again]
            if redo:
                print(f"{c}: {len(redo)} flagged recordings without a still checked again", flush=True)
                for rec in redo:
                    check(rec, [c])
                mods[c].finish(rows_of(c), opts)
        save(c)
        # a still only for a flagged row of the table
        folder = os.path.join(opts.stills_dir, c)
        if os.path.isdir(folder):
            keep = {core.still_name(r["id"]) for r in rows_of(c) if r.get("flags")}
            for f in os.listdir(folder):
                if f.lower().endswith(".jpg") and f not in keep:
                    os.remove(os.path.join(folder, f))

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


def reflag(m, row, rec, opts):
    """A kept row's flags, made again from its columns by the check (when it has flags(); not for a failed row)."""
    if hasattr(m, "flags") and not FAILED & set(row.get("flags") or []):
        row["flags"] = m.flags(row, rec, opts)


def read_earlier(path):
    """Rows of a check's table from an earlier run, flags as a list ([] when there is no table or it cannot be read).
    A row with more or fewer fields than the header loses its inputs, so that it is checked again."""
    if not os.path.exists(path):
        return []
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
    except (OSError, UnicodeDecodeError, csv.Error) as e:
        print(f"cannot read {path} ({type(e).__name__}: {e}): its recordings are checked again")
        return []
    out = []
    for r in rows:
        whole = None not in r and None not in r.values()
        r = {k: ("" if v is None else v) for k, v in r.items() if k is not None}
        if not r.get("id"):
            continue
        if not whole:
            r["inputs"] = ""
        r["flags"] = [f for f in (r.get("flags") or "").split(";") if f]
        r["_whole"] = whole
        out.append(r)
    return out


def newer(made, rec):
    """The table was written after the recording's tracking and video last changed (--adopt)."""
    try:
        return made > max(os.path.getmtime(rec["tracking"]), os.path.getmtime(core.field(rec, "video", "")))
    except (OSError, ValueError):
        return False


if __name__ == "__main__":
    main()
