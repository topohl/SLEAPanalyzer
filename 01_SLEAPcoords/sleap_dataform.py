"""Merge SLEAP animal + geom coordinate exports into DLCAnalyzer-readable CSVs.

The logic lives here rather than in DLCA_DataForm.ipynb so it can be run from a
command line, imported and tested -- and so a stale notebook buffer cannot
silently execute an old version of it.

    python sleap_dataform.py --batch B6 --assay SocP --dry-run
    python sleap_dataform.py --batch B6 --assay SocP --overwrite

Output format is unchanged from the original notebook: three header rows, then
one row per frame of ``frame, x, y, likelihood`` per bodypart. SLEAP's
``analysis_locs`` export carries no confidence column, so ``likelihood`` is 1
throughout.

Two design notes, both learned the hard way:

* The merge is done **in memory**. Writing ~650 MB of intermediate CSVs to the
  network share and reading them straight back produced a truncated output
  (one file came out at 1 kB instead of 12 MB) -- a write-then-read race that
  never appears on a local disk. The intermediates have no downstream
  consumer, so they are not written unless ``--keep-merged`` asks for them.

* Output filenames match the existing corpus exactly, including the legacy
  ``.merged_locs`` infix. Renaming them created a *parallel* set of 57 files
  beside the originals, which the overwrite guard could not see (nothing it
  planned existed yet) and which would have left the downstream staging glob
  picking up both.
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #

BEHAVIOR_ROOT = Path(
    r"S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Raw Data\Behavior"
)

#: Phases per assay. An empty tuple means the assay has a single unphased
#: recording per animal, so files are named ``<CODE>`` and land flat in
#: ``formatted/``.
PHASES: dict[str, tuple[str, ...]] = {
    "SocP": ("HAB", "S1", "S2"),
    "NOR": ("HAB", "NOV"),
    "EPM": (),
    "OFT": (),
}

#: Legacy infix in the existing formatted filenames. Kept so re-running
#: overwrites in place rather than creating a second set of files.
NAME_INFIX = ".merged_locs"

#: CRLF matches the existing corpus, so re-running reproduces those files
#: byte-for-byte rather than differing only by line ending. R's read.csv and
#: DLCAnalyzer accept either; the point is consistency within a file.
LINE_TERM = "\r\n"

#: A phase pair whose mean absolute coordinate difference is below this (px)
#: is almost certainly the same footage twice. Genuine sessions in this dataset
#: differ by 18-201 px; the one known duplicate differed by 1.2e-4.
DUPLICATE_PX = 1.0


# --------------------------------------------------------------------------- #
# filename parsing
# --------------------------------------------------------------------------- #

def parse_key(filename: str, phases: tuple[str, ...]) -> tuple[str, str]:
    """Return ``(code, phase)`` parsed from a SLEAP export filename.

    Matches on content, not on underscore position::

        E9_B6_SocP-animal_3.004_E9_SIS_B6_D3U7_S1.analysis_locs.csv
                                          ^^^^ ^^
                                          code phase

    The original notebook used ``filename.split("_", 7)[-1]``, which on
    ``B6_A3R3_S1.merged_locs`` yields ``"locs"`` -- so every output went to the
    same name and overwrote the last. Counting underscores breaks the moment
    the export prefix changes, and it fails silently. This raises instead.
    """
    stem = Path(filename).name
    if phases:
        pattern = r"_([A-Z0-9]{4})_(" + "|".join(phases) + r")\."
        hits = re.findall(pattern, stem)
        if len(hits) != 1:
            raise ValueError(
                f"expected exactly one <CODE>_<PHASE> in {stem!r}, "
                f"found {len(hits)}: {hits}"
            )
        return hits[0]
    hits = re.findall(r"_([A-Z0-9]{4})\.", stem)
    if len(hits) != 1:
        raise ValueError(
            f"expected exactly one <CODE> in {stem!r}, found {len(hits)}: {hits}"
        )
    return (hits[0], "")


def index_folder(folder: Path, phases: tuple[str, ...]) -> dict[tuple[str, str], Path]:
    """Map ``(code, phase)`` to path for every CSV in *folder*, refusing duplicates."""
    out: dict[tuple[str, str], Path] = {}
    for path in sorted(folder.glob("*.csv")):
        key = parse_key(path.name, phases)
        if key in out:
            raise ValueError(
                f"two files map to {key} in {folder.name}:\n"
                f"  {out[key].name}\n  {path.name}"
            )
        out[key] = path
    return out


def stem_for(code: str, phase: str) -> str:
    return f"{code}_{phase}" if phase else code


def output_path(formatted_root: Path, code: str, phase: str) -> Path:
    """Where the formatted CSV for one recording belongs.

    Matches the existing corpus exactly: ``<CODE>_<PHASE>.merged_locs_formatted.csv``
    inside ``formatted/<PHASE>/``, or flat for unphased assays.
    """
    name = f"{stem_for(code, phase)}{NAME_INFIX}_formatted.csv"
    return (formatted_root / phase / name) if phase else (formatted_root / name)


# --------------------------------------------------------------------------- #
# formatting
# --------------------------------------------------------------------------- #

def dlc_header(bodyparts: list[str]) -> list[list[str]]:
    """The three DLCAnalyzer header rows for a set of bodyparts.

    Row 1 reproduces the ``column_N`` / ``column_N_new`` names the original
    script emitted, so output stays byte-identical to the existing corpus.
    DLCAnalyzer reads rows 2 and 3; row 1 is the scorer row and is not used.
    """
    row1 = ["column_1"]
    for i in range(1, len(bodyparts) + 1):
        row1 += [f"column_{2 * i}", f"column_{2 * i + 1}", f"column_{2 * i + 1}_new"]
    row2 = ["bodyparts"] + [bp for bp in bodyparts for _ in range(3)]
    row3 = ["coords"] + ["x", "y", "likelihood"] * len(bodyparts)
    return [row1, row2, row3]


def bodyparts_of(df: pd.DataFrame, label: str) -> list[str]:
    """Bodypart names, in order, from ``<bodypart>_x`` / ``<bodypart>_y`` columns.

    Reading the count from the header is what the original per-assay magic
    numbers (+20 SocP, +18 OFT, +26 EPM) were standing in for.
    """
    if len(df.columns) % 2:
        raise ValueError(
            f"{label}: odd column count {len(df.columns)}; expected x/y pairs"
        )
    names: list[str] = []
    seen: set[str] = set()
    for col in df.columns:
        bp, _, axis = col.rpartition("_")
        if axis not in ("x", "y"):
            raise ValueError(f"{label}: column {col!r} does not end in _x or _y")
        if bp not in seen:
            seen.add(bp)
            names.append(bp)
    for bp in names:
        for axis in ("x", "y"):
            if f"{bp}_{axis}" not in df.columns:
                raise ValueError(f"{label}: {bp} is missing its _{axis} column")
    return names


def process_one(geom_path: Path, animal_path: Path, dst: Path,
                keep_merged: Path | None = None) -> dict:
    """Merge one geom/animal pair in memory and write the formatted CSV.

    Geometry columns come first, then the animal bodyparts -- the order
    DLCAnalyzer's zone definitions expect.
    """
    geom = pd.read_csv(geom_path)
    animal = pd.read_csv(animal_path)
    if len(geom) != len(animal):
        raise ValueError(
            f"{dst.name}: geom has {len(geom)} rows, animal has {len(animal)}"
        )
    if len(geom) == 0:
        raise ValueError(f"{dst.name}: no data rows in {geom_path.name}")

    merged = pd.concat([geom, animal], axis=1)
    if keep_merged is not None:
        keep_merged.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(keep_merged, index=False)

    bodyparts = bodyparts_of(merged, dst.name)
    body = pd.DataFrame({"frame": range(len(merged))})
    for bp in bodyparts:
        body[f"{bp}_x"] = merged[f"{bp}_x"].values
        body[f"{bp}_y"] = merged[f"{bp}_y"].values
        body[f"{bp}_likelihood"] = 1

    expected = 1 + 3 * len(bodyparts)
    if body.shape[1] != expected:
        raise AssertionError(
            f"{dst.name}: built {body.shape[1]} columns, expected {expected}"
        )

    # One explicit terminator for both the hand-written header and the pandas
    # body. Left implicit, the header gets "\n" while pandas writes "\r\n" on
    # Windows, leaving the file internally inconsistent.
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", newline="") as fh:
        for row in dlc_header(bodyparts):
            fh.write(",".join(row) + LINE_TERM)
        # pandas < 1.5 calls this argument ``line_terminator``; newer releases
        # renamed it to ``lineterminator`` and pandas 2 removed the old alias.
        # Select by capability so the SLEAP Python 3.7 environment (pandas
        # 1.3.5) and current environments produce the same bytes.
        terminator_arg = (
            "lineterminator"
            if "lineterminator" in inspect.signature(body.to_csv).parameters
            else "line_terminator"
        )
        body.to_csv(
            fh, index=False, header=False, **{terminator_arg: LINE_TERM}
        )

    # Read the row count back. On a network share a write can be visible
    # before it is complete, and a short read is how a 12 MB file silently
    # became 1 kB.
    with open(dst, "rb") as fh:
        written = sum(1 for _ in fh) - 3
    if written != len(merged):
        raise IOError(
            f"{dst.name}: wrote {len(merged)} rows but {written} readable back; "
            "the write did not complete"
        )

    return {"bodyparts": len(bodyparts), "cols": expected, "frames": len(merged)}


# --------------------------------------------------------------------------- #
# pipeline
# --------------------------------------------------------------------------- #

def run(batch: str, assay: str, *, root: Path = BEHAVIOR_ROOT,
        overwrite: bool = False, dry_run: bool = False,
        formatted_dir: Path | None = None, keep_merged: bool = False,
        verbose: bool = True) -> pd.DataFrame:
    """Merge, format and verify one batch/assay. Returns the per-file report."""
    if assay not in PHASES:
        raise ValueError(f"unknown assay {assay!r}; add it to PHASES")
    phases = PHASES[assay]

    sleap = root / batch / assay / "SLEAP"
    animal, geom = sleap / "animal", sleap / "geom"
    for d in (animal, geom):
        if not d.is_dir():
            raise FileNotFoundError(f"missing input folder: {d}")
    formatted_root = formatted_dir or (sleap / "formatted")
    merged_root = (sleap / "merged") if keep_merged else None

    say = print if verbose else (lambda *a, **k: None)
    say(f"{batch} / {assay}")
    say(f"  animal : {len(list(animal.glob('*.csv'))):3d} csv")
    say(f"  geom   : {len(list(geom.glob('*.csv'))):3d} csv")
    say(f"  phases : {phases or '(none)'}")
    say(f"  output : {formatted_root}")

    animal_idx = index_folder(animal, phases)
    geom_idx = index_folder(geom, phases)
    only_animal = sorted(set(animal_idx) - set(geom_idx))
    only_geom = sorted(set(geom_idx) - set(animal_idx))
    paired = sorted(set(animal_idx) & set(geom_idx))

    say(f"\npaired: {len(paired)}   animal-only: {len(only_animal)}   "
        f"geom-only: {len(only_geom)}")
    for k in only_animal:
        say(f"  !! no geom for   {k}  ({animal_idx[k].name})")
    for k in only_geom:
        say(f"  !! no animal for {k}  ({geom_idx[k].name})")

    # ---- plan every path before writing anything --------------------------
    plan: dict[tuple[str, str], Path] = {}
    for code, phase in paired:
        dst = output_path(formatted_root, code, phase)
        if dst in plan.values():
            raise RuntimeError(f"name collision: ({code},{phase}) maps to {dst}")
        plan[(code, phase)] = dst

    existing = sorted(d for d in plan.values() if d.exists())

    # Anything already in formatted/ that this run will NOT replace. A leftover
    # from an older naming scheme would be picked up twice by the downstream
    # staging glob, so say so rather than leave it to be discovered later.
    planned = set(plan.values())
    orphans = [p for p in formatted_root.rglob("*.csv") if p not in planned] \
        if formatted_root.exists() else []
    if orphans:
        say(f"\n  !! {len(orphans)} file(s) in formatted/ will NOT be replaced by "
            "this run and may be stale:")
        for p in orphans[:10]:
            say(f"       {p.relative_to(formatted_root)}")
        if len(orphans) > 10:
            say(f"       ... and {len(orphans) - 10} more")

    # A dry run reports and exits cleanly; it must not trip the overwrite guard,
    # since its whole purpose is to tell you what a real run would do.
    if dry_run:
        say(f"\n[dry run] would write {len(plan)} file(s); "
            f"{len(existing)} already exist and would be replaced.")
        if existing and not overwrite:
            say("           a real run needs --overwrite. Nothing changed.")
        else:
            say("           nothing changed.")
        return pd.DataFrame(
            [{"code": c, "phase": p, "out": str(d)} for (c, p), d in sorted(plan.items())]
        )

    if existing and not overwrite:
        raise RuntimeError(
            f"{len(existing)} output file(s) already exist, e.g. {existing[0].name}.\n"
            "Pass --overwrite if replacing them is intended."
        )

    # ---- merge in memory and format --------------------------------------
    rows, problems = [], []
    for key in paired:
        code, phase = key
        km = (merged_root / f"{stem_for(code, phase)}.merged.csv") if merged_root else None
        try:
            info = process_one(geom_idx[key], animal_idx[key], plan[key], km)
        except (ValueError, IOError, AssertionError) as exc:
            problems.append(f"{key}: {exc}")
            continue
        rows.append({"code": code, "phase": phase, **info, "out": plan[key].name})

    for p in problems:
        say(f"  !! {p}")
    report = pd.DataFrame(rows).sort_values(["phase", "code"]).reset_index(drop=True)
    say(f"\nformatted {len(report)} of {len(paired)} file(s)")
    if len(report):
        n_bp = sorted(report["bodyparts"].unique())
        n_fr = sorted(report["frames"].unique())
        say(f"  bodyparts per file: {n_bp}"
            + ("  <-- INCONSISTENT" if len(n_bp) > 1 else ""))
        say(f"  frames per file:    {n_fr}"
            + ("  <-- differing lengths" if len(n_fr) > 1 else ""))

    ok = {k: plan[k] for k in plan if plan[k].exists()}
    verify(ok, phases, say)
    return report


def verify(plan: dict[tuple[str, str], Path], phases: tuple[str, ...], say=print) -> None:
    """Confirm outputs are distinct, and flag duplicated recordings."""
    written = sorted(plan.values())
    if not written:
        say("\nnothing written; skipping verification")
        return

    by_hash: dict[str, list[str]] = {}
    for p in written:
        by_hash.setdefault(hashlib.sha256(p.read_bytes()).hexdigest(), []).append(p.name)
    collisions = {h: v for h, v in by_hash.items() if len(v) > 1}

    say(f"\n{len(written)} file(s) written, {len(by_hash)} distinct by content")
    if collisions:
        say("  !! byte-identical outputs:")
        for names in collisions.values():
            say(f"     {names}")
    else:
        say("  all outputs distinct")

    if not phases:
        return

    # Two exports of one recording agree to ~1e-4 px; two genuine sessions
    # differ by tens of pixels. This is how the OR643 duplication was caught --
    # see SISanalyzer/exp9_publication/data/DATA_ISSUES.md.
    by_code: dict[str, dict[str, Path]] = {}
    for (code, phase), dst in plan.items():
        by_code.setdefault(code, {})[phase] = dst

    flagged, unknown = [], []
    say("\nper-animal agreement between phases (mean |diff|, px):")
    for code, d in sorted(by_code.items()):
        for p1, p2 in combinations(sorted(d), 2):
            a = pd.read_csv(d[p1], skiprows=3, header=None).iloc[:, 1:]
            b = pd.read_csv(d[p2], skiprows=3, header=None).iloc[:, 1:]
            n = min(len(a), len(b))
            diff = np.abs(a.iloc[:n].values - b.iloc[:n].values)
            # nanmean, not mean: the geom columns carry NaN where a landmark
            # was not detected, and a plain mean returns NaN. Worse, NaN <
            # DUPLICATE_PX is False, so a NaN silently made this check pass and
            # report "no duplicates" for every pair.
            m = float(np.nanmean(diff)) if np.any(~np.isnan(diff)) else float("nan")
            dup = (m < DUPLICATE_PX) if not np.isnan(m) else False
            if np.isnan(m):
                unknown.append((code, p1, p2))
            if dup:
                flagged.append((code, p1, p2, m))
            say(f"  {code}  {p1} vs {p2}: {m:12.6f}"
                + ("  <-- SAME RECORDING?" if dup else ""))

    if flagged:
        say(f"\n!! {len(flagged)} phase pair(s) look like the same recording:")
        for code, p1, p2, m in flagged:
            say(f"   {code} {p1}/{p2}  mean |diff| = {m:.2e} px")
        say("   Check the source videos before analysing these.")
    else:
        say("\nno duplicated recordings detected")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Merge SLEAP animal + geom exports into DLCAnalyzer CSVs.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--batch", required=True, help="cohort, e.g. B6")
    ap.add_argument("--assay", required=True, choices=sorted(PHASES))
    ap.add_argument("--root", type=Path, default=BEHAVIOR_ROOT,
                    help="behaviour data root")
    ap.add_argument("--formatted-dir", type=Path, default=None,
                    help="write elsewhere instead of <SLEAP>/formatted (for testing)")
    ap.add_argument("--overwrite", action="store_true",
                    help="replace existing formatted files")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be written, change nothing")
    ap.add_argument("--keep-merged", action="store_true",
                    help="also write the intermediate merged CSVs (debugging only)")
    args = ap.parse_args(argv)
    try:
        run(args.batch, args.assay, root=args.root,
            formatted_dir=args.formatted_dir, overwrite=args.overwrite,
            dry_run=args.dry_run, keep_merged=args.keep_merged)
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
