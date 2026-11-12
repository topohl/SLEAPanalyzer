"""Re-run ONE re-tracked recording through the whole SLEAPanalyzer chain, in one command.

    h5 export -> locs.csv (geom + animal) -> merged + formatted CSV -> staged copy -> scoring of that batch -> before/after

Plan only by default: nothing is written. --apply does it. Nothing is deleted: every file that is replaced moves to
``_superseded_<date>/`` (h5 and locs files, the merged csv, staged copies) or ``formatted/old/`` (formatted files; the
audit script 07 skips only that folder name); for EPM, OFT and SocP a second file of the same name archived on the same
day is kept as ``<name>__2.csv``. Anything already up to date is left alone, so a second run changes nothing (EPM,
OFT and SocP: unless the scoring of the installed file is still due; it then scores).

NOR (novel-object trial)
    python tools/retrack_recording.py --batch B5 --assay NOR --code H1E2 --phase NOV --static-landmark objR            # plan
    python tools/retrack_recording.py --batch B5 --assay NOR --code H1E2 --phase NOV --static-landmark objR --apply    # do it

EPM (one trial per animal, no --phase; staged with the neck alias)
    python tools/retrack_recording.py --batch B6 --assay EPM --code R8I2
    python tools/retrack_recording.py --batch B6 --assay EPM --code R8I2 --apply --baseline

OFT (batch 1 names its tracking files by GoPro clip, so B1 needs --clip; --reuse-geom keeps the corners already tracked)
    python tools/retrack_recording.py --batch B1 --assay OFT --code X1H0 --clip GX010804 --reuse-geom
    python tools/retrack_recording.py --batch B1 --assay OFT --code X1H0 --clip GX010804 --reuse-geom --apply --baseline

SocP (--phase HAB, S1 or S2; HAB is scored with socp_hab_all.yaml, S1 and S2 with socp_all.yaml)
    python tools/retrack_recording.py --batch B2 --assay SocP --code G6K2 --phase S2
    python tools/retrack_recording.py --batch B2 --assay SocP --code G6K2 --phase S2 --apply --baseline

NOR, as before: reads the ONE ``*_<CODE>_<PHASE>.analysis.h5`` in geom/ and animal/ (names containing ``_old`` are
ignored; two matches stop it), writes ``..._<CODE>_<PHASE>.analysis_locs.csv`` next to them (``--static-landmark`` holds
a hand-labelled static object at its first labelled position), installs ``formatted/<CODE>_<PHASE>.merged_locs_formatted.csv``
(``sleap_dataform.process_one``; a legacy-labelled file of the same animal such as ``H1E2_NOR...`` moves to formatted/old/),
copies it to ``<config input_dir>/<batch>/<CODE>_<PHASE>.csv``, backs up the batch's output, runs
``DLCA_NOR v1.2.1.R`` for that batch only and reports which animals' rows changed (all others must be identical).

EPM, OFT and SocP:
  * Where the new export goes. Put the arena (geom) and animal exports in
    ``Raw Data/Behavior/<batch>/<assay>/SLEAP/geom/`` and ``.../animal/`` under SLEAP's new default name, which carries
    the video name (never the name of the export in place: SLEAP would overwrite it before it is archived). It must
    contain the animal code (SocP: ``_<CODE>_<PHASE>``), end in ``.analysis.h5`` and not contain ``_old``. Or leave
    the exports anywhere and pass ``--geom-h5``/``--animal-h5``; they are then copied into those folders (under a name
    with ``_<CODE>[_<PHASE>]`` added if theirs would not be found again). The plan prints the folders, an example name,
    and exports of the recording lying in SLEAP/ or SLEAP/predictions, where the tool does not read them.
  * Which export is which. In place: it rebuilds its locs file, or it has a locs file of its own that it does not
    rebuild, or it is older than the chain's locs files ("in place, does not rebuild": used only with ``--rebuild`` or
    when named with ``--geom-h5``/``--animal-h5``). New: no locs file of its own and newer than the chain's. Replaced
    ``.analysis.h5`` and ``_locs.csv`` files move to ``_superseded_<date>/``; other files of the recording stay.
  * Checks, whenever the staged file would change: the stored node names and a single track
    (``sleap_h5.check_export``); a new export's file name and stored video_path name the animal (SocP with the phase;
    OFT B1 the clip) and no other animal (``--accept-name`` overrides); OFT B1: the merged file of ``--clip`` formats
    to the animal's formatted file, so the clip is its own (``--accept-clip``); the length against the protocol
    (light-check rule: 5 min EPM, 10 min OFT and SocP, within 2 s at 29.97 frames/s; ``--accept-length``); geom and
    animal of equal length; the same trial twice as 31_light_checks.R finds it: a byte copy of any staged file of the
    assay, or the body centre within 0.15 cm of any staged file of the assay and batch (every animal and phase) at an
    offset of up to 60 s (``--accept-duplicate``). Reported without stopping: the arena's side ratio (SocP, OFT), the
    arm width through the 60 cm calibration (EPM), how far the arena points moved against the animal's other phases.
  * Files written: the locs csv next to each h5 (``sleap_h5.extract_locs`` with the batch's labels), the merged csv in
    SLEAP/ and the formatted csv (``sleap_dataform.format_merged``: the route that reproduces the existing files byte
    for byte; OFT B1/B2 at 10 significant digits, as written in 2023), both under the names this animal already has;
    every other formatted file that staging would take for the animal moves to formatted/old/ (SocP:
    formatted/old/<PH>_old/ where it exists, and the flat copy in formatted/ of B3/B5 to formatted/old/). The staged
    copy is written as ``08_stage_all_batches.R`` writes it (EPM: with the neck alias, by 08's own ``copy_with_neck``
    through Rscript; ``--apply`` without Rscript needs ``--python-neck-copy``). SocP: the side file
    ``novelLoc<phase>.txt`` must already list the animal; it is checked, not changed.
  * Scoring: a temporary copy of the config with ``batches: [<batch>]`` (SocP also ``phases: [<phase>]``, and
    ``qc_decisions_file`` made absolute). The batch's output is backed up outside the output tree first; OFT's
    all-batches folder is put back after each one-batch run, and only a failures.csv written by the run is read.
    Landmark repair, head restore and the on-cup decisions (SocP) run inside the scoring, so re-scoring redoes them.
    Scoring is due when the staged file in place is not the one the output was scored from: the md5 its
    ``run_manifest.yaml`` records, or (OFT, no manifest) the md5 the tool keeps in ``retrack_scored.tsv`` in the
    output folder; a re-run then scores without a new export. If the scoring code changed since the current output was
    made, or this cannot be told (OFT), every animal of the batch would change: ``--baseline`` scores the batch once on
    the files the output was scored from (if the new file is already staged, with its archived predecessor put back
    for that run) and compares against that, or ``--scoring-repo`` scores from a checkout at the old commit. Without
    either, ``--apply`` stops before changing anything (``--skip-analysis`` installs the files and does not score).
  * Steps the tool does not do are printed at the end with the exact command: review of the light-check decision
    rows of this recording (keyed to the staged file's md5), the on-cup decision rows (SocP; their times refer to the
    old clip), the retrack and leave-out lists, the video clips, a light check of the whole assay and batch (31 looks
    for the same trial twice only among the recordings it selects), the release and the SISanalyzer ingest.
  * If the exports in place do not rebuild the files in the chain (an interrupted run), the plan lists what differs;
    ``--apply --rebuild`` carries it out without a new export.

Testing on copies: ``--raw-root`` (or ``--root``), ``--stage-root``, ``--output-root`` and ``--qc-root`` point the tool,
and the scoring through the temporary config, at other folders.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "01_SLEAPcoords"))
import sleap_dataform as sd  # noqa: E402
import sleap_h5 as sh  # noqa: E402

CONFIG_DIR = REPO / "validation" / "exp9_boris" / "config"
DEFAULT_CONFIG = CONFIG_DIR / "nor_all.yaml"
NOR_SCRIPT = REPO / "02_SLEAPanalzyer" / "DLCA_NOR v1.2.1.R"
STAGE_SCRIPT = Path("validation") / "exp9_boris" / "scripts" / "08_stage_all_batches.R"
LIGHT_CHECKS = REPO / "validation" / "exp9_boris" / "scripts" / "31_light_checks.R"
RELEASE_SCRIPT = REPO / "tools" / "build_exp9_behavior_bundle.ps1"
DECISIONS = CONFIG_DIR / "light_check_decisions.csv"
QC_ROOT = Path(r"S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\Behavior\SLEAP\qc\light_checks")
DEFAULT_RSCRIPT = Path(r"C:\Users\topohl\AppData\Local\Programs\R\R-4.5.1\bin\x64\Rscript.exe")
TODAY = dt.date.today().isoformat()

VIDEO_FPS = 30000 / 1001          # the cameras' rate, as in 31_light_checks.R
LENGTH_LIMIT_S = 2.0              # 31_light_checks.R LIMIT$length_s
DUPLICATE_CM = 0.15               # 31_light_checks.R LIMIT$duplicate_cm
DUPLICATE_LAG = 1800              # frames (60 s), 31_light_checks.R LIMIT$duplicate_lag
DUPLICATE_MIN_OVERLAP = 300       # frames
ASPECT_LIMIT = 0.03               # 31_light_checks.R LIMIT$aspect
ARM_CM = (4.5, 5.5)               # 31_light_checks.R LIMIT$arm_cm
SCORING_CODE = ["02_SLEAPanalzyer/DLCAnalyzer_Functions_final.R", "02_SLEAPanalzyer/Behavioral_Metrics_Phase1.R",
                "02_SLEAPanalzyer/core"]
ALIASES = {"M8P1": "M8P3"}        # 08_stage_all_batches.R


# --------------------------------------------------------------------------- #
# shared helpers
# --------------------------------------------------------------------------- #

def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


def config_value(config: Path, key: str) -> str:
    return text_value(config.read_text(encoding="utf-8"), key, config.name)


def text_value(text: str, key: str, label: str = "config", required: bool = True) -> str | None:
    m = re.search(rf"^{key}:\s*(.+?)\s*$", text, re.M)
    if not m:
        if required:
            raise SystemExit(f"{label} has no '{key}:' line")
        return None
    return m.group(1).strip("'\"")


def one_h5(folder: Path, code: str, phase: str, explicit: Path | None) -> Path:
    if explicit:
        return explicit
    hits = sorted(p for p in folder.glob(f"*_{code}_{phase}.analysis.h5") if "_old" not in p.name)
    if len(hits) != 1:
        raise SystemExit(f"expected exactly one *_{code}_{phase}.analysis.h5 in {folder}, found {[p.name for p in hits]}; "
                         "name it with --geom-h5 / --animal-h5")
    return hits[0]


def install(new: Path, dst: Path, archive_dir: Path, apply: bool, log: list[str]) -> bool:
    """Put *new* at *dst*; a different existing file moves to *archive_dir* first. True if anything changed."""
    if dst.exists() and sha(dst) == sha(new):
        log.append(f"  unchanged   {dst}")
        return False
    log.append(f"  {'install' if not dst.exists() else 'replace '}   {dst}" + (f"   (old -> {archive_dir})" if dst.exists() else ""))
    if apply:
        if dst.exists():
            archive_dir.mkdir(parents=True, exist_ok=True)
            target = archive_dir / dst.name
            if target.exists():
                raise SystemExit(f"{target} already exists; refusing to overwrite an archived file")
            shutil.move(str(dst), str(target))
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = dst.with_name(dst.name + ".installing")
        shutil.copyfile(new, tmp)
        if sha(tmp) != sha(new):
            tmp.unlink()
            raise SystemExit(f"hash mismatch writing {dst}")
        os.replace(tmp, dst)
    return True


def read_rows(csv_path: Path) -> dict[str, dict[str, str]]:
    with open(csv_path, newline="", encoding="utf-8") as fh:
        return {r["Code"]: r for r in csv.DictReader(fh)}


def changed_cells(before: dict, after: dict, ignore: tuple[str, ...] = ()) -> dict:
    out: dict = {}
    for code in sorted(set(before) | set(after)):
        b, a = before.get(code), after.get(code)
        if b is None or a is None:
            out[code] = ["<row missing>"]
            continue
        cols = []
        for c in a:
            if c in ignore:
                continue
            x, y = b.get(c, ""), a[c]
            try:
                same = abs(float(x) - float(y)) <= 1e-9
            except ValueError:
                same = x == y
            if not same:
                cols.append(c)
        if cols:
            out[code] = cols
    return out


def run_scoring(rscript: Path, script: Path, config: Path, cwd: Path) -> subprocess.CompletedProcess:
    """Run one assay script with SLEAP_ANALYZER_CONFIG set to *config*."""
    env = dict(os.environ, SLEAP_ANALYZER_CONFIG=str(config))
    return subprocess.run([str(rscript), str(script)], cwd=cwd, env=env, capture_output=True, text=True)


#: The scoring call; the tests replace it with a stub that never runs R.
SCORER = run_scoring


def retarget(text: str, moves: list[tuple[str, str]]) -> str:
    """Point every path of a config text that starts with one root at another (longest root first)."""
    for old, new in sorted(moves, key=lambda m: -len(m[0])):
        if old != new:
            text = text.replace(old, new)
    return text


# --------------------------------------------------------------------------- #
# NOR (unchanged chain)
# --------------------------------------------------------------------------- #

def run_nor(args: argparse.Namespace) -> int:
    config = args.config or DEFAULT_CONFIG
    sleap = args.raw_root / args.batch / args.assay / "SLEAP"
    geom_dir, animal_dir, fmt_dir = sleap / "geom", sleap / "animal", sleap / "formatted"
    input_dir, output_dir = Path(config_value(config, "input_dir")), Path(config_value(config, "output_dir"))
    moves = []
    if args.stage_root:
        moves.append((input_dir.parent.as_posix(), args.stage_root.as_posix()))
        input_dir = args.stage_root / input_dir.name
    if args.output_root:
        moves.append((output_dir.parent.as_posix(), args.output_root.as_posix()))
        output_dir = args.output_root / output_dir.name
    stage_dir, out_dir = input_dir / args.batch, output_dir / args.batch
    arch = f"_superseded_{TODAY}"
    log: list[str] = []
    mode = "APPLY" if args.apply else "PLAN (nothing is written; add --apply)"
    print(f"{args.batch} {args.assay} {args.code} {args.phase}: {mode}")

    geom_h5 = one_h5(geom_dir, args.code, args.phase, args.geom_h5)
    animal_h5 = one_h5(animal_dir, args.code, args.phase, args.animal_h5)
    print(f"  geom h5   {geom_h5.name}\n  animal h5 {animal_h5.name}")
    if args.static_landmark:
        print(f"  static    {', '.join(args.static_landmark)}")

    work = Path(tempfile.mkdtemp(prefix="retrack_"))
    try:
        # ---- extract and format into a scratch folder first: nothing live is touched until all of it worked
        geom_stem = geom_h5.name[: -len(".analysis.h5")]
        geom_csv = work / f"{geom_stem}.analysis_locs.csv"
        animal_name = (geom_stem.replace("-geom.", "-animal.") if "-geom." in geom_stem else animal_h5.name[: -len(".analysis.h5")]) + ".analysis_locs.csv"
        animal_csv = work / animal_name
        g, notes = sh.extract_locs(geom_h5, "geom", args.assay, tuple(args.static_landmark))
        a, _ = sh.extract_locs(animal_h5, "animal", args.assay)
        for n in notes:
            print("  note     ", n)
        g.to_csv(geom_csv, index=False)
        a.to_csv(animal_csv, index=False)
        fmt_name = f"{args.code}_{args.phase}{sd.NAME_INFIX}_formatted.csv"
        fmt_new = work / fmt_name
        info = sd.process_one(geom_csv, animal_csv, fmt_new)
        print(f"  formatted {info['frames']} frames, {info['bodyparts']} landmarks, {int(g.isna().sum().sum() + a.isna().sum().sum())} NaN after filling")

        # ---- install, archiving what is replaced
        changed = False
        changed |= install(geom_csv, geom_dir / geom_csv.name, geom_dir / arch, args.apply, log)
        changed |= install(animal_csv, animal_dir / animal_csv.name, animal_dir / arch, args.apply, log)
        legacy = [p for p in sorted(fmt_dir.glob(f"{args.code}_*{sd.NAME_INFIX}_formatted.csv"))
                  if p.name != fmt_name and not p.name.startswith(f"{args.code}_HAB.")]
        for p in legacy:                                    # would be staged ahead of the new file (alphabetical, first wins)
            log.append(f"  supersede   {p}   (-> formatted/old/)")
            if args.apply:
                (fmt_dir / "old").mkdir(exist_ok=True)
                if (fmt_dir / "old" / p.name).exists():
                    raise SystemExit(f"{fmt_dir / 'old' / p.name} already exists")
                shutil.move(str(p), str(fmt_dir / "old" / p.name))
            changed = True
        changed |= install(fmt_new, fmt_dir / fmt_name, fmt_dir / "old", args.apply, log)
        staged = stage_dir / f"{args.code}_{args.phase}.csv"
        changed |= install(fmt_new, staged, stage_dir.parent / arch, args.apply, log)
        print("\n".join(log))
        if not changed:
            print("everything is already up to date")
            if args.apply and not args.rerun_analysis:
                return 0                                    # nothing new to analyse; --rerun-analysis forces it
        if args.skip_analysis or not args.apply:
            if not args.apply:
                print(f"  then      back up {out_dir}, run {NOR_SCRIPT.name} for {args.batch}, report changed animals")
            return 0

        # ---- analysis for this batch only, with a backup of the previous output
        # outside the assay output root: the release scripts read <root>/NOR/<batch>/..., and nothing may glob into a backup
        before_dir = out_dir.parent.parent.parent / f"_superseded_{TODAY}_{args.code}" / f"{args.assay}_{args.batch}_before"
        if out_dir.exists() and not before_dir.exists():
            before_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(out_dir, before_dir)
            print(f"  backed up {out_dir} -> {before_dir}")
        before = read_rows(before_dir / "combined_output.csv") if (before_dir / "combined_output.csv").exists() else {}
        cfg = work / "nor_batch.yaml"
        cfg.write_text(retarget(re.sub(r"^batches:.*$", f"batches: [{args.batch}]", config.read_text(encoding="utf-8"), flags=re.M), moves),
                       encoding="utf-8")
        print(f"  running   {NOR_SCRIPT.name} for {args.batch} ...")
        res = SCORER(args.rscript, NOR_SCRIPT, cfg, REPO)
        if res.returncode != 0:
            print(res.stdout[-1500:], res.stderr[-1500:], file=sys.stderr)
            raise SystemExit(f"{NOR_SCRIPT.name} failed (exit {res.returncode}); the previous output is in {before_dir}")
        after = read_rows(out_dir / "combined_output.csv")
        diff = changed_cells(before, after)
        print(f"  analysis  {len(after)} animals; rows that changed against the backup: {dict((k, len(v)) for k, v in diff.items()) or 'none'}")
        others = [k for k in diff if k != args.code]
        if others:
            print(f"  !! animals other than {args.code} changed: {others} -- check before using this output", file=sys.stderr)
            return 2
        row = after.get(args.code, {})
        print(f"  {args.code}: validFraction {row.get('validFraction')}, qcPass {row.get('qcPass')}, contactNov {row.get('contactNov')}, "
              f"contactFam {row.get('contactFam')}, latency {row.get('latency')}")
        print("next: build a release (tools/build_exp9_behavior_bundle.ps1), then the SISanalyzer ingest")
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


# --------------------------------------------------------------------------- #
# EPM, OFT, SocP: what differs between the assays and batches
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Table:
    """A scored-output table compared before/after, keyed by *key* (the first key column starts with the code)."""
    name: str
    key: tuple[str, ...]
    ignore: tuple[str, ...] = ()


#: Columns of the open-field tables that are z-scores within the batch: they move for every animal when one changes.
OFT_COHORT = ("z_center_time_percent", "z_center_entries", "z_center_latency_inverse", "z_mean_wall_distance",
              "oft_center_exploration_score_cohort_z_experimental", "oft_center_exploration_cohort_n")


@dataclass(frozen=True)
class Spec:
    assay: str
    phases: tuple[str, ...]           # () = one recording per animal
    script: str                       # in 02_SLEAPanalzyer/
    configs: dict                     # phase ("" if none) -> config file in validation/exp9_boris/config/
    expected_s: float                 # protocol length (31_light_checks.R EXPECT)
    arena_cm: tuple | None            # (width, height) of the arena polygon tl-tr-br-bl; None = EPM (tl-bl = 60 cm)
    tables: tuple
    neck: bool = False                # staged with the neck alias of spine1 (08_stage_all_batches.R)
    shows: tuple = ()                 # columns printed for the re-tracked animal


SPECS: dict[str, Spec] = {
    "EPM": Spec("EPM", (), "DLCA_EPM v1.0.0.R", {"": "epm_all.yaml"}, 300, None,
                (Table("Report.csv", ("file",)), Table("tracking_qc.csv", ("file",))), neck=True,
                shows=("frames", "duration_s", "bodycentre.total.time", "bodycentre.raw.distance", "nose.dip",
                       "valid_fraction", "qcPass")),
    "OFT": Spec("OFT", (), "DLCA_OFT v1.2.0.R", {"": "oft_all.yaml"}, 600, (49, 49),
                (Table("combined_summary.csv", ("Code",), OFT_COHORT),
                 Table("combined_enhanced_metrics.csv", ("Code",), OFT_COHORT + ("enhanced_qc_flag",)),
                 Table("combined_first_last_epoch_metrics.csv", ("Code", "epoch")),
                 Table("combined_bins.csv", ("Code", "bin_s")),
                 Table("combined_qc.csv", ("file", "stage", "point"))),
                shows=("total_time_s", "distance_cm", "center_time_percent", "center_entries", "corner_time_percent",
                       "immobility_time_s", "mean_wall_distance_cm", "qc_flag")),
    "SocP": Spec("SocP", ("HAB", "S1", "S2"), "DLCA_SocP v.1.1.0.R",
                 {"HAB": "socp_hab_all.yaml", "S1": "socp_all.yaml", "S2": "socp_all.yaml"}, 600, (44, 24),
                 (Table("combined_output.csv", ("Code",)), Table("tracking_qc.csv", ("file",))),
                 shows=("frames", "totalTime", "validFraction", "contactNovel", "contactFamiliar", "contactLeft",
                        "contactRight", "novelLoc")),
}

#: Bodypart names the formatted files use where the batch's locs files say otherwise (open field B1/B2).
FORMAT_RENAME = {("OFT", "B1"): {"center": "bodycentre"}, ("OFT", "B2"): {"center": "bodycentre"}}

#: Significant digits of the formatted files where the batch was written rounded (open field B1/B2, 2023 writer).
FORMAT_DIGITS = {("OFT", "B1"): 10, ("OFT", "B2"): 10}

#: Where the scoring writes no run_manifest.yaml (open field), the tool keeps the md5 of the staged file it scored here,
#: in the batch's output folder.
MARKER = "retrack_scored.tsv"

#: The trial names of the retrack and leave-out lists, as phases.
TRIAL_PHASE = {"habituation": "HAB", "hab": "HAB", "phase 1": "S1", "s1": "S1", "phase 2": "S2", "s2": "S2"}


def default_formatted(assay: str, batch: str, code: str, phase: str, g: str) -> str:
    """The formatted name of the batch's convention, used only when the animal has no formatted file yet."""
    if assay == "SocP":
        return f"{code}_{phase}{sd.NAME_INFIX}_formatted.csv"
    if assay == "EPM":
        return {"B1": f"{code}.merged_locs_formatted.csv", "B2": f"{g}.merged_locs_formatted.csv",
                "B3": f"E9_SIS_B3_EPM_B3_{code}_formatted.csv", "B4": f"{g}.merged_locs_formatted.csv"
                }.get(batch, f"{code}_EPM.merged_locs_formatted.csv")
    return {"B1": f"E9_B1_OFT_{code}.merged_locs_formatted.csv", "B2": f"E9_SIS_B2_{code}.merged_locs_formatted.csv",
            "B3": f"{code}.merged_locs_formatted.csv"}.get(batch, f"{code}_OFT.merged_locs_formatted.csv")


def default_merged(assay: str, batch: str, code: str, phase: str, g: str) -> tuple[str, bool]:
    """The merged-csv name of the batch's convention, and whether it follows the geom locs stem ({g})."""
    if assay == "SocP":
        return (f"{g}.merged_locs.csv", True) if batch in ("B3", "B5") else (f"{batch}_{code}_{phase}.merged_locs.csv", False)
    if assay == "EPM":
        return (f"{code}.merged_locs.csv", False) if batch == "B1" else (f"{g}.merged_locs.csv", True)
    if batch in ("B1", "B2"):                                # E9_B1_OFT_<clip> / E9_SIS_B2_<CODE>: the geom locs stem
        return f"{g}.merged_locs.csv", True
    return (f"{code}.merged_locs.csv", False) if batch == "B3" else (f"{code}_OFT.merged_locs.csv", False)


def load_codes(stage_root: Path) -> list[str]:
    p = stage_root / "animalIDCode.txt"
    if not p.exists():
        return []
    with open(p, newline="", encoding="utf-8") as fh:
        return [r["Code"].strip() for r in csv.DictReader(fh, delimiter="\t") if r.get("Code")]


def resolve_code(stem: str, codes: list[str], target: str) -> str | None:
    """The animal 08_stage_all_batches.R stages a file as (its resolve_code: substring, lookalikes, aliases)."""
    pool = codes or [target]
    hit = [c for c in pool if c in stem]
    if len(hit) == 1:
        return hit[0]
    if len(hit) > 1:
        return None
    alt = stem.replace("0", "O").replace("l", "I")
    hit = [c for c in pool if c in alt]
    if len(hit) == 1:
        return hit[0]
    alias = [ALIASES[t] for t in re.findall(r"[A-Za-z0-9]{4}", stem) if t in ALIASES]
    return alias[0] if len(alias) == 1 else None


def token(word: str) -> re.Pattern:
    """*word* as a whole name part (not inside a longer run of letters and digits)."""
    return re.compile(rf"(?<![A-Za-z0-9]){re.escape(word)}(?![A-Za-z0-9])")


def when(p: Path) -> str:
    return dt.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M")


def free_path(p: Path) -> Path:
    """*p*, or the first of name__2, name__3, ... that does not exist yet (files and folders alike)."""
    n, out = 1, p
    while out.exists():
        n += 1
        out = p.with_name(f"{p.stem}__{n}{p.suffix}")
    return out


@dataclass
class Layout:
    """Every path of one recording along the chain, for one batch's naming."""
    assay: str
    batch: str
    code: str
    phase: str
    clip: str | None
    raw_root: Path
    stage_root: Path
    codes: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.sleap = self.raw_root / self.batch / self.assay / "SLEAP"
        self.dirs = {"geom": self.sleap / "geom", "animal": self.sleap / "animal"}
        self.fmt_dir = self.sleap / "formatted" / self.phase if self.phase else self.sleap / "formatted"
        self.fmt_old = self.sleap / "formatted" / "old"
        tail = [self.batch] + ([self.phase] if self.phase else [])
        self.stage_dir = self.stage_root.joinpath(self.assay, *tail)
        self.staged = self.stage_dir / (f"{self.code}_{self.phase}.csv" if self.phase else f"{self.code}.csv")
        self.stage_arch = self.stage_root.joinpath(self.assay, f"_superseded_{TODAY}", *tail)
        who = re.escape(self.code) + (f"|{re.escape(self.clip)}" if self.clip else "")
        if self.phase:
            self.pattern = re.compile(rf"(?<![A-Za-z0-9]){re.escape(self.code)}(?:_GX\d+)?_{self.phase}(?![A-Za-z0-9])")
        else:
            self.pattern = re.compile(rf"(?<![A-Za-z0-9])(?:{who})(?![A-Za-z0-9])")

    @property
    def label(self) -> str:
        return f"{self.batch} {self.assay} {self.code}" + (f" {self.phase}" if self.phase else "")

    @property
    def by_clip(self) -> bool:
        """Open field batch 1: the locs and merged files are named by clip, not by export."""
        return self.assay == "OFT" and self.batch == "B1"

    @property
    def tag(self) -> str:
        """What a file name of this recording must contain: <CODE>_<PH>, <CODE> or (OFT B1) <CODE>_<clip>."""
        return self.code + (f"_{self.phase}" if self.phase else "") + (f"_{self.clip}" if self.by_clip else "")

    def mine(self, folder: Path, suffix: str = "") -> list[Path]:
        """Files of this recording directly in *folder* (no subfolders)."""
        if not folder.is_dir():
            return []
        return sorted(p for p in folder.iterdir() if p.is_file() and p.name.endswith(suffix) and self.pattern.search(p.name))

    def candidates(self, kind: str) -> list[Path]:
        return [p for p in self.mine(self.dirs[kind], ".analysis.h5") if "_old" not in p.name]

    def chain_locs(self, kind: str) -> list[Path]:
        return [p for p in self.mine(self.dirs[kind], "_locs.csv") if "_old" not in p.name]

    def locs_path(self, kind: str, h5_name: str) -> Path:
        if self.by_clip:                                    # named by clip, not by the h5
            return self.dirs[kind] / f"E9_B1_OFT_{self.clip}.analysis_locs.csv"
        return self.dirs[kind] / (h5_name[: -len(".h5")] + "_locs.csv")

    def install_name(self, h5: Path) -> str:
        """The name an export is installed under: its own if the tool finds it again by that name, else one with the tag."""
        name = h5.name
        if name.endswith(".analysis.h5") and "_old" not in name and self.pattern.search(name):
            return name
        base = name[: -len(".analysis.h5")] if name.endswith(".analysis.h5") else re.sub(r"\.h5$", "", name)
        return f"{base.replace('_old', '-old')}_{self.tag}.analysis.h5"

    def archived_exports(self) -> list[Path]:
        """Exports and locs files of this recording the tool archived earlier (a re-track was installed before; open field
        batch 1 keeps its old exports in predictions/, so there only the archived locs files show it)."""
        out = []
        for d in self.dirs.values():
            for a in sorted(d.glob("_superseded_*")) if d.is_dir() else []:
                out += self.mine(a, ".analysis.h5") + self.mine(a, "_locs.csv")
        return out

    def strays(self) -> list[Path]:
        """Exports of this recording anywhere under SLEAP/ outside geom/, animal/ and formatted/ (SLEAP's export dialog
        saves next to the labels file: SLEAP/projects, SLEAP/predictions ...), archives left out."""
        out = []
        for p in sorted(self.sleap.rglob("*.analysis.h5")):
            parts = p.relative_to(self.sleap).parts
            if parts[0] in ("geom", "animal", "formatted") or any(x == "old" or x.startswith("_superseded_") for x in parts[:-1]):
                continue
            if self.pattern.search(p.name):
                out.append(p)
        return out

    def formatted_existing(self) -> list[Path]:
        if not self.fmt_dir.is_dir():
            return []
        return [p for p in sorted(self.fmt_dir.glob("*.csv")) if resolve_code(p.stem, self.codes, self.code) == self.code]

    def formatted_target(self, g: str) -> Path:
        existing = self.formatted_existing()
        default = self.fmt_dir / default_formatted(self.assay, self.batch, self.code, self.phase, g)
        if default in existing or not existing:
            return default
        return existing[0]                                  # the one staging takes today: keep its name

    @property
    def fmt_archive(self) -> Path:
        """SocP keeps replaced phase files in formatted/old/<PH>_old/ where that folder exists (B1, B2)."""
        if self.phase:
            d = self.fmt_old / f"{self.phase}_old"
            if d.is_dir():
                return d
        return self.fmt_old

    def flat_copies(self) -> list[Path]:
        """SocP B3/B5: copies of the recording's formatted file directly in formatted/, beside formatted/<PH>/."""
        return self.mine(self.sleap / "formatted", ".csv") if self.phase else []

    def merged_existing(self) -> list[Path]:
        return [p for p in self.mine(self.sleap, ".csv")]

    def merged_target(self, g: str) -> Path:
        """Batches that name the merged csv after the geom locs follow the new export; the others keep the name."""
        existing = self.merged_existing()
        name, follows_stem = default_merged(self.assay, self.batch, self.code, self.phase, g)
        default = self.sleap / name
        if default in existing or follows_stem or not existing:
            return default
        named = [p for p in existing if ".merged_locs" in p.name]
        return (named or existing)[0]

    def side_files(self) -> tuple[Path, Path] | None:
        if self.assay != "SocP":
            return None
        side = "S1" if self.phase == "HAB" else self.phase          # 08: HAB is staged with the S1 side file
        return self.raw_root / self.batch / "SocP" / f"novelLoc{side}.txt", self.stage_dir / f"novelLoc{self.phase}.txt"


def side_value(path: Path, code: str) -> str | None:
    if not path.exists():
        return None
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("Code") or "").strip() == code:
                return (r.get("NovelLoc") or "").strip()
    return None


# --------------------------------------------------------------------------- #
# reading tracking files for the checks
# --------------------------------------------------------------------------- #

def read_dlc(path: Path) -> tuple[list[str], np.ndarray]:
    """Bodypart name per column and the numbers of a formatted (DLCAnalyzer) csv."""
    with open(path, newline="", encoding="utf-8") as fh:
        fh.readline()
        parts = fh.readline().rstrip("\r\n").split(",")
    return parts, pd.read_csv(path, skiprows=3, header=None).to_numpy(float)


def point(parts: list[str], data: np.ndarray, name: str) -> np.ndarray:
    idx = [i for i, p in enumerate(parts) if p == name]
    if len(idx) < 2:
        raise KeyError(name)
    return data[:, idx[:2]]


def medians(parts: list[str], data: np.ndarray, names: list[str]) -> dict[str, np.ndarray]:
    return {n: np.nanmedian(point(parts, data, n), axis=0) for n in names if n in parts}


def polygon_area(pts: list[np.ndarray]) -> float:
    x, y = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
    return float(abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2)


def px_per_cm(spec: Spec, parts: list[str], data: np.ndarray) -> float:
    m = medians(parts, data, ["tl", "tr", "br", "bl"])
    if spec.arena_cm is None:                               # EPM: tl-bl is the 60 cm span along one arm edge
        return float(np.hypot(*(m["tl"] - m["bl"])) / 60)
    w, h = spec.arena_cm
    return float(np.sqrt(polygon_area([m["tl"], m["tr"], m["br"], m["bl"]]) / (w * h)))


def body_cm(path: Path, spec: Spec) -> np.ndarray:
    """The body-centre track of a staged file in cm, at the file's own scale (as 31_light_checks.R calibrates it).

    Reads only the body centre and the corners."""
    with open(path, newline="", encoding="utf-8") as fh:
        fh.readline()
        parts = fh.readline().rstrip("\r\n").split(",")
    need = ["bodycentre"] + (["tl", "bl"] if spec.arena_cm is None else ["tl", "tr", "br", "bl"])
    cols = sorted(i for n in need for i in [j for j, p in enumerate(parts) if p == n][:2])
    data = pd.read_csv(path, skiprows=3, header=None, usecols=cols).to_numpy(float)
    sub = [parts[i] for i in cols]
    return point(sub, data, "bodycentre") / px_per_cm(spec, sub, data)


def median_gap(a: np.ndarray, b: np.ndarray, lag: int) -> float:
    """Median distance between a(t) and b(t + lag) over the frames both have (inf below DUPLICATE_MIN_OVERLAP)."""
    s, e = max(0, -lag), min(len(a), len(b) - lag)
    if e - s < DUPLICATE_MIN_OVERLAP:
        return float("inf")
    d = np.hypot(*(a[s:e] - b[s + lag:e + lag]).T)
    d = d[np.isfinite(d)]
    return float(np.median(d)) if d.size else float("inf")


def closest_offset(a: np.ndarray, b: np.ndarray) -> tuple[float, int]:
    """Smallest median distance between tracks *a* and *b* (a(t) against b(t + offset), offsets up to DUPLICATE_LAG
    frames), with the limits of 31_light_checks.R. The offset is found by FFT, as there, but as the smallest mean
    squared distance over the frames both have (correlation plus the sums of squares of the overlap) rather than the
    highest correlation, which can sit far off the true offset for slowly wandering tracks; the median distance is
    then taken within 3 frames of it, as there. Missing points take the column median."""
    def fill(m: np.ndarray) -> np.ndarray:
        m = np.array(m, dtype=float)
        for k in range(2):
            v = m[:, k]
            ok = np.isfinite(v)
            v[~ok] = np.median(v[ok]) if ok.any() else 0.0
        return m
    fa, fb = fill(a), fill(b)
    centre = fa.mean(axis=0)                                # one shift for both keeps the distances
    fa, fb = fa - centre, fb - centre
    na, nb = len(fa), len(fb)
    n = 1 << (na + nb).bit_length()
    cc = sum(np.fft.irfft(np.fft.rfft(fb[:, k], n) * np.conj(np.fft.rfft(fa[:, k], n)), n) for k in range(2))
    lags = np.arange(-DUPLICATE_LAG, DUPLICATE_LAG + 1)        # cc[lag] = sum over t of a(t) b(t + lag)
    sa = np.concatenate([[0.0], np.cumsum((fa ** 2).sum(axis=1))])
    sb = np.concatenate([[0.0], np.cumsum((fb ** 2).sum(axis=1))])
    s = np.maximum(0, -lags)
    e = np.minimum(na, nb - lags)
    overlap = e - s
    ok = overlap >= DUPLICATE_MIN_OVERLAP
    msd = np.full(len(lags), np.inf)
    sl, el, ll = s[ok], e[ok], lags[ok]
    msd[ok] = (sa[el] - sa[sl] + sb[el + ll] - sb[sl + ll] - 2 * cc[ll % n]) / overlap[ok]
    if not ok.any():
        return float("inf"), 0
    best = int(lags[np.argmin(msd)])
    near = [best + d for d in range(-3, 4)]
    gaps = [median_gap(a, b, lag) for lag in near]
    i = int(np.argmin(gaps))
    return gaps[i], near[i]


# --------------------------------------------------------------------------- #
# staging (08_stage_all_batches.R's copy, EPM with the neck alias)
# --------------------------------------------------------------------------- #

NECK_R = r"""
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
for (e in parse(file = args[1])) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("copy_with_neck"))) eval(e)
}
if (!exists("copy_with_neck")) stop("copy_with_neck not found in ", args[1])
res <- copy_with_neck(args[2], args[3])
cat(res, "\n")
if (!grepl("^copied", res)) quit(status = 1)
"""


def neck_copy_python(src: Path, dst: Path) -> str:
    """copy_with_neck in Python: the same columns, numbers at 15 significant digits (R's fwrite rounds the last digit
    its own way, so this is equal to 1e-11 px but not byte-identical to 08's file)."""
    with open(src, newline="", encoding="utf-8") as fh:
        head = [fh.readline().rstrip("\r\n").split(",") for _ in range(3)]
    if "neck" in head[1]:
        shutil.copyfile(src, dst)
        return "copied (had neck)"
    idx = [i for i, p in enumerate(head[1]) if p == "spine1"]
    if len(idx) != 3 or [head[2][i].lower() for i in idx] != ["x", "y", "likelihood"]:
        raise SystemExit(f"{src.name}: no usable spine1 triplet for the neck alias")
    body = pd.read_csv(src, skiprows=3, header=None)
    body = pd.concat([body, body.iloc[:, idx].set_axis(range(body.shape[1], body.shape[1] + 3), axis=1)], axis=1)
    with open(dst, "w", newline="") as fh:
        fh.write(",".join(head[0] + ["column_neck_x", "column_neck_y", "column_neck_lh"]) + "\r\n")
        fh.write(",".join(head[1] + ["neck"] * 3) + "\r\n")
        fh.write(",".join(head[2] + ["x", "y", "likelihood"]) + "\r\n")
        body.to_csv(fh, index=False, header=False, float_format="%.15g", lineterminator="\r\n")
    return "copied (+neck, Python)"


#: The EPM staging copy; None = 08's copy_with_neck through Rscript (Python in a plan if Rscript is missing). Tests
#: replace it.
NECK_COPIER = None


def neck_copy(src: Path, dst: Path, rscript: Path, stage_script: Path, notes: list[str],
              apply: bool = False, allow_python: bool = False) -> str:
    if NECK_COPIER is not None:
        return NECK_COPIER(src, dst)
    if not rscript.exists() or not stage_script.exists():
        missing = rscript if not rscript.exists() else stage_script
        if apply and not allow_python:
            raise SystemExit(f"{missing} not found: the EPM staged copy is written by 08's copy_with_neck through Rscript. "
                             "Point --rscript at Rscript.exe, or add --python-neck-copy to write it in Python (equal to "
                             "1e-11 px, not byte-identical to what 08 writes)")
        notes.append(f"{missing} not found: staged copy written in Python (equal to 1e-11 px, not byte-identical to "
                     "what 08 writes); --apply needs Rscript or --python-neck-copy")
        return neck_copy_python(src, dst)
    prog = dst.parent / "copy_with_neck.R"
    prog.write_text(NECK_R, encoding="utf-8")
    res = subprocess.run([str(rscript), str(prog), stage_script.as_posix(), src.as_posix(), dst.as_posix()],
                         capture_output=True, text=True)
    if res.returncode != 0 or not dst.exists():
        raise SystemExit(f"copy_with_neck failed: {res.stdout[-500:]} {res.stderr[-800:]}")
    return res.stdout.strip()


# --------------------------------------------------------------------------- #
# EPM, OFT, SocP: the run
# --------------------------------------------------------------------------- #

@dataclass
class Action:
    verb: str               # "install", "replace", "unchanged" or "move"
    dst: Path
    src: Path | None = None
    archive: Path | None = None
    target: Path | None = None   # where the old file goes (set by check_actions)

    def line(self) -> str:
        if self.verb == "move":
            return f"  move        {self.dst}   (-> {self.archive})"
        if self.verb == "replace":
            return f"  replace     {self.dst}   (old -> {self.archive})"
        return f"  {self.verb:<11} {self.dst}"


def plan_file(new: Path, dst: Path, archive: Path) -> Action:
    if dst.exists():
        return Action("unchanged" if sha(dst) == sha(new) else "replace", dst, new, archive)
    return Action("install", dst, new)


def check_actions(actions: list[Action]) -> None:
    """Decide every archive name before anything moves; an archived file is never overwritten (a second one the same
    day is numbered: name__2.csv)."""
    taken: set[Path] = set()
    for a in actions:
        if a.verb in ("replace", "move"):
            target, n = a.archive / a.dst.name, 1
            while target.exists() or target in taken:
                n += 1
                target = a.archive / f"{a.dst.stem}__{n}{a.dst.suffix}"
            a.target = target
            taken.add(target)


def do_actions(actions: list[Action]) -> None:
    for a in actions:
        if a.verb in ("replace", "move"):
            a.archive.mkdir(parents=True, exist_ok=True)
            if a.target is None or a.target.exists():
                raise SystemExit(f"no free archive name for {a.dst}")
            shutil.move(str(a.dst), str(a.target))
        if a.verb in ("install", "replace"):
            a.dst.parent.mkdir(parents=True, exist_ok=True)
            tmp = a.dst.with_name(a.dst.name + ".installing")
            shutil.copyfile(a.src, tmp)
            if sha(tmp) != sha(a.src):
                tmp.unlink()
                raise SystemExit(f"hash mismatch writing {a.dst}")
            os.replace(tmp, a.dst)


def put_file(src: Path, dst: Path) -> None:
    """Copy *src* over *dst* through a temporary name, checking the bytes."""
    tmp = dst.with_name(dst.name + ".installing")
    shutil.copyfile(src, tmp)
    if sha(tmp) != sha(src):
        tmp.unlink()
        raise SystemExit(f"hash mismatch writing {dst}")
    os.replace(tmp, dst)


def same_table(df: pd.DataFrame, path: Path) -> bool:
    """*df* equals the csv at *path* to 1e-9 (columns, shape and missing values alike)."""
    if not path.exists():
        return False
    old = pd.read_csv(path)
    if list(old.columns) != list(df.columns) or old.shape != df.shape:
        return False
    a, b = df.to_numpy(float), old.to_numpy(float)
    if not np.array_equal(np.isnan(a), np.isnan(b)):
        return False
    d = np.abs(a - b)
    return bool(np.nanmax(d, initial=0.0) <= 1e-9)


def count_rows(path: Path) -> int:
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


def h5_info(path: Path) -> tuple[int, str]:
    with h5py.File(path, "r") as f:
        frames = int(f["tracks"].shape[-1])
        vp = f["video_path"][()] if "video_path" in f else b""
    return frames, vp.decode(errors="replace") if isinstance(vp, bytes) else str(vp)


@dataclass
class Exports:
    """What each kind (geom, animal) brings: a new export, the export in place, or (geom) the locs reused."""
    new: dict = field(default_factory=dict)        # kind -> h5 path
    current: dict = field(default_factory=dict)    # kind -> h5 paths in place that rebuild their locs file
    stale: dict = field(default_factory=dict)      # kind -> h5 paths in place that do not rebuild it
    reuse: dict = field(default_factory=dict)      # kind -> locs path
    tables: dict = field(default_factory=dict)     # h5 path -> extracted locs table
    errors: dict = field(default_factory=dict)     # h5 path -> why it could not be read

    def inplace(self, kind: str) -> list[Path]:
        return self.current.get(kind) or self.stale.get(kind, [])


def find_exports(lay: Layout, args: argparse.Namespace) -> Exports:
    """Sort the exports of the recording in geom/ and animal/.

    In place: it rebuilds the locs file of the chain (current), or it has a locs file of its own that it does not
    rebuild, or it is older than the chain's locs files (stale: rebuilt from only with --rebuild). New: no locs file of
    its own and newer than the chain's locs files (open field B1 names its locs by clip, so there an export never has
    its own)."""
    ex = Exports()
    explicit = {"geom": args.geom_h5, "animal": args.animal_h5}
    for kind in ("geom", "animal"):
        if kind == "geom" and args.reuse_geom:
            locs = lay.chain_locs(kind)
            if len(locs) != 1:
                raise SystemExit(f"--reuse-geom: expected one locs file of {lay.code} in {lay.dirs['geom']}, found {[p.name for p in locs]}")
            ex.reuse[kind] = locs[0]
            continue
        cands = lay.candidates(kind)
        if explicit[kind] is not None:
            p = explicit[kind].resolve()
            if not p.exists():
                raise SystemExit(f"--{kind}-h5 {p} does not exist")
            ex.new[kind] = p
            cands = [c for c in cands if c.resolve() != p]
        newest = max((q.stat().st_mtime for q in lay.chain_locs(kind)), default=None)
        cur, new, stale = [], [], []
        for c in cands:
            own = lay.locs_path(kind, c.name)
            df = None
            try:
                df, _ = sh.extract_locs(c, kind, lay.assay, batch=lay.batch)
                ex.tables[c] = df
            except (ValueError, OSError, KeyError) as exc:
                ex.errors[c] = str(exc)
            if df is not None and same_table(df, own):
                cur.append(c)
            elif (lay.by_clip or not own.exists()) and (newest is None or c.stat().st_mtime > newest):
                new.append(c)
            else:
                stale.append(c)
        ex.current[kind], ex.stale[kind] = cur, stale
        if kind in ex.new:
            continue
        if len(new) > 1:
            raise SystemExit(f"{len(new)} new {kind} exports for {lay.label} in {lay.dirs[kind]}: {[p.name for p in new]}; "
                             f"name the one to use with --{kind}-h5 (the others are then moved to _superseded_{TODAY}/)")
        if new:
            ex.new[kind] = new[0]
    return ex


def identity_problems(lay: Layout, kind: str, h5: Path, video: str) -> list[str]:
    """The export's file name and stored video_path must name the recording (OFT B1: the clip) and no other animal."""
    vname = re.split(r"[\\/]", video)[-1] if video else ""
    if lay.by_clip:
        what = f"{lay.code} or {lay.clip}"
        name_ok = bool(token(lay.code).search(h5.name) or token(lay.clip).search(h5.name))
        video_ok, vwhat = bool(token(lay.clip).search(vname)), lay.clip
    elif lay.phase:
        what = vwhat = f"{lay.code}_{lay.phase}"
        name_ok, video_ok = bool(lay.pattern.search(h5.name)), bool(lay.pattern.search(vname))
    else:
        what = vwhat = lay.code
        name_ok, video_ok = bool(token(lay.code).search(h5.name)), bool(token(lay.code).search(vname))
    out = []
    if not name_ok:
        out.append(f"the {kind} export's file name {h5.name} does not name {what}")
    if not video_ok:
        out.append(f"the {kind} export's video_path ({video or 'none stored'}) does not name {vwhat}")
    others = [c for c in lay.codes if c != lay.code and (token(c).search(h5.name) or token(c).search(vname))]
    if others:
        out.append(f"the {kind} export names {', '.join(others)} (file {h5.name}, video {vname or '-'})")
    return out


def check_clip(lay: Layout) -> tuple[bool, str]:
    """OFT B1: the merged csv of --clip must format to the animal's formatted file in place, or the clip is another
    animal's (the tool would rewrite that animal's locs and merged files)."""
    fmt = lay.formatted_existing()
    merged = [p for p in lay.merged_existing() if token(lay.clip).search(p.name) and ".merged_locs" in p.name]
    if not fmt:
        return True, f"{lay.code} has no formatted file in {lay.fmt_dir}: --clip {lay.clip} cannot be checked against it"
    if not merged:
        return False, (f"no merged csv of clip {lay.clip} in {lay.sleap}: cannot confirm that {lay.clip} is {lay.code}'s "
                       "clip (the pair sheet B1/OFT/E9_SIS_B1_OF_19Nov2022.txt gives it; --accept-clip overrides)")
    with tempfile.TemporaryDirectory(prefix="retrack_clip_") as tmp:
        out = Path(tmp) / fmt[0].name
        sd.format_merged(merged[0], out, FORMAT_RENAME.get((lay.assay, lay.batch)), FORMAT_DIGITS.get((lay.assay, lay.batch)))
        if out.read_bytes() == fmt[0].read_bytes():
            return True, f"{merged[0].name} formats to {fmt[0].name} byte for byte: {lay.clip} is {lay.code}'s clip"
        pa, a = read_dlc(out)
        pb, b = read_dlc(fmt[0])
        if pa == pb and a.shape == b.shape and np.allclose(a, b, rtol=1e-9, atol=1e-6, equal_nan=True):
            return True, f"{merged[0].name} formats to {fmt[0].name} (to 10 significant digits): {lay.clip} is {lay.code}'s clip"
    return False, (f"{merged[0].name} does not format to {lay.code}'s {fmt[0].name}: --clip {lay.clip} is another "
                   "animal's clip, or the chain was edited (the pair sheet B1/OFT/E9_SIS_B1_OF_19Nov2022.txt gives the "
                   "clip; --accept-clip overrides)")


def example_name(lay: Layout, kind: str, current: list[Path]) -> str:
    if current:
        name = current[0].name
        if re.search(r"\.\d{3}_", name):
            return re.sub(r"\.\d{3}_", ".<NNN>_", name, count=1)
        return f"<labels file>.<NNN>_{name}"                # never the name in place itself
    if lay.by_clip:
        return f"<labels file>.<NNN>_{lay.code}_{lay.clip}.analysis.h5"
    return f"<labels file>.<NNN>_<video name with _{lay.code}{'_' + lay.phase if lay.phase else ''}>.analysis.h5"


def print_strays(lay: Layout) -> None:
    newest = max((q.stat().st_mtime for k in lay.dirs for q in lay.chain_locs(k)), default=None)
    for p in lay.strays():
        if newest is None or p.stat().st_mtime > newest:
            print(f"  !! export not in geom/ or animal/ (the tool does not read it): {p}  ({when(p)}, newer than the "
                  "chain's locs files); if it is the new export, move it to geom/ or animal/, or pass --geom-h5/--animal-h5")
        else:
            print(f"  elsewhere {p}  ({when(p)}, older than the chain's locs files: an earlier export, not read; a new "
                  "export goes to geom/ or animal/, or is passed with --geom-h5/--animal-h5)")


def print_waiting(lay: Layout, spec: Spec, ex: Exports, missing: list[str], args: argparse.Namespace,
                  again: bool = False) -> None:
    frames = round(spec.expected_s * VIDEO_FPS)
    must = f"_{lay.code}_{lay.phase}" if lay.phase else f"_{lay.code}" + (f" (or {lay.clip})" if lay.clip else "")
    if again:
        print("  if you re-track again:")
    else:
        print(f"  WAITING FOR THE NEW EXPORT: no new {' and '.join(missing)} export of {lay.label} found.")
    for kind in missing:
        what = "arena (geom)" if kind == "geom" else "animal"
        print(f"    put the {what} export in {lay.dirs[kind]}")
        print(f"      e.g. {example_name(lay, kind, ex.inplace(kind))}")
    print(f"      the name must contain {must}, end in .analysis.h5 and not contain _old; or keep it anywhere and pass "
          f"--geom-h5 / --animal-h5 PATH (it is then copied there)")
    for kind in missing:
        for p in ex.inplace(kind):
            print(f"      do not export under the name of the {kind} export in place ({p.name}): SLEAP would overwrite it "
                  "before the tool can archive it; keep SLEAP's new default name")
    names = {k: sh.H5_NODES[(lay.assay, k)] for k in missing}
    for k, n in names.items():
        print(f"      {k} export: one track; nodes in this order: {', '.join(n)}")
    mates = [p for p in sorted(lay.stage_dir.glob("*.csv")) if p != lay.staged] if lay.stage_dir.is_dir() else []
    mate = f"; {mates[0].stem} of this batch has {count_rows(mates[0]) - 3}" if mates else ""
    print(f"      length: {spec.expected_s:g} s = {frames} frames at 29.97/s (the light checks allow {LENGTH_LIMIT_S:g} s either way){mate}")
    if "geom" in ex.reuse:
        print(f"      --reuse-geom: the animal export must have exactly {count_rows(ex.reuse['geom']) - 1} frames, as the "
              "corners in place; otherwise track the arena on the new clip too")
    if "geom" in missing and "animal" not in missing:
        print("      if the arena was not re-tracked, --reuse-geom keeps the corners in place (same frame count needed)")
    if lay.assay == "OFT" and "geom" in missing and not args.reuse_geom:
        print("      open field: the corners are often not re-tracked; --reuse-geom keeps the current geom locs file")
    print_strays(lay)


def print_current(lay: Layout, ex: Exports) -> None:
    print("  current chain (stays until a new export is installed):")
    for kind in ("geom", "animal"):
        if kind in ex.reuse:
            print(f"    {kind:<9} reused   {ex.reuse[kind]}")
            continue
        for p in ex.current.get(kind, []) + ex.stale.get(kind, []):
            fr, _ = h5_info(p)
            how = "rebuilds its locs file" if p in ex.current.get(kind, []) else "in place, does NOT rebuild its locs file"
            print(f"    {kind:<9} h5       {p}   ({fr} frames, {when(p)}; {how})")
        if len(ex.inplace(kind)) > 1:
            print(f"    {kind:<9}          {len(ex.inplace(kind))} exports of this recording are in place; all of them "
                  f"move to _superseded_{TODAY}/ when the new one is installed")
        for p in lay.mine(lay.dirs[kind], ".csv"):
            print(f"    {kind:<9} locs     {p}")
    for p in lay.merged_existing():
        print(f"    merged             {p}")
    for p in lay.formatted_existing():
        print(f"    formatted          {p}")
    if lay.staged.exists():
        print(f"    staged             {lay.staged}   (md5 {md5(lay.staged)})")
    else:
        print(f"    staged             {lay.staged}   (missing)")


def scoring_setup(lay: Layout, spec: Spec, args: argparse.Namespace, config: Path, cfg_text: str,
                  stage_default: Path, output_default: Path | None) -> dict:
    repo = (args.scoring_repo or REPO).resolve()
    if spec.assay == "OFT":
        sub = text_value(cfg_text, "output_subdir") or "output_v1.2.0"
        out_dirs = [lay.stage_root / "OFT" / lay.batch / sub]
        global_dir = lay.stage_root / "OFT" / "OFT_SLEAP_output_v1.2.0"
        backup_root = lay.stage_root / f"_superseded_{TODAY}_{lay.code}"
    else:
        out_root = args.output_root or output_default
        tail = [lay.batch] + ([lay.phase] if lay.phase else [])
        out_dirs = [out_root.joinpath(spec.assay, *tail)]
        global_dir = None
        backup_root = out_root.parent / f"_superseded_{TODAY}_{lay.code}"
    tag = "_".join([spec.assay, lay.batch] + ([lay.phase] if lay.phase else []))
    text = re.sub(r"^batches:.*$", f"batches: [{lay.batch}]", cfg_text, flags=re.M)
    if lay.phase:
        text = re.sub(r"^phases:.*$", f"phases: [{lay.phase}]", text, flags=re.M)
    moves = [(stage_default.as_posix(), lay.stage_root.as_posix())]
    if output_default is not None and args.output_root:
        moves.append((output_default.as_posix(), args.output_root.as_posix()))
    text = retarget(text, moves)

    def absolute(m: re.Match) -> str:                     # relative paths resolve against the config's own folder
        val = m.group(2).strip().strip("'\"")
        if val in ("", "~") or Path(val).is_absolute():
            return m.group(0)
        return f"{m.group(1)}: {(config.parent / val).resolve().as_posix()}"
    text = re.sub(r"^(qc_decisions_file|zone_file|animal_id_code_file|metadata_dir):\s*(.+?)\s*$", absolute, text, flags=re.M)
    manifest = out_dirs[0] / "run_manifest.yaml"
    return dict(repo=repo, script=repo / "02_SLEAPanalzyer" / spec.script, text=text, out_dirs=out_dirs,
                global_dir=global_dir, backup_root=backup_root, tag=tag, manifest=manifest, marker=out_dirs[0] / MARKER)


def code_status(sc: dict, spec: Spec) -> tuple[str, str]:
    """('same' | 'changed' | 'unknown', explanation) for the scoring code against the code of the current output."""
    m = sc["manifest"]
    if not m.exists():
        return "unknown", f"no run_manifest.yaml in {m.parent}: the code version of the current output is unknown"
    text = m.read_text(encoding="utf-8", errors="replace")
    commit = re.search(r"^\s*commit:\s*([0-9a-f]{7,40})", text, re.M)
    stamp_at = re.search(r"timestamp_utc:\s*'?([^'\n]+)", text)
    if not commit:
        return "unknown", f"{m} names no commit"
    c = commit.group(1)
    paths = [f"02_SLEAPanalzyer/{spec.script}"] + SCORING_CODE
    res = subprocess.run(["git", "-C", str(sc["repo"]), "diff", "--name-only", c, "--", *paths],
                         capture_output=True, text=True)
    stamp = f"{stamp_at.group(1) if stamp_at else '?'}, commit {c[:7]}"
    if res.returncode != 0:
        return "unknown", f"current output made {stamp}; git cannot compare it ({res.stderr.strip()[:200]})"
    files = [ln for ln in res.stdout.splitlines() if ln.strip()]
    if not files:
        return "same", f"current output made {stamp}; the scoring code in {sc['repo']} is the same"
    return "changed", (f"current output made {stamp}; the scoring code in {sc['repo']} differs in {len(files)} file(s): "
                       + ", ".join(Path(f).name for f in files[:6]) + (" ..." if len(files) > 6 else ""))


#: Function returning the code status; the tests replace it.
CODE_STATUS = code_status


def manifest_md5(manifest: Path, staged: Path) -> str | None:
    """The md5 run_manifest.yaml records for the input *staged* (matched by file and folder name), or None."""
    text = manifest.read_text(encoding="utf-8", errors="replace")
    for m in re.finditer(r"^\s*-\s*path:\s*(.+?)\s*$((?:\n(?!\s*-\s*path:).*)*)", text, re.M):
        p = Path(m.group(1).strip("'\""))
        if p.name == staged.name and p.parent.name == staged.parent.name:
            hit = re.search(r"^\s*md5:\s*'?([0-9a-f]{32})", m.group(2), re.M)
            return hit.group(1) if hit else None
    return None


def read_marker(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["file"]: r["md5"] for r in csv.DictReader(fh, delimiter="\t")}


def write_marker(path: Path, name: str, digest: str, note: str) -> None:
    rows = []
    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("file") != name]
    rows.append({"file": name, "md5": digest, "written": dt.datetime.now().isoformat(timespec="seconds"), "note": note})
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, ["file", "md5", "written", "note"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def scored_state(sc: dict, lay: Layout) -> tuple[str | None, str]:
    """(md5 of the staged file the current output was scored from, where that is recorded).

    From the output's run_manifest.yaml ("" when it does not list the file), else from the tool's marker
    (open field); (None, "") when neither tells."""
    if sc["manifest"].exists():
        return manifest_md5(sc["manifest"], lay.staged) or "", str(sc["manifest"])
    rows = read_marker(sc["marker"])
    if lay.staged.name in rows:
        return rows[lay.staged.name], str(sc["marker"])
    return None, ""


def scoring_note(lay: Layout, sc: dict, scored: tuple[str | None, str], staged_md5: str | None) -> tuple[bool, str]:
    """(scoring due?, explanation) for the staged file in place."""
    digest, where = scored
    if digest is None:
        return False, (f"which staged file the current output was scored from is not recorded (no run_manifest.yaml in "
                       f"{sc['out_dirs'][0]} and no {MARKER} yet; the tool writes one when it installs a file there)")
    if staged_md5 is None:
        return False, f"no staged file in place ({lay.staged})"
    if digest == staged_md5:
        return False, f"the current output was scored from the staged file in place ({where})"
    told = f"records md5 {digest[:8]}" if digest and digest != "-" else "does not list it"
    return True, (f"scoring of the installed file is still due: {where} {told} for {lay.staged.name}, the staged file "
                  f"in place is {staged_md5[:8]}")


def previous_staged(lay: Layout, digest: str) -> Path | None:
    """The archived staged copy with md5 *digest* (the file the current output was scored from), newest archive first."""
    tail = [lay.batch] + ([lay.phase] if lay.phase else [])
    same = re.compile(rf"^{re.escape(lay.staged.stem)}(__\d+)?{re.escape(lay.staged.suffix)}$")
    for d in sorted(lay.stage_root.joinpath(lay.assay).glob("_superseded_*"), reverse=True):
        folder = d.joinpath(*tail)
        if folder.is_dir():
            for p in sorted(folder.iterdir(), reverse=True):
                if same.match(p.name) and md5(p) == digest:
                    return p
    return None


def read_keyed(path: Path, table: Table) -> dict:
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {tuple(r.get(k, "") for k in table.key): r for r in csv.DictReader(fh)}


def snapshot(out_dir: Path, spec: Spec) -> dict:
    return {t.name: read_keyed(out_dir / t.name, t) for t in spec.tables}


def compare(before: dict, after: dict, spec: Spec, code: str) -> tuple[dict, set]:
    per_table, others = {}, set()
    for t in spec.tables:
        diff = changed_cells(before.get(t.name, {}), after.get(t.name, {}), t.ignore)
        codes = sorted({k[0][:4] for k in diff})
        per_table[t.name] = (len(after.get(t.name, {})), codes, diff)
        others |= {c for c in codes if c != code}
    return per_table, others


def tool_command(argv: list[str], add: list[str] = (), drop: tuple[str, ...] = ("--apply",)) -> str:
    keep = [a for a in argv if a not in drop]
    return subprocess.list2cmdline([sys.executable, str(Path(__file__).resolve())] + keep + list(add))


def row_phase(r: dict) -> str | None:
    """The phase of a retrack or leave-out row: its phase column, else its trial (habituation, phase 1, phase 2)."""
    for col in ("phase", "trial"):
        if col in r:
            v = (r.get(col) or "").strip()
            return TRIAL_PHASE.get(v.lower(), v)
    return None


def print_followups(lay: Layout, spec: Spec, args: argparse.Namespace, new_md5: str | None,
                    qc_file: Path | None = None) -> None:
    print("  then, by hand:")
    ph = lay.phase or lay.assay
    cur = md5(lay.staged) if lay.staged.exists() else None
    rows, related = [], []
    if args.decisions.exists():
        with open(args.decisions, newline="", encoding="utf-8") as fh:
            for i, r in enumerate(csv.DictReader(fh), start=2):
                if r.get("assay") != lay.assay or r.get("batch") != lay.batch:
                    continue
                if r.get("phase") == ph and r.get("Code") == lay.code:
                    rows.append((i, r))
                elif (r.get("Code") == lay.code and "duplicate" in (r.get("flags") or "")) or \
                        (r.get("Code") != lay.code and lay.code in (r.get("finding") or "")):
                    related.append((i, r))      # e.g. the other copy of a duplicated trial
    if rows:
        print(f"    light_check_decisions.csv ({args.decisions}) has {len(rows)} row(s) for {lay.assay} {ph} {lay.batch} {lay.code}, keyed to the staged md5:")
        for i, r in rows:
            k = r.get("tracking_md5", "")
            state = "= staged file now" if k == cur else ("= new staged file" if k == new_md5 else "matches neither")
            print(f"      line {i}: {r.get('flags')} -> {r.get('outcome')}  (md5 {k[:8]}, {state})")
        print("      after the new file is staged and 31 has run, these show as 'tracking changed': review, then update or drop them")
    if related:
        print("    related rows, whose flag may no longer be raised after the re-track (review them too):")
        for i, r in related:
            print(f"      line {i}: {r.get('assay')} {r.get('phase')} {r.get('batch')} {r.get('Code')}: {r.get('flags')} -> {r.get('outcome')}")
    if qc_file is not None and qc_file.exists() and lay.phase:
        with open(qc_file, newline="", encoding="utf-8-sig") as fh:
            hits = [(i, r) for i, r in enumerate(csv.DictReader(fh), start=2)
                    if (r.get("file") or "").strip() == f"{lay.code}_{lay.phase}"]
        if hits:
            print(f"    {qc_file.name} ({qc_file}) has {len(hits)} row(s) for {lay.code}_{lay.phase}; their times refer to the "
                  "old clip: check them on the new one, then keep, change or drop them (the scoring applies them):")
            for i, r in hits:
                print(f"      line {i}: {r.get('target')} {r.get('start_s')}-{r.get('end_s')} s -> {r.get('decision')}: "
                      f"{(r.get('reason') or '')[:90]}")
    qc = args.qc_root
    for name in ("retrack_list.csv", "leave_out_list.csv"):
        p = qc / name
        if not p.exists():
            continue
        with open(p, newline="", encoding="utf-8-sig") as fh:
            hits = [r for r in csv.DictReader(fh) if (r.get("code") or "").strip() == lay.code
                    and (r.get("test") or "").strip() == lay.assay
                    and not (lay.phase and row_phase(r) in spec.phases and row_phase(r) != lay.phase)]
        for r in hits:
            what = r.get("problem") or r.get("reason") or ""
            print(f"    {name}: {r.get('animal')} {lay.code} {lay.assay} {r.get('trial')}: {what[:110]}"
                  + (f" [status: {r.get('status')}]" if r.get("status") else "") + " -> close or update it after the re-track")
    print(f"    video QC (32_video_qc_manifest.R) wants one clip per code and trial in the batch's cropped/ folder: "
          f"move a replaced clip to cropped/old/ (the tool does not touch videos)")
    only = f"{lay.assay} {ph} {lay.batch} {lay.code}"
    batch_sel = f"^{lay.assay} [A-Za-z0-9]+ {lay.batch} "
    print(f"    light check (writes to a temporary folder). 31 compares for the same trial twice only among the recordings "
          f"it selects, so a one-recording selection ('{only}') does not check duplicates; this one takes all of "
          f"{lay.assay} {lay.batch}" + (", every phase" if spec.phases else "") + ":")
    print(f"      $env:EXP9_LIGHT_ONLY = '{batch_sel}'; & '{args.rscript}' '{LIGHT_CHECKS}'; Remove-Item Env:EXP9_LIGHT_ONLY")
    print(f"    after all re-tracks, a release (runs 05, 09-13 and the QC; about 30 min more without -SkipQC):")
    print(f"      powershell -NoProfile -ExecutionPolicy Bypass -File '{RELEASE_SCRIPT}' -ReleaseName {TODAY}_exp9_all-batches_canonical-retrack")
    print("    then the SISanalyzer ingest of that release")


def print_scoring_plan(lay: Layout, spec: Spec, sc: dict, status: tuple[str, str], args: argparse.Namespace,
                       argv: list[str], config: Path, scored_text: str = "") -> None:
    cfg_note = f"batches: [{lay.batch}]" + (f", phases: [{lay.phase}]" if lay.phase else "")
    print(f"  scoring   {sc['script'].name} on {lay.batch}" + (f" {lay.phase}" if lay.phase else "")
          + f" only: temporary copy of {config} with {cfg_note}")
    for d in sc["out_dirs"]:
        print(f"            output {d}  (backed up to {sc['backup_root'] / (sc['tag'] + '_before')} first)")
    if sc["global_dir"] is not None:
        print(f"            the one-batch run rewrites {sc['global_dir']} with {lay.batch} only: it is backed up and put back afterwards")
    if spec.assay == "SocP":
        print("            landmark repair, head restore and the on-cup decisions (qc_decisions_socp.csv) run inside this scoring")
    if scored_text:
        print(f"            {scored_text}")
    print(f"            {status[1]}")
    if status[0] == "same":
        print("            -> compared against the backup: only this animal may change")
        return
    if args.baseline:
        print("            -> --baseline: the batch is scored once on the files the current output was scored from, then "
              "again with the new file; only this animal may change between the two")
        return
    print("            -> every animal of the batch would change. Choose one, then add it to the --apply command:")
    print(f"               --baseline (scores the batch twice):  {tool_command(argv, ['--apply', '--baseline'])}")
    m = sc["manifest"]
    commit = None
    if m.exists():
        hit = re.search(r"^\s*commit:\s*([0-9a-f]{7,40})", m.read_text(encoding="utf-8", errors="replace"), re.M)
        commit = hit.group(1)[:7] if hit else None
    if commit:
        wt = REPO.parent / f"SLEAPanalyzer_scoring_{commit}"
        print(f"               or score with the code of {commit}:  git -C '{REPO}' worktree add '{wt}' {commit}")
        print(f"                  then:  {tool_command(argv, ['--apply', '--scoring-repo', str(wt)])}")
    print("            (without one of them, --apply stops before changing anything; --skip-analysis installs without scoring)")


def run_assay(args: argparse.Namespace, argv: list[str]) -> int:
    spec = SPECS[args.assay]
    if args.static_landmark:
        raise SystemExit(f"--static-landmark is a NOR option: the {args.assay} files interpolate every geom point")
    phase = args.phase or ""
    if spec.phases and phase not in spec.phases:
        raise SystemExit(f"{args.assay} needs --phase {' / '.join(spec.phases)}")
    if not spec.phases and phase not in ("", args.assay):
        raise SystemExit(f"{args.assay} has one trial per animal; leave out --phase")
    phase = phase if spec.phases else ""
    if args.assay == "OFT" and args.batch == "B1" and not args.clip:
        raise SystemExit("open field batch 1 names its tracking files by GoPro clip: pass --clip (X1H0 = GX010804; "
                         "the pair sheet B1/OFT/E9_SIS_B1_OF_19Nov2022.txt or the h5 video_path give the clip)")
    if args.reuse_geom and (args.geom_h5 is not None):
        raise SystemExit("--reuse-geom and --geom-h5 exclude each other")

    if args.config is not None:
        config = Path(os.path.abspath(args.config))
        if not config.exists():
            raise SystemExit(f"--config {config} does not exist")
    else:
        config = (args.scoring_repo or REPO) / CONFIG_DIR.relative_to(REPO) / spec.configs[phase]
        if not config.exists():                             # a scoring checkout older than the config
            config = CONFIG_DIR / spec.configs[phase]
    cfg_text = config.read_text(encoding="utf-8")
    if args.assay == "OFT":
        stage_default, output_default = Path(text_value(cfg_text, "behavior_root", config.name)).parent, None
        if args.output_root:
            raise SystemExit("--output-root: the open-field output lives in the staging tree; use --stage-root")
    else:
        stage_default = Path(text_value(cfg_text, "input_dir", config.name)).parent
        output_default = Path(text_value(cfg_text, "output_dir", config.name)).parent
    qc_value = text_value(cfg_text, "qc_decisions_file", config.name, required=False)
    qc_file = (config.parent / qc_value).resolve() if qc_value and qc_value != "~" else None
    stage_root = args.stage_root or stage_default
    lay = Layout(args.assay, args.batch, args.code, phase, args.clip, args.raw_root, stage_root, load_codes(stage_root))
    mode = "APPLY" if args.apply else "PLAN (nothing is written; add --apply)"
    print(f"{lay.label}: {mode}")
    if not lay.sleap.is_dir():
        raise SystemExit(f"{lay.sleap} does not exist")
    if args.assay == "OFT" and args.batch == "B2":
        print("  note      the 2023 B2 animal locs files left a column empty when its first or last frame was empty "
              "(4 animals, columns the score does not use); a new export is filled like every other batch")

    problems: list[str] = []
    notes: list[str] = []
    if lay.by_clip:
        ok, why = check_clip(lay)
        print(f"  clip      {why}" if ok else f"  !! clip   {why}")
        if not ok and not args.accept_clip:
            problems.append(why)
        elif not ok:
            notes.append(f"clip check overridden (--accept-clip): {why}")

    ex = find_exports(lay, args)
    for p, err in ex.errors.items():
        print(f"  !! {p.name} cannot be read as a {args.assay} export: {err}")
    for k in ("geom", "animal"):
        for p in ex.stale.get(k, []):
            if ex.current.get(k):
                todo = (f"--apply --rebuild keeps {ex.current[k][0].name} (it rebuilds the chain) and archives this file; "
                        f"to use this file instead, pass it with --{k}-h5 PATH")
            else:
                todo = f"--apply --rebuild rebuilds the chain from it; --{k}-h5 PATH installs it as the new export"
            print(f"  !! {k} export in place does not rebuild its locs file: {p} ({when(p)}); its own locs file differs "
                  "from what it gives (exported again under the same name?) or it is older than the chain's locs files. "
                  + todo)
    have = {k: (k in ex.new or k in ex.reuse or len(ex.inplace(k)) == 1) for k in ("geom", "animal")}
    new_kinds = [k for k in ("geom", "animal") if k in ex.new]
    work = Path(tempfile.mkdtemp(prefix="retrack_"))
    try:
        sc = scoring_setup(lay, spec, args, config, cfg_text, stage_default, output_default)
        status = CODE_STATUS(sc, spec)
        scored = scored_state(sc, lay)
        staged_md5 = md5(lay.staged) if lay.staged.exists() else None
        due, scored_text = scoring_note(lay, sc, scored, staged_md5)
        missing = [k for k in ("geom", "animal") if k not in ex.new and k not in ex.reuse]
        if new_kinds and missing:
            # one kind arrived, the other did not: never pair a new export with an old one by accident
            for k in new_kinds:
                print(f"  new {k:<6} {ex.new[k]}")
            print_waiting(lay, spec, ex, missing, args)
            return 0
        if not all(have.values()):
            print_current(lay, ex)
            print_waiting(lay, spec, ex, [k for k in ("geom", "animal") if not have[k]], args)
            print_scoring_plan(lay, spec, sc, status, args, argv, config, scored_text)
            print_followups(lay, spec, args, None, qc_file)
            return 0

        use = {k: ex.new.get(k) or (ex.inplace(k)[0] if k not in ex.reuse else None) for k in ("geom", "animal")}
        if not new_kinds:
            print(f"  no new export found: checking that the exports in place rebuild the files in the chain")
        h5_dst: dict[str, Path] = {}
        for k in ("geom", "animal"):
            if k in ex.reuse:
                print(f"  {k:<9} reuse  {ex.reuse[k]}")
                continue
            fr, vp = h5_info(use[k])
            tag = "NEW" if k in ex.new else "in place"
            print(f"  {k:<9} {tag:<8} {use[k]}   ({fr} frames; video {vp or '?'})")
            h5_dst[k] = lay.dirs[k] / lay.install_name(use[k])
            if h5_dst[k].name != use[k].name:
                print(f"            copied in as {h5_dst[k].name} (the name must contain {lay.tag} for the tool to find it again)")
            if k in ex.new or use[k] in ex.stale.get(k, []):    # a new export, or one used through --rebuild
                for why in identity_problems(lay, k, use[k], vp):
                    if args.accept_name:
                        notes.append(f"{why} (accepted: --accept-name)")
                    else:
                        problems.append(f"{why}: is it the right export? (--accept-name installs it anyway)")

        # ---- build everything in the work folder; nothing live is touched until all of it worked
        locs: dict[str, Path] = {}
        tables: dict[str, pd.DataFrame] = {}
        for k in ("geom", "animal"):
            if k in ex.reuse:
                locs[k] = ex.reuse[k]
                tables[k] = pd.read_csv(ex.reuse[k])
                continue
            df = ex.tables.get(use[k])
            if df is None:
                try:
                    df, _ = sh.extract_locs(use[k], k, lay.assay, batch=lay.batch)
                except (ValueError, OSError, KeyError) as exc:
                    raise SystemExit(f"{use[k]}: {exc}")
            tables[k] = df
            locs[k] = work / k / lay.locs_path(k, h5_dst[k].name).name
            locs[k].parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(locs[k], index=False)
        n_geom, n_animal = len(tables["geom"]), len(tables["animal"])
        if n_geom != n_animal:
            raise SystemExit(f"geom has {n_geom} frames and animal {n_animal}: export both from the same clip"
                             + (" (the reused corners belong to the old clip: track the arena on the new clip, or drop "
                                "--reuse-geom)" if ex.reuse else ""))
        g_stem = locs["geom"].name[: -len(".analysis_locs.csv")] if locs["geom"].name.endswith(".analysis_locs.csv") else locs["geom"].stem
        merged_dst = lay.merged_target(g_stem)
        merged = work / merged_dst.name
        pd.concat([pd.read_csv(locs["geom"]), pd.read_csv(locs["animal"])], axis=1).to_csv(merged, index=False)
        fmt_dst = lay.formatted_target(g_stem)
        fmt_new = work / "formatted" / fmt_dst.name
        info = sd.format_merged(merged, fmt_new, FORMAT_RENAME.get((lay.assay, lay.batch)),
                                FORMAT_DIGITS.get((lay.assay, lay.batch)))
        staged_new = work / "staged" / lay.staged.name
        staged_new.parent.mkdir(parents=True, exist_ok=True)
        if spec.neck:
            how = neck_copy(fmt_new, staged_new, args.rscript, REPO / STAGE_SCRIPT, notes, args.apply, args.python_neck_copy)
        else:
            shutil.copyfile(fmt_new, staged_new)
            how = "copied"
        nan = int(tables["geom"].isna().sum().sum() + tables["animal"].isna().sum().sum())
        print(f"  formatted {info['frames']} frames, {info['bodyparts']} landmarks, {nan} NaN after filling; staged copy: {how}")
        stage_changes = staged_md5 is None or sha(lay.staged) != sha(staged_new)

        # ---- checks on the new file: the protective ones whenever the staged file would change
        seconds = n_animal / VIDEO_FPS
        length_ok = abs(seconds - spec.expected_s) <= LENGTH_LIMIT_S
        print(f"  length    {n_animal} frames = {seconds:.1f} s at 29.97/s (protocol {spec.expected_s:g} s): "
              + ("ok" if length_ok else "OUTSIDE the light-check limit" + (" (--accept-length)" if args.accept_length else "")))
        if not length_ok and stage_changes and not args.accept_length:
            problems.append(f"length {seconds:.1f} s is not {spec.expected_s:g} +/- {LENGTH_LIMIT_S:g} s; check the cut "
                            "(add --accept-length to install it anyway)")
        parts, data = read_dlc(staged_new)
        scale = px_per_cm(spec, parts, data)
        m = medians(parts, data, ["tl", "tr", "br", "bl"])
        if spec.arena_cm is None:
            arm = float(np.hypot(*(m["tl"] - m["tr"]))) / scale
            print(f"  arena     {scale:.2f} px/cm (tl-bl = 60 cm); arm width tl-tr {arm:.2f} cm"
                  + ("" if ARM_CM[0] <= arm <= ARM_CM[1] else f"  !! outside {ARM_CM[0]}-{ARM_CM[1]} cm"))
        else:
            wpx = (np.hypot(*(m["tl"] - m["tr"])) + np.hypot(*(m["bl"] - m["br"]))) / 2
            hpx = (np.hypot(*(m["tl"] - m["bl"])) + np.hypot(*(m["tr"] - m["br"]))) / 2
            ratio = (wpx / hpx) / (spec.arena_cm[0] / spec.arena_cm[1])
            print(f"  arena     {scale:.2f} px/cm; side ratio {ratio:.3f} of the expected"
                  + ("" if abs(ratio - 1) <= ASPECT_LIMIT else f"  !! more than {ASPECT_LIMIT:.0%} off"))
        if stage_changes:
            check_copies(lay, spec, staged_new, parts, data, scale, args, problems, notes)
        side = lay.side_files()
        if side:
            raw_side, staged_side = side
            v_raw, v_staged = side_value(raw_side, lay.code), side_value(staged_side, lay.code)
            print(f"  side file {staged_side.name}: {lay.code} {v_staged or 'MISSING'} (staged), {v_raw or 'MISSING'} ({raw_side})")
            if v_staged is None or v_staged != v_raw:
                problems.append(f"{staged_side} does not list {lay.code} as the raw {raw_side.name} does: copy it as 08 "
                                f"does (Copy-Item '{raw_side}' '{staged_side}') after checking the other animals' sides")
        for n in notes:
            print(f"  note      {n}")

        # ---- what would change
        actions: list[Action] = []
        left: list[Path] = []
        for k in ("geom", "animal"):
            if k in ex.reuse:
                continue
            kdir, arch = lay.dirs[k], lay.dirs[k] / f"_superseded_{TODAY}"
            if use[k].resolve() != h5_dst[k].resolve():
                actions.append(plan_file(use[k], h5_dst[k], arch))
            locs_dst = kdir / locs[k].name
            for p in lay.mine(kdir):
                if p.resolve() in (h5_dst[k].resolve(), locs_dst.resolve()) or "_old" in p.name:
                    continue
                if p.name.endswith(".analysis.h5") or p.name.endswith("_locs.csv"):
                    actions.append(Action("move", p, archive=arch))
                else:
                    left.append(p)
            actions.append(plan_file(locs[k], locs_dst, arch))
        root_arch = lay.sleap / f"_superseded_{TODAY}"
        merged_action = plan_file(merged, merged_dst, root_arch)
        if merged_action.verb != "unchanged":               # other merged copies (B1, B3) have no reader: moved
            actions += [Action("move", p, archive=root_arch)  # only with the one they stand beside
                        for p in lay.merged_existing() if p != merged_dst]
        actions.append(merged_action)
        fmt_moves = [Action("move", p, archive=lay.fmt_archive) for p in lay.formatted_existing() if p != fmt_dst]
        fmt_action = plan_file(fmt_new, fmt_dst, lay.fmt_archive)
        actions += fmt_moves + [fmt_action]
        if fmt_moves or fmt_action.verb != "unchanged":
            actions += [Action("move", p, archive=lay.fmt_old) for p in lay.flat_copies()]
        stage_action = plan_file(staged_new, lay.staged, lay.stage_arch)
        changed = any(a.verb != "unchanged" for a in actions + [stage_action])
        print("\n".join(a.line() for a in actions + [stage_action]))
        for p in left:
            print(f"  left      {p}   (not an .analysis.h5 or _locs.csv file: stays where it is)")
        new_md5 = md5(staged_new)
        print(f"  staged md5 {staged_md5 or '-'} -> {new_md5}")
        if problems:
            for p in problems:
                print(f"  !! {p}", file=sys.stderr)

        if not new_kinds and not changed:
            names = ", ".join(f"{k} {use[k].name} ({when(use[k])})" for k in ("geom", "animal") if k not in ex.reuse)
            print(f"  in place: {names}; they rebuild every file in the chain"
                  + ("; scoring of the installed file is still due" if due else ": nothing to do"))
            print_waiting(lay, spec, ex, ["geom", "animal"] if "geom" not in ex.reuse else ["animal"], args,
                          again=bool(lay.archived_exports()))
        if not args.apply:
            print_scoring_plan(lay, spec, sc, status, args, argv, config, scored_text if not stage_changes else "")
            if not new_kinds and changed:
                print("  !! no new export, but the files above differ from what the exports in place give; "
                      "--apply --rebuild would carry it out")
            if new_kinds or changed:                            # (otherwise print_waiting listed them)
                print_strays(lay)
            print_followups(lay, spec, args, new_md5, qc_file)
            return 0

        # ---- apply
        if problems:
            raise SystemExit("not applied: " + "; ".join(problems))
        stage_changed = stage_action.verb != "unchanged"
        if not changed and not due and not args.rerun_analysis:
            print("everything is already up to date")
            return 0
        if changed and not new_kinds and not args.rebuild:
            raise SystemExit("no new export found; add --rebuild to rebuild the chain from the exports in place")
        score = not args.skip_analysis and (stage_changed or due or args.rerun_analysis)
        if score and status[0] != "same" and not args.baseline:
            print_scoring_plan(lay, spec, sc, status, args, argv, config, scored_text if not stage_changed else "")
            print("  stopped before changing anything: choose --baseline or --scoring-repo (or --skip-analysis)", file=sys.stderr)
            return 3
        previous = None
        if score and args.baseline and not stage_changed and due:
            previous = previous_staged(lay, scored[0])
            if previous is None:
                raise SystemExit(f"--baseline scores the batch on the files the current output was scored from, but the "
                                 f"staged copy with md5 {scored[0] or '?'} is not among the archived ones in "
                                 f"{lay.stage_root / lay.assay}/_superseded_*/: score with --scoring-repo instead")
        baseline = score and args.baseline and (stage_changed or due or status[0] != "same")
        check_actions(actions + [stage_action])
        do_actions(actions)
        backups: dict[Path, Path] = {}

        def install_staged() -> None:
            if stage_action.verb == "unchanged":
                return
            if not sc["manifest"].exists() and sc["out_dirs"][0].is_dir() and lay.staged.name not in read_marker(sc["marker"]):
                write_marker(sc["marker"], lay.staged.name, staged_md5 or "-",
                             "the staged file in place when the tool first replaced it (taken as the one scored)")
            do_actions([stage_action])

        def backup() -> None:
            for d in sc["out_dirs"] + ([sc["global_dir"]] if sc["global_dir"] else []):
                name = sc["tag"] + "_before" if d in sc["out_dirs"] else f"{spec.assay}_all_batches_before"
                dst = sc["backup_root"] / name
                backups[d] = dst
                if d.exists() and not dst.exists():
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copytree(d, dst)
                    print(f"  backed up {d} -> {dst}")

        cfg = work / f"{sc['tag']}.yaml"
        cfg.write_text(sc["text"], encoding="utf-8")

        def run(label: str) -> dict:
            fails = sc["global_dir"] / "failures.csv" if sc["global_dir"] is not None else None
            if fails is not None and fails.exists():          # read only a failures.csv this run writes
                aside = free_path(sc["backup_root"] / f"failures_before_{lay.batch}_{label}_run.csv")
                aside.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(fails), str(aside))
                print(f"  moved     {fails} -> {aside} (an earlier run's)")
            print(f"  running   {sc['script'].name} for {lay.batch}" + (f" {lay.phase}" if lay.phase else "") + f" ({label}) ...")
            res = SCORER(args.rscript, sc["script"], cfg, sc["repo"])
            try:
                if res.returncode != 0:
                    print(res.stdout[-1500:], res.stderr[-1500:], file=sys.stderr)
                    raise SystemExit(f"{sc['script'].name} failed (exit {res.returncode}); the previous output is in "
                                     f"{sc['backup_root']}. The files installed so far stay; once the cause is fixed, run "
                                     "the same command again: the plan shows what is left (with --rebuild if files of "
                                     "the chain still differ, or the scoring that is still due)")
                if fails is not None and fails.exists():
                    text = fails.read_text(encoding="utf-8", errors="replace")
                    if lay.code in text:
                        raise SystemExit(f"{lay.code} is in {fails}, written by this run: the scoring did not take the "
                                         "file (the scoring of the installed file is still due)")
                    print(f"  note      {fails} lists other files of this run")
                return snapshot(sc["out_dirs"][0], spec)
            finally:
                restore_global(label)

        def restore_global(label: str) -> None:
            g = sc["global_dir"]
            if g is None or g not in backups:
                return
            written = free_path(sc["backup_root"] / f"{spec.assay}_all_batches_written_by_{lay.batch}_{label}_run")
            if g.exists():
                shutil.move(str(g), str(written))
            if backups[g].exists():
                shutil.copytree(backups[g], g)
                print(f"  put back  {g} (the {lay.batch}-only version is in {written})")
            elif written.exists():                              # there was none before the run: none after it either
                print(f"  moved     the {lay.batch}-only {g.name} written by the run to {written} (there was none before)")

        if not score:
            install_staged()
            print("  scoring not run" + (" (--skip-analysis)" if args.skip_analysis else "")
                  + ("; the scoring of the installed file is still due" if stage_changed or due else ""))
            print_followups(lay, spec, args, new_md5, qc_file)
            return 0
        backup()
        before = snapshot(backups[sc["out_dirs"][0]], spec) if backups[sc["out_dirs"][0]].exists() else {}
        if args.baseline and not baseline:
            print(f"  note      nothing new to score and the scoring code is the same: no baseline run, compared "
                  f"against {backups[sc['out_dirs'][0]]}")
        if baseline and previous is not None:
            # the installed file is not scored yet: score the batch once with the file it replaced, then put it back
            keep = free_path(sc["backup_root"] / f"{sc['tag']}_staged_during_baseline" / lay.staged.name)
            keep.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(lay.staged, keep)
            print(f"  baseline  with {previous} in place of {lay.staged} (kept in {keep} until it is back)")
            put_file(previous, lay.staged)
            try:
                before = run("baseline")
            finally:
                put_file(keep, lay.staged)
        elif baseline:
            before = run("baseline")
        if baseline:
            snap = free_path(sc["backup_root"] / f"{sc['tag']}_baseline")
            if sc["out_dirs"][0].exists():
                shutil.copytree(sc["out_dirs"][0], snap)
        install_staged()
        after = run("after")
        if not sc["manifest"].exists() and sc["out_dirs"][0].is_dir():
            write_marker(sc["marker"], lay.staged.name, md5(lay.staged), "scored by the tool's 'after' run")
        per_table, others = compare(before, after, spec, lay.code)
        for name, (n, codes, _) in per_table.items():
            print(f"  {name:<38} {n} rows; animals whose rows changed: {codes or 'none'}")

        def row_of(snap: dict) -> dict:                     # the animal's row, first two tables together
            out: dict = {}
            for t in reversed(spec.tables[:2]):
                out.update(next((r for k, r in snap.get(t.name, {}).items() if k[0][:4] == lay.code), {}))
            return out
        row, old = row_of(after), row_of(before)
        shown = [f"{c} {old.get(c, '-')} -> {row[c]}" for c in spec.shows if c in row]
        print(f"  {lay.code}: " + "; ".join(shown))
        if others:
            print(f"  !! animals other than {lay.code} changed: {sorted(others)} -- check before using this output", file=sys.stderr)
            return 2
        print_followups(lay, spec, args, new_md5, qc_file)
        return 0
    finally:
        shutil.rmtree(work, ignore_errors=True)


def check_copies(lay: Layout, spec: Spec, staged_new: Path, parts: list[str], data: np.ndarray, scale: float,
                 args: argparse.Namespace, problems: list[str], notes: list[str]) -> None:
    """The same trial twice, as 31_light_checks.R finds it: a byte copy of any staged file of the assay, or the body
    centre within DUPLICATE_CM of any staged file of the same assay and batch (every phase, every animal) at an offset
    of up to DUPLICATE_LAG frames."""
    assay_dir = lay.stage_root / lay.assay
    globs = [f"B*/{ph}/*.csv" for ph in spec.phases] or ["B*/*.csv"]     # not the open-field output folder in B*/
    staged_all = [p for g in globs for p in sorted(assay_dir.glob(g)) if p != lay.staged]
    size, digest = staged_new.stat().st_size, sha(staged_new)
    copies = [p for p in staged_all if p.stat().st_size == size and sha(p) == digest]
    if copies and not args.accept_duplicate:
        problems.append(f"the new file is byte-identical to {[str(p) for p in copies]} (add --accept-duplicate if intended)")
    elif copies:
        notes.append(f"byte-identical to {[str(p) for p in copies]} (accepted)")
    batch_dir = assay_dir / lay.batch
    others = [p for p in staged_all if batch_dir in p.parents]
    body = point(parts, data, "bodycentre") / scale
    closest = None
    for other in others:
        same_animal = lay.phase and other.name == f"{lay.code}_{other.parent.name}.csv"
        if same_animal:
            oparts, odata = read_dlc(other)
            ob = point(oparts, odata, "bodycentre") / px_per_cm(spec, oparts, odata)
        else:
            try:
                ob = body_cm(other, spec)
            except (KeyError, ValueError) as exc:
                notes.append(f"{other.name}: not compared ({exc})")
                continue
        dist, lag = closest_offset(body, ob)
        what = f"{other.parent.name}/{other.stem}" if lay.phase else other.stem
        if same_animal:
            om = medians(oparts, odata, ["tl", "tr", "br", "bl", "socl", "socr"])
            shift = max(float(np.hypot(*(om[k] - mv))) for k, mv in medians(parts, data, list(om)).items()) / scale
            print(f"  vs {other.parent.name:<6} body centre at best {dist:.2f} cm apart (offset {lag} frames); "
                  f"arena points moved up to {shift:.2f} cm")
        if closest is None or dist < closest[0]:
            closest = (dist, lag, what)
        if dist < DUPLICATE_CM:
            msg = (f"the new {lay.phase} track coincides with {other.parent.name} at an offset of {lag} frames"
                   if same_animal else f"the new track coincides with {what} at an offset of {lag} frames")
            msg += f" (median {dist:.3f} cm): the same trial twice"
            if args.accept_duplicate:
                notes.append(msg + " (accepted)")
            else:
                problems.append(msg + " (add --accept-duplicate if intended)")
    if closest is not None:
        print(f"  copies    body centre against {len(others)} staged {lay.assay} {lay.batch} files"
              + (" (every phase)" if lay.phase else "") + f": closest {closest[2]} at {closest[0]:.2f} cm "
              f"(offset {closest[1]} frames; the same trial twice below {DUPLICATE_CM} cm)")


# --------------------------------------------------------------------------- #
# command line
# --------------------------------------------------------------------------- #

def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--batch", required=True)
    ap.add_argument("--assay", required=True, choices=["NOR", "EPM", "OFT", "SocP"])
    ap.add_argument("--code", required=True, help="animal code, e.g. H1E2")
    ap.add_argument("--phase", choices=sorted({p for v in sd.PHASES.values() for p in v}),
                    help="NOR: NOV (novel object) or HAB; SocP: HAB, S1 or S2; EPM and OFT: none")
    ap.add_argument("--static-landmark", nargs="*", default=[], metavar="NODE", help="NOR: geom landmarks held at their first labelled position, e.g. objR")
    ap.add_argument("--root", "--raw-root", dest="raw_root", type=Path, default=sd.BEHAVIOR_ROOT, help="Raw Data/Behavior root")
    ap.add_argument("--stage-root", type=Path, help="staging root instead of the config's (exp9_sleap_staging)")
    ap.add_argument("--output-root", type=Path, help="scored-output root instead of the config's (sleap_output_all)")
    ap.add_argument("--qc-root", type=Path, default=QC_ROOT, help="folder of retrack_list.csv and leave_out_list.csv")
    ap.add_argument("--decisions", type=Path, default=DECISIONS, help=argparse.SUPPRESS)
    ap.add_argument("--config", type=Path, help="scoring config (default: the assay's *_all.yaml in validation/exp9_boris/config)")
    ap.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    ap.add_argument("--geom-h5", type=Path, help="the new arena export, if not (or not alone) in SLEAP/geom/")
    ap.add_argument("--animal-h5", type=Path, help="the new animal export, if not (or not alone) in SLEAP/animal/")
    ap.add_argument("--reuse-geom", action="store_true", help="EPM/OFT/SocP: keep the geom locs in place (arena not re-tracked)")
    ap.add_argument("--clip", help="OFT B1: the GoPro clip the animal's files are named by, e.g. GX010804")
    ap.add_argument("--accept-length", action="store_true", help="install a recording outside the protocol length")
    ap.add_argument("--accept-duplicate", action="store_true", help="install a recording that looks like a copy")
    ap.add_argument("--accept-name", action="store_true",
                    help="install a new export whose file name or video_path does not name the animal (OFT B1: the clip)")
    ap.add_argument("--accept-clip", action="store_true",
                    help="OFT B1: go on although the merged file of --clip does not give the animal's formatted file")
    ap.add_argument("--python-neck-copy", action="store_true",
                    help="EPM with --apply: write the staged copy in Python when Rscript or 08 is not found")
    ap.add_argument("--baseline", action="store_true", help="score the batch once before installing, compare against that")
    ap.add_argument("--scoring-repo", type=Path, help="score from this checkout (e.g. a worktree at the commit of the current output)")
    ap.add_argument("--rebuild", action="store_true", help="with --apply and no new export: rebuild the chain from the exports in place")
    ap.add_argument("--skip-analysis", action="store_true", help="stop after staging")
    ap.add_argument("--rerun-analysis", action="store_true", help="run the analysis even if no file changed")
    ap.add_argument("--apply", action="store_true", help="make the changes (default: print the plan)")
    args = ap.parse_args(argv)
    if args.assay == "NOR":
        if args.phase not in sd.PHASES["NOR"]:
            ap.error("NOR needs --phase NOV or HAB")
        extra = [o for o, v in (("--reuse-geom", args.reuse_geom), ("--clip", args.clip), ("--baseline", args.baseline),
                                ("--scoring-repo", args.scoring_repo), ("--rebuild", args.rebuild),
                                ("--accept-length", args.accept_length), ("--accept-duplicate", args.accept_duplicate),
                                ("--accept-name", args.accept_name), ("--accept-clip", args.accept_clip),
                                ("--python-neck-copy", args.python_neck_copy)) if v]
        if extra:
            ap.error(f"{', '.join(extra)}: options for EPM, OFT and SocP")
        return run_nor(args)
    return run_assay(args, argv)


if __name__ == "__main__":
    raise SystemExit(main())
