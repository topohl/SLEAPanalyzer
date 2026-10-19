"""Re-run ONE re-tracked recording through the whole SLEAPanalyzer chain, in one command.

    h5 export -> locs.csv (geom + animal) -> formatted CSV -> staged copy -> assay analysis -> before/after report

    python tools/retrack_recording.py --batch B5 --assay NOR --code H1E2 --phase NOV --static-landmark objR            # plan only
    python tools/retrack_recording.py --batch B5 --assay NOR --code H1E2 --phase NOV --static-landmark objR --apply    # do it

What it does with --apply, and how nothing is lost:
  * reads the ONE ``*_<CODE>_<PHASE>.analysis.h5`` in geom/ and animal/ (names containing ``_old`` are ignored; two matches stop it),
  * writes ``..._<CODE>_<PHASE>.analysis_locs.csv`` next to them (``sleap_h5.py``; ``--static-landmark`` holds a hand-labelled static
    object at its first labelled position instead of interpolating it),
  * formats and installs ``formatted/<CODE>_<PHASE>.merged_locs_formatted.csv`` (``sleap_dataform.process_one``),
  * moves every file it replaces into ``_superseded_<date>/`` (locs csv, staged copy) or ``formatted/old/`` (formatted files; the audit script
    07 ignores only that folder name), including a legacy-labelled file of the same animal such as ``H1E2_NOR...`` that staging would otherwise
    pick up first,
  * copies the formatted file to ``<config input_dir>/<batch>/<CODE>_<PHASE>.csv``,
  * backs up the batch's analysis output, runs ``02_SLEAPanalyzer/DLCA_NOR v1.2.1.R`` for that batch only, and reports which animals' rows
    changed (all others must be identical).
NOR only: the staging names of the other assays differ. Anything already up to date is left alone, so a second run changes nothing.
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
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "01_SLEAPcoords"))
import sleap_dataform as sd  # noqa: E402
import sleap_h5 as sh  # noqa: E402

DEFAULT_CONFIG = REPO / "validation" / "exp9_boris" / "config" / "nor_all.yaml"
NOR_SCRIPT = REPO / "02_SLEAPanalzyer" / "DLCA_NOR v1.2.1.R"
DEFAULT_RSCRIPT = Path(r"C:\Users\topohl\AppData\Local\Programs\R\R-4.5.1\bin\x64\Rscript.exe")
TODAY = dt.date.today().isoformat()


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def config_value(config: Path, key: str) -> str:
    m = re.search(rf"^{key}:\s*(.+?)\s*$", config.read_text(encoding="utf-8"), re.M)
    if not m:
        raise SystemExit(f"{config.name} has no '{key}:' line")
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


def changed_cells(before: dict[str, dict[str, str]], after: dict[str, dict[str, str]]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for code in sorted(set(before) | set(after)):
        b, a = before.get(code), after.get(code)
        if b is None or a is None:
            out[code] = ["<row missing>"]
            continue
        cols = []
        for c in a:
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


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--batch", required=True)
    ap.add_argument("--assay", required=True, choices=["NOR"])
    ap.add_argument("--code", required=True, help="animal code, e.g. H1E2")
    ap.add_argument("--phase", required=True, choices=sd.PHASES["NOR"], help="NOV for the novel-object trial")
    ap.add_argument("--static-landmark", nargs="*", default=[], metavar="NODE", help="geom landmarks held at their first labelled position, e.g. objR")
    ap.add_argument("--root", type=Path, default=sd.BEHAVIOR_ROOT, help="Raw Data/Behavior root")
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="NOR config; its input_dir/output_dir define staging and output")
    ap.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    ap.add_argument("--geom-h5", type=Path)
    ap.add_argument("--animal-h5", type=Path)
    ap.add_argument("--skip-analysis", action="store_true", help="stop after staging")
    ap.add_argument("--rerun-analysis", action="store_true", help="run the analysis even if no file changed")
    ap.add_argument("--apply", action="store_true", help="make the changes (default: print the plan)")
    args = ap.parse_args(argv)

    sleap = args.root / args.batch / args.assay / "SLEAP"
    geom_dir, animal_dir, fmt_dir = sleap / "geom", sleap / "animal", sleap / "formatted"
    stage_dir = Path(config_value(args.config, "input_dir")) / args.batch
    out_dir = Path(config_value(args.config, "output_dir")) / args.batch
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
        cfg.write_text(re.sub(r"^batches:.*$", f"batches: [{args.batch}]", args.config.read_text(encoding="utf-8"), flags=re.M), encoding="utf-8")
        env = dict(os.environ, SLEAP_ANALYZER_CONFIG=str(cfg))
        print(f"  running   {NOR_SCRIPT.name} for {args.batch} ...")
        res = subprocess.run([str(args.rscript), str(NOR_SCRIPT)], cwd=REPO, env=env, capture_output=True, text=True)
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


if __name__ == "__main__":
    raise SystemExit(main())
