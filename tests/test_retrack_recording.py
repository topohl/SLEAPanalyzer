"""Tests for tools/retrack_recording.py on synthetic exports in temporary roots.

The scoring is replaced by a stub that writes small tables (and, as the R scripts do outside the open field, a
run_manifest.yaml) from the staged files, and the EPM neck copy by its Python version: nothing here runs R. The stub
refuses to write outside the test's temporary folder. Run from this folder: python -m unittest test_retrack_recording
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "01_SLEAPcoords"))
sys.path.insert(0, str(REPO / "tools"))

import retrack_recording as rr  # noqa: E402
import sleap_dataform as sd  # noqa: E402
import sleap_h5 as sh  # noqa: E402

TODAY = rr.TODAY
N = 600                                                     # frames per synthetic recording

ARENA = {
    "EPM": {"tl": (400, 100), "tr": (475, 100), "ctr": (475, 520), "rt": (900, 520), "rb": (900, 595),
            "cbr": (475, 595), "br": (475, 1000), "bl": (400, 1000), "cbl": (400, 595), "lb": (0, 595),
            "lt": (0, 520), "ctl": (400, 520)},
    "OFT": {"tl": (100, 100), "tr": (900, 100), "bl": (100, 900), "br": (900, 900)},
    "SocP": {"tl": (100, 100), "tr": (900, 100), "br": (900, 536), "bl": (100, 536), "socl": (200, 318), "socr": (800, 318)},
    "NOR": {"tl": (100, 100), "tr": (900, 100), "br": (900, 900), "bl": (100, 900), "objL": (300, 500), "objR": (700, 500)},
}
CENTRE = {"EPM": (437, 557), "OFT": (500, 500), "SocP": (500, 318), "NOR": (500, 500)}


def stored(assay: str, kind: str) -> list[str]:
    return sh.H5_NODES.get((assay, kind)) or sh.NODES[(assay, kind)]


def video_for(name: str) -> str:
    """The video path SLEAP stores for an export of this name: the video name follows the labels file and <NNN>."""
    stem = name[: -len(".analysis.h5")] if name.endswith(".analysis.h5") else name
    return "//share/cropped/" + re.sub(r"^.*\.\d+_", "", stem) + ".mp4"


def write_h5(path: Path, xy: np.ndarray, names: list[str], video: str | None = None, tracks: int = 1) -> Path:
    """A SLEAP analysis file from (frame, node, xy) coordinates."""
    t = np.repeat(xy.transpose(2, 1, 0)[None], tracks, axis=0)          # (track, xy, node, frame)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f["tracks"] = t
        f["node_names"] = np.array([n.encode() for n in names])
        f["video_path"] = video or video_for(path.name)
    return path


def geom_xy(assay: str, n: int = N, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    pts = np.array(list(ARENA[assay].values()), dtype=float)
    xy = pts[None] + rng.normal(0, 0.3, (n, len(pts), 2))
    xy[3:6, 0] = np.nan                                     # a short gap in one corner
    return xy


def animal_xy(assay: str, seed: int, n: int = N) -> np.ndarray:
    rng = np.random.default_rng(seed)
    k = len(sh.NODES[(assay, "animal")])
    path = np.cumsum(rng.normal(0, 3, (n, 2)), axis=0) + np.array(CENTRE[assay], dtype=float)
    xy = path[:, None, :] + rng.normal(0, 12, (1, k, 2)) + rng.normal(0, 0.5, (n, k, 2))
    xy[:, 9] = path                                         # node 9 is the body centre
    xy[10:14, 0] = np.nan                                   # an interior gap
    xy[0, 3] = np.nan                                       # an empty first frame
    return xy


def shifted(xy: np.ndarray, frames: int) -> np.ndarray:
    """The same track *frames* frames on (the end repeated to keep the length)."""
    return np.concatenate([xy[frames:], xy[-frames:]])


def tree_state(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


class Chain:
    """Builds the current chain of one recording the way the corpus was built (locs, merged, formatted, staged)."""

    def __init__(self, tmp: Path) -> None:
        self.raw, self.stage, self.out, self.qc = tmp / "raw", tmp / "stage", tmp / "out", tmp / "qc"
        self.qc.mkdir(parents=True)
        self.codes: list[str] = []

    def export(self, assay, batch, kind, name, xy, video=None, where=None) -> Path:
        d = where or (self.raw / batch / assay / "SLEAP" / kind)
        return write_h5(d / name, xy, stored(assay, kind), video)

    def locs(self, assay, batch, kind, h5: Path, locs_name: str | None = None) -> Path:
        df, _ = sh.extract_locs(h5, kind, assay, batch=batch)
        p = self.raw / batch / assay / "SLEAP" / kind / (locs_name or h5.name[:-3] + "_locs.csv")
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, index=False)
        return p

    def finish(self, assay, batch, code, phase, g_locs: Path, a_locs: Path, merged_name, fmt_name) -> Path:
        sleap = self.raw / batch / assay / "SLEAP"
        merged = sleap / merged_name
        pd.concat([pd.read_csv(g_locs), pd.read_csv(a_locs)], axis=1).to_csv(merged, index=False)
        fmt_dir = sleap / "formatted" / phase if phase else sleap / "formatted"
        fmt = fmt_dir / fmt_name
        sd.format_merged(merged, fmt, rr.FORMAT_RENAME.get((assay, batch)), rr.FORMAT_DIGITS.get((assay, batch)))
        tail = [batch] + ([phase] if phase else [])
        staged = self.stage.joinpath(assay, *tail, f"{code}_{phase}.csv" if phase else f"{code}.csv")
        staged.parent.mkdir(parents=True, exist_ok=True)
        if rr.SPECS[assay].neck:
            rr.neck_copy_python(fmt, staged)
        else:
            staged.write_bytes(fmt.read_bytes())
        if code not in self.codes:
            self.codes.append(code)
        (self.stage / "animalIDCode.txt").write_text("ID\tCode\n" + "".join(f"X{i}\t{c}\n" for i, c in enumerate(self.codes)))
        return fmt

    def args(self, *extra: str) -> list[str]:
        return ["--raw-root", str(self.raw), "--stage-root", str(self.stage), "--qc-root", str(self.qc),
                "--rscript", str(self.raw / "no_R.exe"), *extra]


class FakeScorer:
    """Writes the tables the before/after report reads, as a function of the staged files only.

    Outside the open field it also writes run_manifest.yaml with the md5 of every input, as the R scripts do (without
    a commit, so the code version stays unknown). It writes only inside *root*."""

    def __init__(self, root: Path | None = None) -> None:
        self.root = root.resolve() if root else None
        self.calls: list[str] = []
        self.bump_on_call: int | None = None
        self.bump_code: str | None = None
        self.fail_on_call: int | None = None
        self.failures_on_call: int | None = None             # open field: list this code in failures.csv
        self.failures_code: str | None = None

    def guard(self, p: Path) -> Path:
        p = Path(p).resolve()
        if self.root is not None and self.root not in p.parents:
            raise AssertionError(f"the scoring stub was pointed outside the test folder: {p}")
        return p

    def __call__(self, rscript, script, config, cwd):
        text = Path(config).read_text(encoding="utf-8")
        self.calls.append(text)
        if len(self.calls) == self.fail_on_call:
            return subprocess.CompletedProcess([], 1, "", "R error")
        val = lambda k: re.search(rf"^{k}:\s*(.+?)\s*$", text, re.M).group(1).strip("'\"")  # noqa: E731
        batches = [b.strip() for b in val("batches").strip("[]").split(",")]
        bump = (len(self.calls) == self.bump_on_call)
        for b in batches:
            if "DLCA_OFT" in script.name:
                self.oft(Path(val("behavior_root")), b, bump)
            elif "DLCA_SocP" in script.name:
                for ph in [p.strip() for p in val("phases").strip("[]").split(",")]:
                    self.table(Path(val("input_dir")) / b / ph, Path(val("output_dir")) / b / ph, "combined_output.csv", bump)
            elif "DLCA_EPM" in script.name:
                self.table(Path(val("input_dir")) / b, Path(val("output_dir")) / b, "Report.csv", bump, key="file")
            else:
                self.table(Path(val("input_dir")) / b, Path(val("output_dir")) / b, "combined_output.csv", bump)
        return subprocess.CompletedProcess([], 0, "", "")

    def value(self, f: Path, bump: bool) -> float:
        parts, data = rr.read_dlc(f)
        v = float(np.nanmean(rr.point(parts, data, "bodycentre")[:, 0]))
        return v + (1.0 if bump and f.name[:4] == self.bump_code else 0.0)

    def table(self, indir: Path, outdir: Path, name: str, bump: bool, key: str = "Code") -> None:
        outdir = self.guard(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        files = sorted(indir.glob("*.csv"))
        rows = [{key: (f.name if key == "file" else f.name[:4]), "frames": len(rr.read_dlc(f)[1]), "value": self.value(f, bump)}
                for f in files]
        pd.DataFrame(rows).to_csv(outdir / name, index=False)
        if name != "combined_output.csv" or "SocP" in str(outdir):
            pd.DataFrame([{"file": r.get("file", r.get("Code", "") + ".csv"), "frames": r["frames"]} for r in rows]).to_csv(
                outdir / "tracking_qc.csv", index=False)
        lines = ["code:", "  dirty_working_tree: no", "inputs:"]
        for f in files:
            lines += [f"  - path: {f.as_posix()}", "    exists: yes", f"    md5: {md5(f)}"]
        (outdir / "run_manifest.yaml").write_text("\n".join(lines) + "\n", encoding="utf-8")

    def oft(self, root: Path, b: str, bump: bool) -> None:
        out = self.guard(root / b / "output_v1.2.0")
        out.mkdir(parents=True, exist_ok=True)
        files = sorted((root / b).glob("*.csv"))
        vals = [self.value(f, bump) for f in files]
        z = (np.array(vals) - np.mean(vals)) / (np.std(vals) or 1)
        rows = [{"file": f.stem, "Code": f.stem[:4], "distance_cm": v, "z_center_time_percent": zz}
                for f, v, zz in zip(files, vals, z)]
        pd.DataFrame(rows).to_csv(out / "combined_summary.csv", index=False)
        g = self.guard(root / "OFT_SLEAP_output_v1.2.0")
        g.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(g / "all_batches_summary.csv", index=False)
        (g / f"only_{b}.txt").write_text("written by a one-batch run")
        if len(self.calls) == self.failures_on_call:
            pd.DataFrame([{"file": f"{self.failures_code}.csv", "error": "no frames"}]).to_csv(g / "failures.csv", index=False)


def run_tool(argv: list[str]) -> tuple[int, str]:
    buf, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
        try:
            code = rr.main(argv)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
            err.write(str(exc.code))
    return code, buf.getvalue() + err.getvalue()


class Stubs(unittest.TestCase):
    def setUp(self) -> None:
        self.saved = (rr.SCORER, rr.NECK_COPIER, rr.CODE_STATUS)
        self.tmpdir = tempfile.TemporaryDirectory()
        self.tmp = Path(self.tmpdir.name)
        rr.NECK_COPIER = rr.neck_copy_python
        self.scorer = FakeScorer(self.tmp)
        rr.SCORER = self.scorer

    def tearDown(self) -> None:
        rr.SCORER, rr.NECK_COPIER, rr.CODE_STATUS = self.saved
        self.tmpdir.cleanup()


class FakeScorerTest(Stubs):
    def test_the_stub_writes_only_inside_the_test_folder(self) -> None:
        cfg = self.tmp / "cfg.yaml"
        cfg.write_text(f"input_dir: {(self.tmp / 'in').as_posix()}\noutput_dir: {Path(tempfile.gettempdir()).as_posix()}/elsewhere\n"
                       "batches: [B6]\n")
        with self.assertRaises(AssertionError):
            self.scorer(None, Path("DLCA_EPM v1.0.0.R"), cfg, None)


# --------------------------------------------------------------------------- #
# EPM
# --------------------------------------------------------------------------- #

class EpmTest(Stubs):
    def build(self) -> Chain:
        c = Chain(self.tmp)
        for i, code in enumerate(["A2C8", "R8I2"]):
            stem = f"{i:03d}_E9_B6_{code}_EPM" if code == "A2C8" else f"013_E9_B6_{code}_EPM"
            g = c.export("EPM", "B6", "geom", f"E9_B6_EPM-geom.{stem}.analysis.h5", geom_xy("EPM", seed=i))
            a = c.export("EPM", "B6", "animal", f"E9_B6_EPM-animal.{stem}.analysis.h5", animal_xy("EPM", seed=10 + i))
            c.finish("EPM", "B6", code, "", c.locs("EPM", "B6", "geom", g), c.locs("EPM", "B6", "animal", a),
                     f"E9_B6_EPM-geom.{stem}.merged_locs.csv", f"{code}_EPM.merged_locs_formatted.csv")
        self.scorer(None, Path("DLCA_EPM v1.0.0.R"), self.cfg(c), None)     # the scored output that exists today
        self.scorer.calls.clear()
        return c

    def cfg(self, c: Chain) -> Path:
        p = self.tmp / "epm_cfg.yaml"
        p.write_text(f"input_dir: {(c.stage / 'EPM').as_posix()}\noutput_dir: {(c.out / 'EPM').as_posix()}\nbatches: [B6]\n")
        return p

    def add_new_export(self, c: Chain, n: int = N, seed: int = 99, video: str | None = None) -> None:
        c.export("EPM", "B6", "geom", "E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.analysis.h5", geom_xy("EPM", n, seed=5), video=video)
        c.export("EPM", "B6", "animal", "E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5", animal_xy("EPM", seed, n), video=video)

    def base(self, c: Chain, *extra: str) -> list[str]:
        return c.args("--batch", "B6", "--assay", "EPM", "--code", "R8I2", "--output-root", str(c.out), *extra)

    def test_with_only_the_export_in_place_the_plan_says_where_the_new_one_goes(self) -> None:
        c = self.build()
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c))
        self.assertEqual(code, 0, out)
        self.assertIn("rebuild every file in the chain", out)
        self.assertIn("WAITING FOR THE NEW EXPORT", out)
        self.assertIn(str(c.raw / "B6" / "EPM" / "SLEAP" / "geom"), out)
        self.assertIn("E9_B6_EPM-geom.<NNN>_E9_B6_R8I2_EPM.analysis.h5", out)
        self.assertIn("do not export under the name of the geom export in place (E9_B6_EPM-geom.013_E9_B6_R8I2_EPM.analysis.h5)", out)
        self.assertIn("the current output was scored from the staged file in place", out)
        self.assertIn(str(rr.CONFIG_DIR / "epm_all.yaml"), out)                  # the full config path
        self.assertEqual(tree_state(self.tmp), before)
        self.assertEqual(self.scorer.calls, [])

    def test_new_export_is_installed_old_files_are_kept_and_only_this_animal_changes(self) -> None:
        c = self.build()
        self.add_new_export(c)
        base = self.base(c, "--accept-length")
        before = tree_state(self.tmp)
        code, out = run_tool(base)
        self.assertEqual(code, 0, out)
        self.assertEqual(tree_state(self.tmp), before)          # a plan writes nothing
        self.assertIn("--baseline", out)
        self.assertIn("copies    body centre against 1 staged EPM B6 files: closest A2C8", out)
        self.assertIn("'^EPM [A-Za-z0-9]+ B6 '", out)            # the light check takes the whole assay and batch

        code, out = run_tool(base + ["--apply", "--baseline"])
        self.assertEqual(code, 0, out)
        sleap = c.raw / "B6" / "EPM" / "SLEAP"
        arch = f"_superseded_{TODAY}"
        for kind in ("geom", "animal"):
            self.assertTrue((sleap / kind / arch / f"E9_B6_EPM-{kind}.013_E9_B6_R8I2_EPM.analysis.h5").exists())
            self.assertTrue((sleap / kind / arch / f"E9_B6_EPM-{kind}.013_E9_B6_R8I2_EPM.analysis_locs.csv").exists())
            self.assertTrue((sleap / kind / f"E9_B6_EPM-{kind}.019_E9_B6_R8I2_EPM.analysis_locs.csv").exists())
        self.assertTrue((sleap / arch / "E9_B6_EPM-geom.013_E9_B6_R8I2_EPM.merged_locs.csv").exists())
        self.assertTrue((sleap / "E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.merged_locs.csv").exists())
        self.assertTrue((sleap / "formatted" / "old" / "R8I2_EPM.merged_locs_formatted.csv").exists())
        self.assertEqual(sorted(p.name for p in (sleap / "formatted").glob("*.csv")),
                         ["A2C8_EPM.merged_locs_formatted.csv", "R8I2_EPM.merged_locs_formatted.csv"])
        self.assertTrue((c.stage / "EPM" / arch / "B6" / "R8I2.csv").exists())
        parts, data = rr.read_dlc(c.stage / "EPM" / "B6" / "R8I2.csv")
        self.assertTrue(np.array_equal(rr.point(parts, data, "neck"), rr.point(parts, data, "spine1")))
        self.assertEqual(len(self.scorer.calls), 2)              # baseline, then after
        self.assertIn("batches: [B6]", self.scorer.calls[0])
        self.assertIn("animals whose rows changed: ['R8I2']", out)
        self.assertTrue((self.tmp / f"_superseded_{TODAY}_R8I2" / "EPM_B6_before" / "Report.csv").exists())
        self.assertTrue((self.tmp / f"_superseded_{TODAY}_R8I2" / "EPM_B6_baseline" / "Report.csv").exists())

        # a second run finds the new export in place and changes nothing
        state = tree_state(self.tmp)
        code, out = run_tool(base + ["--apply", "--baseline"])
        self.assertEqual(code, 0, out)
        self.assertIn("up to date", out)
        self.assertEqual(tree_state(self.tmp), state)
        self.assertEqual(len(self.scorer.calls), 2)
        # and the plan names the export in place instead of waiting for one
        code, out = run_tool(base)
        self.assertNotIn("WAITING FOR THE NEW EXPORT", out)
        self.assertIn("in place: geom E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.analysis.h5 (", out)
        self.assertIn("if you re-track again:", out)

    def test_without_a_baseline_an_unknown_code_version_stops_before_any_change(self) -> None:
        c = self.build()
        self.add_new_export(c)
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--accept-length", "--apply"))
        self.assertEqual(code, 3, out)
        self.assertIn("stopped before changing anything", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_unchanged_scoring_code_compares_against_the_backup(self) -> None:
        rr.CODE_STATUS = lambda sc, spec: ("same", "same code")
        c = self.build()
        self.add_new_export(c)
        code, out = run_tool(self.base(c, "--accept-length", "--apply"))
        self.assertEqual(code, 0, out)
        self.assertEqual(len(self.scorer.calls), 1)
        self.assertIn("animals whose rows changed: ['R8I2']", out)

    def test_another_animal_changing_is_reported(self) -> None:
        self.scorer.bump_on_call, self.scorer.bump_code = 2, "A2C8"
        c = self.build()
        self.add_new_export(c)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--baseline"))
        self.assertEqual(code, 2, out)
        self.assertIn("animals other than R8I2 changed: ['A2C8']", out)

    def test_length_outside_the_protocol_is_not_installed(self) -> None:
        c = self.build()
        self.add_new_export(c)
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertNotEqual(code, 0)
        self.assertIn("is not 300 +/- 2 s", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_new_exports_of_different_lengths_are_refused(self) -> None:
        c = self.build()
        c.export("EPM", "B6", "geom", "E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.analysis.h5", geom_xy("EPM", N, seed=5))
        c.export("EPM", "B6", "animal", "E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5", animal_xy("EPM", 99, N + 30))
        code, out = run_tool(self.base(c))
        self.assertNotEqual(code, 0)
        self.assertIn("export both from the same clip", out)

    def test_only_the_animal_export_arrived_so_the_tool_waits(self) -> None:
        c = self.build()
        c.export("EPM", "B6", "animal", "E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5", animal_xy("EPM", 99))
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--apply"))
        self.assertEqual(code, 0, out)
        self.assertIn("no new geom export", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_two_new_exports_need_a_choice(self) -> None:
        c = self.build()
        self.add_new_export(c)
        c.export("EPM", "B6", "geom", "E9_B6_EPM-geom.020_E9_B6_R8I2_EPM.analysis.h5", geom_xy("EPM", seed=6))
        code, out = run_tool(self.base(c))
        self.assertNotEqual(code, 0)
        self.assertIn("--geom-h5", out)

    def test_an_export_with_another_node_order_is_refused(self) -> None:
        c = self.build()
        self.add_new_export(c)
        names = list(stored("EPM", "animal"))
        names[1], names[2] = names[2], names[1]
        write_h5(c.raw / "B6" / "EPM" / "SLEAP" / "animal" / "E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5",
                 animal_xy("EPM", 99), names)
        code, out = run_tool(self.base(c, "--accept-length"))
        self.assertIn("cannot be read as a EPM export", out)
        self.assertNotEqual(code, 0)

    def test_files_that_differ_from_the_exports_in_place_need_rebuild(self) -> None:
        c = self.build()
        fmt = c.raw / "B6" / "EPM" / "SLEAP" / "formatted" / "R8I2_EPM.merged_locs_formatted.csv"
        fmt.write_bytes(fmt.read_bytes().replace(b"\r\n0,", b"\r\n0,0", 1))        # an edited formatted file
        args = self.base(c)
        code, out = run_tool(args)
        self.assertEqual(code, 0, out)
        self.assertIn("--apply --rebuild would carry it out", out)
        before = tree_state(self.tmp)
        code, out = run_tool(args + ["--apply", "--skip-analysis"])
        self.assertNotEqual(code, 0)
        self.assertEqual(tree_state(self.tmp), before)
        code, out = run_tool(args + ["--apply", "--skip-analysis", "--rebuild"])
        self.assertEqual(code, 0, out)
        self.assertTrue((fmt.parent / "old" / fmt.name).exists())
        self.assertEqual(self.scorer.calls, [])

    def test_the_checks_also_run_when_rebuilding_changes_the_staged_file(self) -> None:
        c = self.build()
        staged = c.stage / "EPM" / "B6" / "R8I2.csv"
        staged.write_bytes(staged.read_bytes().replace(b"\r\n0,", b"\r\n0,0", 1))  # an edited staged file
        args = self.base(c, "--apply", "--skip-analysis", "--rebuild")
        before = tree_state(self.tmp)
        code, out = run_tool(args)
        self.assertNotEqual(code, 0)
        self.assertIn("is not 300 +/- 2 s", out)                 # the length check is not only for new exports
        self.assertIn("copies    body centre against 1 staged EPM B6 files", out)
        self.assertEqual(tree_state(self.tmp), before)
        code, out = run_tool(args + ["--accept-length"])
        self.assertEqual(code, 0, out)

    def test_a_second_retrack_the_same_day_keeps_both_old_versions(self) -> None:
        c = self.build()
        base = self.base(c, "--accept-length", "--apply", "--skip-analysis")
        self.add_new_export(c, seed=99)
        self.assertEqual(run_tool(base)[0], 0)
        sleap = c.raw / "B6" / "EPM" / "SLEAP"
        first = (sleap / "formatted" / "R8I2_EPM.merged_locs_formatted.csv").read_bytes()
        for kind in ("geom", "animal"):                          # a corrected export under the same name
            write_h5(sleap / kind / f"E9_B6_EPM-{kind}.019_E9_B6_R8I2_EPM.analysis.h5",
                     geom_xy("EPM", seed=8) if kind == "geom" else animal_xy("EPM", 98), stored("EPM", kind))
        code, out = run_tool(base)                               # in place but not rebuilding: not taken as new
        self.assertNotEqual(code, 0)
        self.assertIn("export in place does not rebuild its locs file", out)
        self.assertIn("add --rebuild", out)
        code, out = run_tool(base + ["--rebuild"])
        self.assertEqual(code, 0, out)
        old = sleap / "formatted" / "old"
        self.assertEqual(sorted(p.name for p in old.iterdir()),
                         ["R8I2_EPM.merged_locs_formatted.csv", "R8I2_EPM.merged_locs_formatted__2.csv"])
        self.assertEqual((old / "R8I2_EPM.merged_locs_formatted__2.csv").read_bytes(), first)
        self.assertTrue((c.stage / "EPM" / f"_superseded_{TODAY}" / "B6" / "R8I2__2.csv").exists())

    def test_an_export_older_than_the_chain_is_not_taken_as_new(self) -> None:
        c = self.build()
        self.add_new_export(c)
        for kind in ("geom", "animal"):
            p = c.raw / "B6" / "EPM" / "SLEAP" / kind / f"E9_B6_EPM-{kind}.019_E9_B6_R8I2_EPM.analysis.h5"
            os.utime(p, (1.6e9, 1.6e9))                         # 2020: older than the chain's locs files
        code, out = run_tool(self.base(c, "--accept-length"))
        self.assertEqual(code, 0, out)
        self.assertNotIn(" NEW ", out)
        self.assertIn("geom export in place does not rebuild its locs file", out)

    def test_a_failed_scoring_run_is_finished_with_rebuild(self) -> None:
        c = self.build()
        self.add_new_export(c)
        inner = self.scorer
        rr.SCORER = lambda *a: subprocess.CompletedProcess([], 1, "", "R error")
        base = self.base(c, "--accept-length")
        code, out = run_tool(base + ["--apply", "--baseline"])
        self.assertNotEqual(code, 0)
        self.assertIn("--rebuild", out)
        staged_old = (c.stage / "EPM" / "B6" / "R8I2.csv").read_bytes()
        rr.SCORER = inner
        code, out = run_tool(base)
        self.assertIn("--apply --rebuild would carry it out", out)
        code, out = run_tool(base + ["--apply", "--baseline", "--rebuild"])
        self.assertEqual(code, 0, out)
        self.assertNotEqual((c.stage / "EPM" / "B6" / "R8I2.csv").read_bytes(), staged_old)
        self.assertIn("animals whose rows changed: ['R8I2']", out)

    def test_a_failed_after_run_leaves_the_scoring_due_and_the_next_run_scores(self) -> None:
        c = self.build()
        self.add_new_export(c)
        self.scorer.fail_on_call = 2                             # the baseline works, the after run fails
        base = self.base(c, "--accept-length")
        code, out = run_tool(base + ["--apply", "--baseline"])
        self.assertNotEqual(code, 0)
        staged = c.stage / "EPM" / "B6" / "R8I2.csv"
        installed = staged.read_bytes()
        self.scorer.fail_on_call = None
        code, out = run_tool(base)
        self.assertEqual(code, 0, out)
        self.assertIn("scoring of the installed file is still due", out)
        self.assertNotIn("--apply --rebuild would carry it out", out)
        before = tree_state(self.tmp)
        code, out = run_tool(base + ["--apply"])                 # the code version rules hold for this scoring too
        self.assertEqual(code, 3, out)
        self.assertEqual(tree_state(self.tmp), before)
        code, out = run_tool(base + ["--apply", "--baseline"])
        self.assertEqual(code, 0, out)
        self.assertIn(f"baseline  with {c.stage / 'EPM' / f'_superseded_{TODAY}' / 'B6' / 'R8I2.csv'}", out)
        self.assertIn("animals whose rows changed: ['R8I2']", out)
        self.assertEqual(staged.read_bytes(), installed)         # the installed file is back after the baseline run
        self.assertEqual(len(self.scorer.calls), 4)
        state = tree_state(self.tmp)
        code, out = run_tool(base + ["--apply", "--baseline"])
        self.assertIn("everything is already up to date", out)
        self.assertEqual(tree_state(self.tmp), state)

    def test_an_install_without_scoring_is_scored_by_the_next_run(self) -> None:
        rr.CODE_STATUS = lambda sc, spec: ("same", "same code")
        c = self.build()
        self.add_new_export(c)
        base = self.base(c, "--accept-length")
        self.assertEqual(run_tool(base + ["--apply", "--skip-analysis"])[0], 0)
        self.assertEqual(self.scorer.calls, [])
        code, out = run_tool(base + ["--apply"])
        self.assertEqual(code, 0, out)
        self.assertIn("scoring of the installed file is still due", out)
        self.assertEqual(len(self.scorer.calls), 1)
        self.assertIn("animals whose rows changed: ['R8I2']", out)
        code, out = run_tool(base + ["--apply"])
        self.assertIn("everything is already up to date", out)
        self.assertEqual(len(self.scorer.calls), 1)

    def test_rerun_analysis_scores_without_a_new_export(self) -> None:
        rr.CODE_STATUS = lambda sc, spec: ("same", "same code")
        c = self.build()
        code, out = run_tool(self.base(c, "--apply", "--rerun-analysis"))
        self.assertEqual(code, 0, out)
        self.assertEqual(len(self.scorer.calls), 1)
        self.assertIn("animals whose rows changed: none", out)

    def test_rerun_analysis_with_another_code_version_needs_a_baseline(self) -> None:
        c = self.build()
        code, out = run_tool(self.base(c, "--apply", "--rerun-analysis"))
        self.assertEqual(code, 3, out)
        self.assertEqual(self.scorer.calls, [])
        code, out = run_tool(self.base(c, "--apply", "--rerun-analysis", "--baseline"))
        self.assertEqual(code, 0, out)
        self.assertEqual(len(self.scorer.calls), 2)
        self.assertIn("animals whose rows changed: none", out)

    def test_scoring_from_another_checkout_uses_its_script_and_config(self) -> None:
        other = self.tmp / "checkout"
        (other / "02_SLEAPanalzyer").mkdir(parents=True)
        (other / "02_SLEAPanalzyer" / "DLCA_EPM v1.0.0.R").write_text("# old version")
        cfg = other / "validation" / "exp9_boris" / "config"
        cfg.mkdir(parents=True)
        (cfg / "epm_all.yaml").write_text((rr.CONFIG_DIR / "epm_all.yaml").read_text(encoding="utf-8") + "\n# old copy\n")
        rr.CODE_STATUS = lambda sc, spec: ("same", f"checkout {sc['repo']}") if sc["repo"] == other.resolve() else ("changed", "x")
        seen = []
        inner = self.scorer
        rr.SCORER = lambda rscript, script, config, cwd: (seen.append((script, cwd)), inner(rscript, script, config, cwd))[1]
        c = self.build()
        self.add_new_export(c)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--scoring-repo", str(other)))
        self.assertEqual(code, 0, out)
        self.assertEqual(seen, [(other.resolve() / "02_SLEAPanalzyer" / "DLCA_EPM v1.0.0.R", other.resolve())])
        self.assertIn("# old copy", inner.calls[0])

    def test_a_config_that_does_not_exist_is_refused(self) -> None:
        c = self.build()
        missing = self.tmp / "configs" / "epm_typo.yaml"
        code, out = run_tool(self.base(c, "--config", str(missing)))
        self.assertNotEqual(code, 0)
        self.assertIn(f"--config {os.path.abspath(missing)} does not exist", out)
        code, out = run_tool(self.base(c, "--config", "epm_typo_relative.yaml"))   # a relative path: printed in full
        self.assertNotEqual(code, 0)
        self.assertIn(f"--config {os.path.abspath('epm_typo_relative.yaml')} does not exist", out)

    def test_an_export_given_by_path_is_copied_in(self) -> None:
        c = self.build()
        elsewhere = self.tmp / "exports"
        g = c.export("EPM", "B6", "geom", "E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.analysis.h5", geom_xy("EPM", seed=5), where=elsewhere)
        a = c.export("EPM", "B6", "animal", "E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5", animal_xy("EPM", 99), where=elsewhere)
        code, out = run_tool(self.base(c, "--accept-length", "--geom-h5", str(g), "--animal-h5", str(a), "--apply", "--skip-analysis"))
        self.assertEqual(code, 0, out)
        geom = c.raw / "B6" / "EPM" / "SLEAP" / "geom"
        self.assertEqual(sorted(p.name for p in geom.glob("*R8I2*")),
                         ["E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.analysis.h5", "E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.analysis_locs.csv"])
        self.assertTrue(g.exists())                             # the source is copied, not moved
        self.assertEqual(self.scorer.calls, [])

    def test_an_export_named_otherwise_needs_accept_name_and_is_copied_in_under_a_name_the_tool_finds(self) -> None:
        c = self.build()
        elsewhere = self.tmp / "exports"
        g = c.export("EPM", "B6", "geom", "arena_new.analysis.h5", geom_xy("EPM", seed=5), video="//share/E9_B6_R8I2_EPM.mp4", where=elsewhere)
        a = c.export("EPM", "B6", "animal", "mouse_new.analysis.h5", animal_xy("EPM", 99), video="//share/E9_B6_R8I2_EPM.mp4", where=elsewhere)
        args = self.base(c, "--accept-length", "--geom-h5", str(g), "--animal-h5", str(a), "--apply", "--skip-analysis")
        before = tree_state(self.tmp)
        code, out = run_tool(args)
        self.assertNotEqual(code, 0)
        self.assertIn("file name arena_new.analysis.h5 does not name R8I2", out)
        self.assertEqual(tree_state(self.tmp), before)
        code, out = run_tool(args + ["--accept-name"])
        self.assertEqual(code, 0, out)
        geom = c.raw / "B6" / "EPM" / "SLEAP" / "geom"
        self.assertTrue((geom / "arena_new_R8I2.analysis.h5").exists())
        self.assertTrue((geom / "arena_new_R8I2.analysis_locs.csv").exists())
        state = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--skip-analysis"))   # found again by name
        self.assertEqual(code, 0, out)
        self.assertIn("in place: geom arena_new_R8I2.analysis.h5", out)
        self.assertEqual(tree_state(self.tmp), state)

    def test_an_export_of_another_animal_is_refused(self) -> None:
        c = self.build()
        self.add_new_export(c, video="//share/cropped/E9_B6_A2C8_EPM.mp4")
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--skip-analysis"))
        self.assertNotEqual(code, 0)
        self.assertIn("video_path (//share/cropped/E9_B6_A2C8_EPM.mp4) does not name R8I2", out)
        self.assertIn("export names A2C8", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_a_copy_of_another_animals_trial_is_refused(self) -> None:
        c = self.build()
        c.export("EPM", "B6", "geom", "E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.analysis.h5", geom_xy("EPM", seed=5))
        c.export("EPM", "B6", "animal", "E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5", shifted(animal_xy("EPM", 10), 30))
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--skip-analysis"))
        self.assertNotEqual(code, 0)
        self.assertIn("the new track coincides with A2C8 at an offset of 30 frames", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_only_exports_and_locs_files_are_archived_and_stray_exports_are_listed(self) -> None:
        c = self.build()
        geom = c.raw / "B6" / "EPM" / "SLEAP" / "geom"
        (geom / "E9_B6_EPM-geom.013_E9_B6_R8I2_EPM.notes.txt").write_text("a note")
        pred = c.raw / "B6" / "EPM" / "SLEAP" / "predictions"
        write_h5(pred / "E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5", animal_xy("EPM", 99), stored("EPM", "animal"))
        code, out = run_tool(self.base(c))
        self.assertIn(f"export not in geom/ or animal/ (the tool does not read it): {pred / 'E9_B6_EPM-animal.019_E9_B6_R8I2_EPM.analysis.h5'}", out)
        self.assertIn("newer than the chain's locs files", out)
        self.add_new_export(c)
        code, out = run_tool(self.base(c, "--accept-length"))     # listed also when a new export is found
        self.assertEqual(out.count("export not in geom/ or animal/"), 1, out)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--skip-analysis"))
        self.assertEqual(code, 0, out)
        self.assertIn(f"left      {geom / 'E9_B6_EPM-geom.013_E9_B6_R8I2_EPM.notes.txt'}", out)
        self.assertTrue((geom / "E9_B6_EPM-geom.013_E9_B6_R8I2_EPM.notes.txt").exists())

    def test_apply_without_rscript_needs_the_python_neck_copy_option(self) -> None:
        rr.NECK_COPIER = None                                    # the real route: 08's copy_with_neck through Rscript
        c = self.build()
        self.add_new_export(c)
        code, out = run_tool(self.base(c, "--accept-length"))
        self.assertEqual(code, 0, out)
        self.assertIn("staged copy written in Python", out)
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--skip-analysis"))
        self.assertNotEqual(code, 0)
        self.assertIn("--python-neck-copy", out)
        self.assertEqual(tree_state(self.tmp), before)
        code, out = run_tool(self.base(c, "--accept-length", "--apply", "--skip-analysis", "--python-neck-copy"))
        self.assertEqual(code, 0, out)


class EpmB3Test(Stubs):
    """Batch 3 labels spine1 'neck' in its locs files and names the formatted files E9_SIS_B3_EPM_B3_<CODE>."""

    def test_neck_labels_and_the_batch_names_are_kept(self) -> None:
        c = Chain(self.tmp)
        for i, code in enumerate(["B3N2", "D9Q7"]):
            stem = f"{i:03d}_E9_SIS_B3_{code}"
            g = c.export("EPM", "B3", "geom", f"E9_SIS_B3_EPM-geom.{stem}.analysis.h5", geom_xy("EPM", seed=i))
            a = c.export("EPM", "B3", "animal", f"E9_SIS_B3_EPM-animal.{stem}.analysis.h5", animal_xy("EPM", seed=20 + i))
            c.finish("EPM", "B3", code, "", c.locs("EPM", "B3", "geom", g), c.locs("EPM", "B3", "animal", a),
                     f"E9_SIS_B3_EPM-geom.{stem}.merged_locs.csv", f"E9_SIS_B3_EPM_B3_{code}_formatted.csv")
        sleap = c.raw / "B3" / "EPM" / "SLEAP"
        base = c.args("--batch", "B3", "--assay", "EPM", "--code", "D9Q7", "--output-root", str(c.out), "--accept-length")
        code, out = run_tool(base)
        self.assertEqual(code, 0, out)
        self.assertIn("rebuild every file in the chain", out)    # the neck labels reproduce the files in place
        c.export("EPM", "B3", "geom", "E9_SIS_B3_EPM-geom.020_E9_SIS_B3_D9Q7.analysis.h5", geom_xy("EPM", seed=7))
        c.export("EPM", "B3", "animal", "E9_SIS_B3_EPM-animal.020_E9_SIS_B3_D9Q7.analysis.h5", animal_xy("EPM", 77))
        code, out = run_tool(base + ["--apply", "--skip-analysis"])
        self.assertEqual(code, 0, out)
        locs = pd.read_csv(sleap / "animal" / "E9_SIS_B3_EPM-animal.020_E9_SIS_B3_D9Q7.analysis_locs.csv")
        self.assertEqual(list(locs.columns[12:14]), ["neck_x", "neck_y"])
        self.assertEqual(sorted(p.name for p in (sleap / "formatted").glob("*.csv")),
                         ["E9_SIS_B3_EPM_B3_B3N2_formatted.csv", "E9_SIS_B3_EPM_B3_D9Q7_formatted.csv"])
        self.assertTrue((sleap / "formatted" / "old" / "E9_SIS_B3_EPM_B3_D9Q7_formatted.csv").exists())
        self.assertTrue((sleap / "E9_SIS_B3_EPM-geom.020_E9_SIS_B3_D9Q7.merged_locs.csv").exists())
        parts, _ = rr.read_dlc(c.stage / "EPM" / "B3" / "D9Q7.csv")
        self.assertEqual(parts.count("neck"), 3)
        self.assertNotIn("spine1", parts)                         # staged as the batch had it: neck, no alias added
        self.assertIn("copied (had neck)", out)


# --------------------------------------------------------------------------- #
# SocP
# --------------------------------------------------------------------------- #

class SocpTest(Stubs):
    def build(self) -> Chain:
        c = Chain(self.tmp)
        for i, (code, phase, seed) in enumerate([("G6K2", "S1", 1), ("G6K2", "S2", 2), ("A0S5", "S2", 3)]):
            stem = f"charge4.{i:03d}_E9_SIS_B2_{code}_{phase}"
            g = c.export("SocP", "B2", "geom", f"SocP-geom-{stem}.analysis.h5", geom_xy("SocP", seed=seed))
            a = c.export("SocP", "B2", "animal", f"SocP-animal-{stem}.analysis.h5", animal_xy("SocP", seed))
            c.finish("SocP", "B2", code, phase, c.locs("SocP", "B2", "geom", g), c.locs("SocP", "B2", "animal", a),
                     f"B2_{code}_{phase}.merged_locs.csv", f"{code}_{phase}.merged_locs_formatted.csv")
        for ph in ("S1", "S2"):
            side = "Code\tNovelLoc\nG6K2\tR\nA0S5\tL\n"
            (c.raw / "B2" / "SocP" / f"novelLoc{ph}.txt").write_text(side)
            (c.stage / "SocP" / "B2" / ph / f"novelLoc{ph}.txt").write_text(side)
        cfg = self.tmp / "socp_cfg.yaml"
        cfg.write_text(f"input_dir: {(c.stage / 'SocP').as_posix()}\noutput_dir: {(c.out / 'SocP').as_posix()}\n"
                       "batches: [B2]\nphases: [S1, S2]\n")
        self.scorer(None, Path("DLCA_SocP v.1.1.0.R"), cfg, None)
        self.scorer.calls.clear()
        return c

    def args(self, c: Chain, *extra: str) -> list[str]:
        return c.args("--batch", "B2", "--assay", "SocP", "--code", "G6K2", "--phase", "S2", "--output-root", str(c.out),
                      "--accept-length", *extra)

    def new_export(self, c: Chain, animal: np.ndarray) -> None:
        c.export("SocP", "B2", "geom", "SocP-geom-retrack.000_E9_SIS_B2_G6K2_S2.analysis.h5", geom_xy("SocP", seed=7))
        c.export("SocP", "B2", "animal", "SocP-animal-retrack.000_E9_SIS_B2_G6K2_S2.analysis.h5", animal)

    def test_phase_folders_side_file_and_archives(self) -> None:
        c = self.build()
        self.new_export(c, animal_xy("SocP", 77))
        code, out = run_tool(self.args(c))
        self.assertEqual(code, 0, out)
        self.assertIn("socp_all.yaml", out)
        self.assertIn("side file novelLocS2.txt: G6K2 R", out)
        self.assertIn("vs S1", out)
        self.assertIn("copies    body centre against 2 staged SocP B2 files (every phase)", out)
        self.assertIn("'^SocP [A-Za-z0-9]+ B2 '", out)
        code, out = run_tool(self.args(c, "--apply", "--baseline"))
        self.assertEqual(code, 0, out)
        sleap, arch = c.raw / "B2" / "SocP" / "SLEAP", f"_superseded_{TODAY}"
        self.assertTrue((sleap / "geom" / arch / "SocP-geom-charge4.001_E9_SIS_B2_G6K2_S2.analysis.h5").exists())
        self.assertTrue((sleap / "animal" / "SocP-animal-retrack.000_E9_SIS_B2_G6K2_S2.analysis_locs.csv").exists())
        self.assertTrue((sleap / arch / "B2_G6K2_S2.merged_locs.csv").exists())
        self.assertTrue((sleap / "B2_G6K2_S2.merged_locs.csv").exists())       # B2 keeps its merged name
        self.assertTrue((sleap / "formatted" / "old" / "G6K2_S2.merged_locs_formatted.csv").exists())
        self.assertTrue((c.stage / "SocP" / arch / "B2" / "S2" / "G6K2_S2.csv").exists())
        self.assertTrue((sleap / "geom" / "SocP-geom-charge4.000_E9_SIS_B2_G6K2_S1.analysis.h5").exists())   # S1 untouched
        self.assertIn("phases: [S2]", self.scorer.calls[0])
        self.assertIn("qc_decisions_file: " + (rr.CONFIG_DIR / "qc_decisions_socp.csv").resolve().as_posix(), self.scorer.calls[0])
        self.assertIn("animals whose rows changed: ['G6K2']", out)
        state = tree_state(self.tmp)
        code, out = run_tool(self.args(c, "--apply", "--baseline"))
        self.assertEqual((code, tree_state(self.tmp)), (0, state), out)

    def test_a_copy_of_the_other_phase_at_an_offset_is_refused(self) -> None:
        c = self.build()
        self.new_export(c, shifted(animal_xy("SocP", 1), 41))   # the phase-1 track, 41 frames on
        before = tree_state(self.tmp)
        code, out = run_tool(self.args(c, "--apply", "--baseline"))
        self.assertNotEqual(code, 0)
        self.assertIn("coincides with S1 at an offset of 41 frames", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_a_copy_of_another_animals_trial_is_refused(self) -> None:
        c = self.build()
        self.new_export(c, shifted(animal_xy("SocP", 3), 25))   # A0S5's phase-2 track, 25 frames on
        before = tree_state(self.tmp)
        code, out = run_tool(self.args(c, "--apply", "--baseline"))
        self.assertNotEqual(code, 0)
        self.assertIn("the new track coincides with S2/A0S5_S2 at an offset of 25 frames", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_a_byte_copy_of_another_staged_file_is_refused(self) -> None:
        c = self.build()
        sleap = c.raw / "B2" / "SocP" / "SLEAP"
        for kind in ("geom", "animal"):
            src = sleap / kind / f"SocP-{kind}-charge4.002_E9_SIS_B2_A0S5_S2.analysis.h5"
            (sleap / kind / f"SocP-{kind}-retrack.000_E9_SIS_B2_G6K2_S2.analysis.h5").write_bytes(src.read_bytes())
        code, out = run_tool(self.args(c))
        self.assertIn("byte-identical to", out)
        self.assertIn("A0S5_S2.csv", out)
        self.assertIn("export names A0S5", out)                  # its stored video_path is A0S5's

    def test_side_file_without_the_animal_is_refused(self) -> None:
        c = self.build()
        (c.stage / "SocP" / "B2" / "S2" / "novelLocS2.txt").write_text("Code\tNovelLoc\nA0S5\tL\n")
        self.new_export(c, animal_xy("SocP", 77))
        code, out = run_tool(self.args(c, "--apply", "--baseline"))
        self.assertNotEqual(code, 0)
        self.assertIn("does not list G6K2", out)

    def test_habituation_is_scored_with_the_habituation_config(self) -> None:
        c = self.build()
        code, out = run_tool(c.args("--batch", "B2", "--assay", "SocP", "--code", "G6K2", "--phase", "HAB",
                                    "--output-root", str(c.out)))
        self.assertEqual(code, 0, out)
        self.assertIn("socp_hab_all.yaml", out)
        self.assertIn("no new geom and animal export", out)

    def test_follow_ups_list_the_on_cup_rows_and_only_this_phases_list_rows(self) -> None:
        c = self.build()
        cfg_dir = self.tmp / "cfg"
        cfg_dir.mkdir()
        (cfg_dir / "socp_all.yaml").write_text((rr.CONFIG_DIR / "socp_all.yaml").read_text(encoding="utf-8"), encoding="utf-8")
        (cfg_dir / "qc_decisions_socp.csv").write_text(
            "file,target,start_s,end_s,decision,reason,reviewer\n"
            "G6K2_S2,socr,12.5,14.0,not_on_object,rearing at the wire,video check\n"
            "G6K2_S1,socl,1.0,2.0,not_on_object,other phase,video check\n"
            "A0S5_S2,socl,3.0,4.0,not_on_object,other animal,video check\n", encoding="utf-8")
        (c.qc / "retrack_list.csv").write_text(
            "animal,code,test,trial,problem\n"
            "OR131,G6K2,SocP,phase 2,phase-2 clip cut from phase 1\n"
            "OR131,G6K2,SocP,phase 1,a phase-1 problem\n"
            "OR131,G6K2,SocP,habituation,a habituation problem\n", encoding="utf-8")
        code, out = run_tool(self.args(c, "--config", str(cfg_dir / "socp_all.yaml")))
        self.assertEqual(code, 0, out)
        self.assertIn("qc_decisions_socp.csv", out)
        self.assertIn("has 1 row(s) for G6K2_S2; their times refer to the old clip", out)
        self.assertIn("line 2: socr 12.5-14.0 s -> not_on_object", out)
        self.assertIn("phase-2 clip cut from phase 1", out)
        self.assertNotIn("a phase-1 problem", out)
        self.assertNotIn("a habituation problem", out)


class SocpB1Test(Stubs):
    """Batch 1 names its files <CODE>_GX<clip>_<PH> and keeps replaced phase files in formatted/old/<PH>_old/."""

    def test_gx_names_and_the_phase_archive_folder(self) -> None:
        c = Chain(self.tmp)
        sleap = c.raw / "B1" / "SocP" / "SLEAP"
        for i, (code, gx, seed) in enumerate([("A1L9", "GX0786", 1), ("C2P4", "GX0790", 2)]):
            stem = f"charge1.{i:03d}_{code}_{gx}_S1"
            g = c.export("SocP", "B1", "geom", f"SocP-geom-{stem}.analysis.h5", geom_xy("SocP", seed=seed))
            a = c.export("SocP", "B1", "animal", f"SocP-animal-{stem}.analysis.h5", animal_xy("SocP", seed))
            c.finish("SocP", "B1", code, "S1", c.locs("SocP", "B1", "geom", g), c.locs("SocP", "B1", "animal", a),
                     f"{code}_{gx}_S1.merged_locs.csv", f"{code}_S1.merged_locs_formatted.csv")
            (sleap / f"{code}_S1.merged_locs.csv").write_bytes((sleap / f"{code}_{gx}_S1.merged_locs.csv").read_bytes())
        (sleap / "formatted" / "old" / "S1_old").mkdir(parents=True)
        side = "Code\tNovelLoc\nA1L9\tR\nC2P4\tL\n"
        (c.raw / "B1" / "SocP" / "novelLocS1.txt").write_text(side)
        (c.stage / "SocP" / "B1" / "S1" / "novelLocS1.txt").write_text(side)
        base = c.args("--batch", "B1", "--assay", "SocP", "--code", "A1L9", "--phase", "S1", "--output-root", str(c.out),
                      "--accept-length")
        code, out = run_tool(base)
        self.assertIn("rebuild every file in the chain", out)
        c.export("SocP", "B1", "geom", "SocP-geom-retrack.000_A1L9_GX0786_S1.analysis.h5", geom_xy("SocP", seed=5))
        c.export("SocP", "B1", "animal", "SocP-animal-retrack.000_A1L9_GX0786_S1.analysis.h5", animal_xy("SocP", 55))
        untouched = tree_state(sleap / "formatted" / "S1")
        code, out = run_tool(base + ["--apply", "--skip-analysis"])
        self.assertEqual(code, 0, out)
        arch = f"_superseded_{TODAY}"
        self.assertTrue((sleap / "formatted" / "old" / "S1_old" / "A1L9_S1.merged_locs_formatted.csv").exists())
        self.assertNotEqual((sleap / "formatted" / "S1" / "A1L9_S1.merged_locs_formatted.csv").read_bytes(),
                            (sleap / "formatted" / "old" / "S1_old" / "A1L9_S1.merged_locs_formatted.csv").read_bytes())
        self.assertEqual(tree_state(sleap / "formatted" / "S1")["C2P4_S1.merged_locs_formatted.csv"],
                         untouched["C2P4_S1.merged_locs_formatted.csv"])
        self.assertTrue((sleap / "geom" / arch / "SocP-geom-charge1.000_A1L9_GX0786_S1.analysis.h5").exists())
        self.assertTrue((sleap / "A1L9_GX0786_S1.merged_locs.csv").exists())     # the merged name in place is kept
        self.assertTrue((sleap / arch / "A1L9_S1.merged_locs.csv").exists())     # the second merged copy is archived
        self.assertTrue((c.stage / "SocP" / "B1" / "S1" / "A1L9_S1.csv").exists())
        state = tree_state(self.tmp)
        code, out = run_tool(base + ["--apply", "--skip-analysis"])
        self.assertEqual((code, tree_state(self.tmp)), (0, state), out)


class SocpB3Test(Stubs):
    """Batch 3 keeps a flat copy of each formatted file in formatted/ beside formatted/<PH>/."""

    def test_the_flat_copy_moves_to_old_with_the_phase_file(self) -> None:
        c = Chain(self.tmp)
        sleap = c.raw / "B3" / "SocP" / "SLEAP"
        for i, (code, seed) in enumerate([("B3N2", 1), ("D9Q7", 2)]):
            stem = f"{i:03d}_E9_SIS_B3_{code}_S2"
            g = c.export("SocP", "B3", "geom", f"E9_SIS_B3_SocP-geom.{stem}.analysis.h5", geom_xy("SocP", seed=seed))
            a = c.export("SocP", "B3", "animal", f"E9_SIS_B3_SocP-animal.{stem}.analysis.h5", animal_xy("SocP", seed))
            fmt = c.finish("SocP", "B3", code, "S2", c.locs("SocP", "B3", "geom", g), c.locs("SocP", "B3", "animal", a),
                           f"E9_SIS_B3_SocP-geom.{stem}.merged_locs.csv", f"{code}_S2.merged_locs_formatted.csv")
            (sleap / "formatted" / fmt.name).write_bytes(fmt.read_bytes())
            (sleap / f"B3_{code}_S2.analysis_locs.csv").write_bytes((sleap / f"E9_SIS_B3_SocP-geom.{stem}.merged_locs.csv").read_bytes())
        side = "Code\tNovelLoc\nB3N2\tR\nD9Q7\tL\n"
        (c.raw / "B3" / "SocP" / "novelLocS2.txt").write_text(side)
        (c.stage / "SocP" / "B3" / "S2" / "novelLocS2.txt").write_text(side)
        base = c.args("--batch", "B3", "--assay", "SocP", "--code", "D9Q7", "--phase", "S2", "--output-root", str(c.out),
                      "--accept-length")
        code, out = run_tool(base)
        self.assertIn("rebuild every file in the chain", out)    # the second merged copy and the flat copy are no change
        c.export("SocP", "B3", "geom", "E9_SIS_B3_SocP-geom.020_E9_SIS_B3_D9Q7_S2.analysis.h5", geom_xy("SocP", seed=8))
        c.export("SocP", "B3", "animal", "E9_SIS_B3_SocP-animal.020_E9_SIS_B3_D9Q7_S2.analysis.h5", animal_xy("SocP", 88))
        code, out = run_tool(base)
        self.assertEqual(code, 0, out)
        self.assertIn(f"move        {sleap / 'formatted' / 'D9Q7_S2.merged_locs_formatted.csv'}", out)
        code, out = run_tool(base + ["--apply", "--skip-analysis"])
        self.assertEqual(code, 0, out)
        self.assertEqual(sorted(p.name for p in (sleap / "formatted").glob("*.csv")), ["B3N2_S2.merged_locs_formatted.csv"])
        self.assertEqual(sorted(p.name for p in (sleap / "formatted" / "old").iterdir()),
                         ["D9Q7_S2.merged_locs_formatted.csv", "D9Q7_S2.merged_locs_formatted__2.csv"])
        self.assertTrue((sleap / "E9_SIS_B3_SocP-geom.020_E9_SIS_B3_D9Q7_S2.merged_locs.csv").exists())   # follows the export
        self.assertTrue((sleap / f"_superseded_{TODAY}" / "B3_D9Q7_S2.analysis_locs.csv").exists())
        self.assertTrue((sleap / f"_superseded_{TODAY}" / "E9_SIS_B3_SocP-geom.001_E9_SIS_B3_D9Q7_S2.merged_locs.csv").exists())


# --------------------------------------------------------------------------- #
# OFT
# --------------------------------------------------------------------------- #

class OftTest(Stubs):
    def build(self) -> Chain:
        c = Chain(self.tmp)
        pred = c.raw / "B1" / "OFT" / "SLEAP" / "predictions"
        for code, clip, seed in [("A1L9", "GX010796", 1), ("X1H0", "GX010804", 2)]:
            g = c.export("OFT", "B1", "geom", f"OFT_edges.0_{clip}_cropped.analysis.h5", geom_xy("OFT", seed=seed), where=pred)
            a = c.export("OFT", "B1", "animal", f"E9_B1_{clip}.analysis.h5", animal_xy("OFT", seed), where=pred,
                         video=f"//share/cropped/{code}_{clip}.mp4")
            gl = c.locs("OFT", "B1", "geom", g, f"E9_B1_OFT_{clip}.analysis_locs.csv")
            al = c.locs("OFT", "B1", "animal", a, f"E9_B1_OFT_{clip}.analysis_locs.csv")
            c.finish("OFT", "B1", code, "", gl, al, f"E9_B1_OFT_{clip}.merged_locs.csv", f"E9_B1_OFT_{code}.merged_locs_formatted.csv")
        self.scorer.oft(c.stage / "OFT", "B1", False)
        g = c.stage / "OFT" / "OFT_SLEAP_output_v1.2.0"
        (g / "only_B1.txt").unlink()
        pd.DataFrame([{"Batch": "B2", "Code": "ZZZZ"}]).to_csv(g / "all_batches_summary.csv", index=False)
        self.scorer.calls.clear()
        return c

    def new_export(self, c: Chain, seed: int = 55, video: str = "//share/cropped/X1H0_GX010804.mp4") -> None:
        c.export("OFT", "B1", "animal", "labels_retrack.000_X1H0_GX010804.analysis.h5", animal_xy("OFT", seed), video=video)

    def base(self, c: Chain, *extra: str) -> list[str]:
        return c.args("--batch", "B1", "--assay", "OFT", "--code", "X1H0", "--clip", "GX010804", "--reuse-geom",
                      "--accept-length", *extra)

    def test_b1_needs_the_clip(self) -> None:
        c = self.build()
        code, out = run_tool(c.args("--batch", "B1", "--assay", "OFT", "--code", "X1H0"))
        self.assertNotEqual(code, 0)
        self.assertIn("--clip", out)

    def test_b1_waits_for_the_animal_export_and_shows_the_files_by_clip(self) -> None:
        c = self.build()
        code, out = run_tool(c.args("--batch", "B1", "--assay", "OFT", "--code", "X1H0", "--clip", "GX010804", "--reuse-geom"))
        self.assertEqual(code, 0, out)
        self.assertIn("no new animal export", out)
        self.assertIn("E9_B1_OFT_GX010804.analysis_locs.csv", out)
        self.assertIn("E9_B1_OFT_X1H0.merged_locs_formatted.csv", out)
        self.assertIn("<labels file>.<NNN>_X1H0_GX010804.analysis.h5", out)
        self.assertIn("E9_B1_OFT_GX010804.merged_locs.csv formats to E9_B1_OFT_X1H0.merged_locs_formatted.csv byte for byte", out)
        self.assertIn("predictions", out)                        # the old exports there are listed
        self.assertIn("which staged file the current output was scored from is not recorded", out)

    def test_b1_new_animal_export_with_the_corners_in_place(self) -> None:
        c = self.build()
        self.new_export(c)
        g_before = tree_state(c.raw / "B1" / "OFT" / "SLEAP" / "geom")
        all_before = tree_state(c.stage / "OFT" / "OFT_SLEAP_output_v1.2.0")
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertEqual(code, 0, out)
        sleap, arch = c.raw / "B1" / "OFT" / "SLEAP", f"_superseded_{TODAY}"
        self.assertEqual(tree_state(sleap / "geom"), g_before)                  # corners reused, not touched
        self.assertTrue((sleap / "animal" / arch / "E9_B1_OFT_GX010804.analysis_locs.csv").exists())
        self.assertEqual(pd.read_csv(sleap / "animal" / "E9_B1_OFT_GX010804.analysis_locs.csv").columns[18], "center_x")
        self.assertTrue((sleap / "formatted" / "old" / "E9_B1_OFT_X1H0.merged_locs_formatted.csv").exists())
        parts, _ = rr.read_dlc(sleap / "formatted" / "E9_B1_OFT_X1H0.merged_locs_formatted.csv")
        self.assertIn("bodycentre", parts)
        self.assertNotIn("center", parts)
        self.assertTrue((c.stage / "OFT" / arch / "B1" / "X1H0.csv").exists())
        self.assertEqual(tree_state(c.stage / "OFT" / "OFT_SLEAP_output_v1.2.0"), all_before)   # put back whole
        self.assertIn("animals whose rows changed: ['X1H0']", out)              # the cohort z column is not compared
        self.assertEqual(len(self.scorer.calls), 2)
        marker = rr.read_marker(c.stage / "OFT" / "B1" / "output_v1.2.0" / rr.MARKER)
        self.assertEqual(marker["X1H0.csv"], md5(c.stage / "OFT" / "B1" / "X1H0.csv"))
        state = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertEqual((code, tree_state(self.tmp)), (0, state), out)
        code, out = run_tool(self.base(c))           # installed: not waiting, the old exports sit in predictions/
        self.assertNotIn("WAITING FOR THE NEW EXPORT", out)
        self.assertIn("if you re-track again", out)

    def test_an_export_saved_next_to_the_labels_file_is_listed(self) -> None:
        c = self.build()
        projects = c.raw / "B1" / "OFT" / "SLEAP" / "SLEAP_projects"
        write_h5(projects / "OFT_animals_b1.005_X1H0_GX010804.analysis.h5", animal_xy("OFT", 55), stored("OFT", "animal"))
        code, out = run_tool(c.args("--batch", "B1", "--assay", "OFT", "--code", "X1H0", "--clip", "GX010804", "--reuse-geom"))
        self.assertEqual(code, 0, out)
        self.assertIn(str(projects / "OFT_animals_b1.005_X1H0_GX010804.analysis.h5"), out)

    def test_a_global_folder_the_run_creates_is_not_left_behind(self) -> None:
        c = self.build()
        g = c.stage / "OFT" / "OFT_SLEAP_output_v1.2.0"
        for p in sorted(g.rglob("*"), reverse=True):
            p.unlink() if p.is_file() else p.rmdir()
        g.rmdir()
        self.new_export(c)
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertEqual(code, 0, out)
        self.assertFalse(g.exists(), out)
        self.assertIn("there was none before", out)

    def test_a_wrong_clip_is_refused_before_anything_changes(self) -> None:
        c = self.build()
        self.new_export(c)
        args = c.args("--batch", "B1", "--assay", "OFT", "--code", "X1H0", "--clip", "GX010796", "--reuse-geom",
                      "--accept-length", "--apply", "--baseline")
        before = tree_state(self.tmp)
        code, out = run_tool(args)
        self.assertNotEqual(code, 0)
        self.assertIn("is another animal's clip", out)
        self.assertIn("video_path (//share/cropped/X1H0_GX010804.mp4) does not name GX010796", out)
        self.assertEqual(tree_state(self.tmp), before)
        self.assertEqual(self.scorer.calls, [])

    def test_an_export_of_another_animals_clip_is_refused(self) -> None:
        c = self.build()
        self.new_export(c, video="//share/cropped/A1L9_GX010796.mp4")
        before = tree_state(self.tmp)
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertNotEqual(code, 0)
        self.assertIn("export names A1L9", out)
        self.assertEqual(tree_state(self.tmp), before)

    def test_failures_of_an_earlier_run_are_set_aside_and_a_failure_now_leaves_the_scoring_due(self) -> None:
        c = self.build()
        g = c.stage / "OFT" / "OFT_SLEAP_output_v1.2.0"
        pd.DataFrame([{"file": "X1H0.csv", "error": "an earlier run"}]).to_csv(g / "failures.csv", index=False)
        all_before = tree_state(g)
        self.new_export(c)
        self.scorer.failures_on_call, self.scorer.failures_code = 2, "X1H0"
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertNotEqual(code, 0)
        self.assertIn("moved", out)                               # the earlier failures.csv did not stop the baseline
        self.assertIn("written by this run", out)
        self.assertEqual(tree_state(g), all_before)               # put back whole, failures.csv of before included
        self.assertEqual(len(self.scorer.calls), 2)
        marker = rr.read_marker(c.stage / "OFT" / "B1" / "output_v1.2.0" / rr.MARKER)
        self.assertNotEqual(marker["X1H0.csv"], md5(c.stage / "OFT" / "B1" / "X1H0.csv"))
        self.scorer.failures_on_call = None
        code, out = run_tool(self.base(c))
        self.assertIn("scoring of the installed file is still due", out)
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertEqual(code, 0, out)
        self.assertIn("animals whose rows changed: ['X1H0']", out)
        self.assertEqual(tree_state(g), all_before)               # a second restore the same day leaves no B1-only file
        self.assertTrue((c.stage / f"_superseded_{TODAY}_X1H0" / "OFT_all_batches_written_by_B1_after_run__2").is_dir())
        code, out = run_tool(self.base(c, "--apply", "--baseline"))
        self.assertIn("everything is already up to date", out)

    def test_output_root_is_refused_for_the_open_field(self) -> None:
        c = self.build()
        code, out = run_tool(c.args("--batch", "B1", "--assay", "OFT", "--code", "X1H0", "--clip", "GX010804",
                                    "--output-root", str(c.out)))
        self.assertNotEqual(code, 0)
        self.assertIn("--stage-root", out)


class OftB2Test(Stubs):
    """Batch 2: E9_SIS_B2_<CODE> names for geom, animal and merged; 'center' in the locs; formatted at 10 digits."""

    def build(self) -> Chain:
        c = Chain(self.tmp)
        for code, seed in [("A0S5", 1), ("F9M6", 2)]:
            g = c.export("OFT", "B2", "geom", f"E9_SIS_B2_{code}.analysis.h5", geom_xy("OFT", seed=seed))
            a = c.export("OFT", "B2", "animal", f"E9_SIS_B2_{code}.analysis.h5", animal_xy("OFT", seed))
            c.finish("OFT", "B2", code, "", c.locs("OFT", "B2", "geom", g), c.locs("OFT", "B2", "animal", a),
                     f"E9_SIS_B2_{code}.merged_locs.csv", f"E9_SIS_B2_{code}.merged_locs_formatted.csv")
        self.scorer.oft(c.stage / "OFT", "B2", False)
        self.scorer.calls.clear()
        return c

    def test_the_chain_in_place_is_reproduced_and_a_new_export_follows_the_batch_names(self) -> None:
        rr.CODE_STATUS = lambda sc, spec: ("same", "same code")
        c = self.build()
        sleap = c.raw / "B2" / "OFT" / "SLEAP"
        line = (sleap / "formatted" / "E9_SIS_B2_F9M6.merged_locs_formatted.csv").read_text().splitlines()[3]
        self.assertTrue(all(len(v.replace("-", "").replace(".", "").lstrip("0")) <= 10 for v in line.split(",")))
        base = c.args("--batch", "B2", "--assay", "OFT", "--code", "F9M6", "--accept-length")
        code, out = run_tool(base)
        self.assertEqual(code, 0, out)
        self.assertIn("rebuild every file in the chain", out)    # 10 significant digits reproduce the files in place
        self.assertNotIn("  replace     ", out)
        self.assertIn("e.g. <labels file>.<NNN>_E9_SIS_B2_F9M6.analysis.h5", out)
        c.export("OFT", "B2", "geom", "OFT-geom.000_E9_SIS_B2_F9M6.analysis.h5", geom_xy("OFT", seed=6))
        c.export("OFT", "B2", "animal", "OFT-animal.000_E9_SIS_B2_F9M6.analysis.h5", animal_xy("OFT", 66))
        code, out = run_tool(base + ["--apply"])
        self.assertEqual(code, 0, out)
        self.assertTrue((sleap / "OFT-geom.000_E9_SIS_B2_F9M6.merged_locs.csv").exists())
        self.assertTrue((sleap / f"_superseded_{TODAY}" / "E9_SIS_B2_F9M6.merged_locs.csv").exists())
        self.assertTrue((sleap / "formatted" / "E9_SIS_B2_F9M6.merged_locs_formatted.csv").exists())
        self.assertTrue((sleap / "formatted" / "old" / "E9_SIS_B2_F9M6.merged_locs_formatted.csv").exists())
        self.assertEqual(pd.read_csv(sleap / "animal" / "OFT-animal.000_E9_SIS_B2_F9M6.analysis_locs.csv").columns[18], "center_x")
        self.assertIn("animals whose rows changed: ['F9M6']", out)
        state = tree_state(self.tmp)
        code, out = run_tool(base + ["--apply"])
        self.assertEqual((code, tree_state(self.tmp)), (0, state), out)


# --------------------------------------------------------------------------- #
# NOR keeps its command and files
# --------------------------------------------------------------------------- #

class NorTest(Stubs):
    def build(self) -> Chain:
        c = Chain(self.tmp)
        for code, seed in [("H1E2", 1), ("F8O9", 2)]:
            c.export("NOR", "B5", "geom", f"E9_B5_NOR-geom.005_E9_B5_{code}_NOV.analysis.h5", geom_xy("NOR", seed=seed))
            c.export("NOR", "B5", "animal", f"E9_B5_NOR-animal.005_E9_B5_{code}_NOV.analysis.h5", animal_xy("NOR", seed))
        return c

    def test_plan_then_apply_writes_the_nor_names(self) -> None:
        c = self.build()
        args = c.args("--batch", "B5", "--assay", "NOR", "--code", "H1E2", "--phase", "NOV", "--output-root", str(c.out))
        before = tree_state(self.tmp)
        code, out = run_tool(args)
        self.assertEqual(code, 0, out)
        self.assertEqual(tree_state(self.tmp), before)
        self.assertIn("then      back up", out)
        code, out = run_tool(args + ["--apply"])
        self.assertEqual(code, 0, out)
        sleap = c.raw / "B5" / "NOR" / "SLEAP"
        self.assertTrue((sleap / "geom" / "E9_B5_NOR-geom.005_E9_B5_H1E2_NOV.analysis_locs.csv").exists())
        self.assertTrue((sleap / "animal" / "E9_B5_NOR-animal.005_E9_B5_H1E2_NOV.analysis_locs.csv").exists())
        self.assertTrue((sleap / "formatted" / "H1E2_NOV.merged_locs_formatted.csv").exists())
        self.assertTrue((c.stage / "NOR" / "B5" / "H1E2_NOV.csv").exists())
        self.assertEqual(len(self.scorer.calls), 1)
        self.assertIn("batches: [B5]", self.scorer.calls[0])
        self.assertIn(f"input_dir: {(c.stage / 'NOR').as_posix()}", self.scorer.calls[0])
        self.assertIn("H1E2: validFraction", out)
        state = tree_state(self.tmp)
        code, out = run_tool(args + ["--apply"])
        self.assertEqual((code, tree_state(self.tmp)), (0, state), out)
        self.assertIn("everything is already up to date", out)

    def test_nor_refuses_the_options_of_the_other_assays(self) -> None:
        c = self.build()
        for option in ("--baseline", "--accept-name", "--accept-clip", "--python-neck-copy"):
            code, out = run_tool(c.args("--batch", "B5", "--assay", "NOR", "--code", "H1E2", "--phase", "NOV", option))
            self.assertNotEqual(code, 0, option)


class LayoutTest(unittest.TestCase):
    def test_resolve_code_follows_the_staging_script(self) -> None:
        codes = ["R8I2", "F2I1", "M8P3", "A0S5"]
        self.assertEqual(rr.resolve_code("R8I2_EPM.merged_locs_formatted", codes, "R8I2"), "R8I2")
        self.assertEqual(rr.resolve_code("E9_B2_F2l1_x", codes, "F2I1"), "F2I1")          # lookalike l for I
        self.assertEqual(rr.resolve_code("E9_M8P1_x", codes, "M8P3"), "M8P3")             # alias
        self.assertIsNone(rr.resolve_code("R8I2_A0S5", codes, "R8I2"))                    # ambiguous

    def test_merged_names_follow_the_new_export_only_where_the_batch_names_them_after_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lay = rr.Layout("EPM", "B6", "R8I2", "", None, Path(tmp), Path(tmp) / "s")
            lay.sleap.mkdir(parents=True)
            (lay.sleap / "E9_B6_EPM-geom.013_E9_B6_R8I2_EPM.merged_locs.csv").write_text("x")
            self.assertEqual(lay.merged_target("E9_B6_EPM-geom.019_E9_B6_R8I2_EPM").name,
                             "E9_B6_EPM-geom.019_E9_B6_R8I2_EPM.merged_locs.csv")
            lay = rr.Layout("SocP", "B4", "I4O9", "S2", None, Path(tmp), Path(tmp) / "s")
            lay.sleap.mkdir(parents=True)
            (lay.sleap / "B4_I4O9_S2.analysis_locs.csv").write_text("x")
            self.assertEqual(lay.merged_target("E9_B4_SocP-geom_2.007_E9_SIS_B4_I4O9_S2").name, "B4_I4O9_S2.analysis_locs.csv")

    def test_install_names_carry_the_recording(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            lay = rr.Layout("SocP", "B2", "G6K2", "S2", None, Path(tmp), Path(tmp) / "s")
            self.assertEqual(lay.install_name(Path("SocP-animal.003_E9_SIS_B2_G6K2_S2.analysis.h5")),
                             "SocP-animal.003_E9_SIS_B2_G6K2_S2.analysis.h5")
            self.assertEqual(lay.install_name(Path("mouse_old.h5")), "mouse-old_G6K2_S2.analysis.h5")
            lay = rr.Layout("OFT", "B1", "X1H0", "", "GX010804", Path(tmp), Path(tmp) / "s")
            self.assertEqual(lay.install_name(Path("E9_B1_GX010804.analysis.h5")), "E9_B1_GX010804.analysis.h5")
            self.assertEqual(lay.install_name(Path("new.analysis.h5")), "new_X1H0_GX010804.analysis.h5")

    def test_closest_offset_finds_a_copy_at_an_offset(self) -> None:
        track = animal_xy("SocP", 4, 3000)[:, 9] / 18
        for lag in (0, 41, -700, 1800):
            a = shifted(track, lag) if lag > 0 else (np.concatenate([track[:-lag], track[:lag]]) if lag < 0 else track)
            dist, found = rr.closest_offset(a, track)
            self.assertEqual((round(dist, 9), found), (0.0, lag))
        dist, _ = rr.closest_offset(animal_xy("SocP", 5, 3000)[:, 9] / 18, track)
        self.assertGreater(dist, rr.DUPLICATE_CM)


if __name__ == "__main__":
    unittest.main()
