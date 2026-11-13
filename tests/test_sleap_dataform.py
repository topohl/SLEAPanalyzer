"""Focused tests for the SLEAP-to-DLCAnalyzer CSV formatter."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


MODULE_DIR = Path(__file__).resolve().parents[1] / "01_SLEAPcoords"
sys.path.insert(0, str(MODULE_DIR))

import exp9_paths  # noqa: E402
import sleap_dataform as sd  # noqa: E402


class RootAndOutputPathTest(unittest.TestCase):
    """The tracking root holds <batch>/<assay>/; only SocP puts its phases in subfolders of formatted/."""

    def test_the_default_root_is_the_tracking_root(self) -> None:
        self.assertEqual(sd.BEHAVIOR_ROOT, exp9_paths.TRACKING_ROOT)
        self.assertEqual(exp9_paths.tracking_dir("B5", "NOR", Path("T")), Path("T") / "B5" / "NOR")
        self.assertEqual(exp9_paths.records_dir("B2", "SocP", Path("R")), Path("R") / "B2" / "SocP")

    def test_socp_in_phase_folders_nor_and_unphased_flat(self) -> None:
        f = Path("formatted")
        self.assertEqual(sd.output_path(f, "G6K2", "S2", "SocP"), f / "S2" / "G6K2_S2.merged_locs_formatted.csv")
        self.assertEqual(sd.output_path(f, "H1E2", "NOV", "NOR"), f / "H1E2_NOV.merged_locs_formatted.csv")
        self.assertEqual(sd.output_path(f, "R8I2", "", "EPM"), f / "R8I2.merged_locs_formatted.csv")

    def test_one_merged_and_formatted_name_for_every_batch(self) -> None:
        self.assertEqual(sd.merged_name("R8I2"), "R8I2.merged_locs.csv")
        self.assertEqual(sd.merged_name("A1L9", "NOV"), "A1L9_NOV.merged_locs.csv")
        self.assertEqual(sd.formatted_name("X1H0"), "X1H0.merged_locs_formatted.csv")
        self.assertEqual(sd.formatted_name("G6K2", "S2"), "G6K2_S2.merged_locs_formatted.csv")


class ParseKeyTest(unittest.TestCase):
    """One (code, phase) per file name, however the batch named its exports."""

    NOR, SOCP = sd.PHASES["NOR"], sd.PHASES["SocP"]

    def test_names_of_every_batch(self) -> None:
        cases = [
            ("E9_B6_SocP-animal_3.004_E9_SIS_B6_D3U7_S1.analysis_locs.csv", self.SOCP, ("D3U7", "S1")),
            ("SocP-geom-charge1.000_A1L9_GX0786_HAB.analysis_locs.csv", self.SOCP, ("A1L9", "HAB")),      # B1 SocP clip
            ("NOR_B1_geom_predict.000_A1L9_GX010405_NOV.analysis_locs.csv", self.NOR, ("A1L9", "NOV")),   # B1 NOR geom
            ("A1L9_GX010405_NOV.analysis_locs.csv", self.NOR, ("A1L9", "NOV")),                           # B1 NOR animal
            ("E9_B5_NOR-geom.003_E9_B5_F8O9_NOR.analysis_locs.csv", self.NOR, ("F8O9", "NOV")),           # B5 legacy label
            ("E9_B5_NOR-animal.002_E9_B5_F8O9_HAB.analysis_locs.csv", self.NOR, ("F8O9", "HAB")),
            ("R803_GX010395_NOV.analysis_locs.csv", self.NOR, ("R8O3", "NOV")),                           # misspelt code
            ("SocP-geom-charge5.004_E9_SIS_B2_O6O6_S1.analysis_locs.csv", self.SOCP, ("Q6Q6", "S1")),     # B2 SocP
            ("E9_B4_SocP-animal_1.006_E9_SIS_B4_B7C2_S2.analysis_locs.csv", self.SOCP, ("B7C9", "S2")),   # B4 SocP
            ("E9_B2_SocP-geom_retrack_G6K2.000_E9_SIS_B2_G6K2_S1_1.analysis_locs.csv", self.SOCP, ("G6K2", "S1")),
            ("E9_B4_NOR-animal.035_E9_B4_B4Y7_NOV_1.analysis_locs.csv", self.NOR, ("B4Y7", "NOV")),
            ("E9_B1_EPM-geom.018_E9_B1_EPM_GX010035_M8P1.analysis_locs.csv", (), ("M8P3", "")),
            ("E9_SIS_B2_A0S5.analysis_locs.csv", (), ("A0S5", "")),
            ("E9_B2_EPM-geom.019_E9_SIS_B2_T4W6_1.analysis_locs.csv", (), ("T4W6", "")),
        ]
        for name, phases, want in cases:
            with self.subTest(name=name):
                self.assertEqual(sd.parse_key(name, phases), want)

    def test_unclear_names_are_refused(self) -> None:
        for name, phases in [("E9_B6_SocP_D3U7_NOR.analysis_locs.csv", self.SOCP),   # _NOR is a phase only for NOR
                             ("x_A1L9_S1.y_B2C3_S2.analysis_locs.csv", self.SOCP),    # two keys
                             ("A1L9_GX010405_NOVEL.analysis_locs.csv", self.NOR),
                             ("A1L9_NOV_old.analysis_locs.csv", self.NOR),
                             ("E9_B1_OFT_GX010379.analysis_locs.csv", ())]:          # named by clip: no code
            with self.subTest(name=name), self.assertRaises(ValueError):
                sd.parse_key(name, phases)


class ProcessOneTest(unittest.TestCase):
    def test_writes_dlc_header_body_and_crlf_with_runtime_pandas(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            geom_path = root / "geom.csv"
            animal_path = root / "animal.csv"
            dst = root / "formatted.csv"

            pd.DataFrame(
                {"tl_x": [1.0, 2.0], "tl_y": [3.0, 4.0]}
            ).to_csv(geom_path, index=False)
            pd.DataFrame(
                {"nose_x": [5.0, 6.0], "nose_y": [7.0, 8.0]}
            ).to_csv(animal_path, index=False)

            info = sd.process_one(geom_path, animal_path, dst)

            self.assertEqual(
                info, {"bodyparts": 2, "cols": 7, "frames": 2}
            )
            raw = dst.read_bytes()
            self.assertEqual(raw.count(b"\r\n"), 5)
            self.assertNotIn(b"\n", raw.replace(b"\r\n", b""))
            self.assertNotIn(b"\r", raw.replace(b"\r\n", b""))
            self.assertEqual(
                raw.decode("utf-8").splitlines(),
                [
                    "column_1,column_2,column_3,column_3_new,column_4,column_5,column_5_new",
                    "bodyparts,tl,tl,tl,nose,nose,nose",
                    "coords,x,y,likelihood,x,y,likelihood",
                    "0,1.0,3.0,1,5.0,7.0,1",
                    "1,2.0,4.0,1,6.0,8.0,1",
                ],
            )


class CodeFilterTest(unittest.TestCase):
    """B5 NOR: the legacy ``_NOR`` label is the novel trial; beside a ``_NOV`` file of the animal it is set aside.
    --code formats one animal without parsing the others."""

    def make_tree(self, root: Path, names=("H1E2_NOV", "H1E2_HAB", "F8O9_NOR", "M2L3_NOR", "M2L3_NOV")) -> None:
        for sub, cols in (("geom", {"tl_x": [1.0, 2.0], "tl_y": [3.0, 4.0]}),
                          ("animal", {"nose_x": [5.0, 6.0], "nose_y": [7.0, 8.0]})):
            d = root / "B5" / "NOR" / sub
            d.mkdir(parents=True, exist_ok=True)
            for i, name in enumerate(names):
                shift = 100.0 if name.endswith("_NOV") else 0.0
                pd.DataFrame({k: [v + shift + i for v in vs] for k, vs in cols.items()}).to_csv(
                    d / f"E9_B5_NOR-{sub}.{i:03d}_E9_B5_{name}.analysis_locs.csv", index=False)

    def test_the_legacy_label_is_the_novel_trial(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_tree(root)
            report = sd.run("B5", "NOR", root=root, verbose=False)
            self.assertEqual(sorted(zip(report.code, report.phase)),
                             [("F8O9", "NOV"), ("H1E2", "HAB"), ("H1E2", "NOV"), ("M2L3", "NOV")])
            fmt = root / "B5" / "NOR" / "formatted"
            self.assertEqual(sorted(p.name for p in fmt.glob("*.csv")),
                             ["F8O9_NOV.merged_locs_formatted.csv", "H1E2_HAB.merged_locs_formatted.csv",
                              "H1E2_NOV.merged_locs_formatted.csv", "M2L3_NOV.merged_locs_formatted.csv"])
            rows = (fmt / "M2L3_NOV.merged_locs_formatted.csv").read_text().splitlines()
            self.assertEqual(rows[3], "0,105.0,107.0,1,109.0,111.0,1")      # from the _NOV files, not the _NOR ones

    def test_the_set_aside_files_are_listed_and_two_files_of_one_label_stop_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_tree(root)
            notes: list[str] = []
            idx = sd.index_folder(root / "B5" / "NOR" / "geom", sd.PHASES["NOR"], notes=notes)
            self.assertEqual(idx[("M2L3", "NOV")].name, "E9_B5_NOR-geom.004_E9_B5_M2L3_NOV.analysis_locs.csv")
            self.assertEqual(len(notes), 1)
            self.assertIn("E9_B5_NOR-geom.003_E9_B5_M2L3_NOR.analysis_locs.csv not used", notes[0])
        with tempfile.TemporaryDirectory() as tmp:
            self.make_tree(Path(tmp), ("F8O9_NOR", "F8O9_HAB", "F8O9_NOR"))
            with self.assertRaises(ValueError):
                sd.run("B5", "NOR", root=Path(tmp), dry_run=True, verbose=False)

    def test_with_code_only_that_animal_is_formatted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_tree(root)
            report = sd.run("B5", "NOR", root=root, code="H1E2", verbose=False)
            self.assertEqual(sorted(zip(report.code, report.phase)), [("H1E2", "HAB"), ("H1E2", "NOV")])
            fmt = root / "B5" / "NOR" / "formatted"
            out = sorted(p.relative_to(fmt).as_posix() for p in fmt.rglob("*.csv"))     # flat: no phase folders
            self.assertEqual(out, ["H1E2_HAB.merged_locs_formatted.csv", "H1E2_NOV.merged_locs_formatted.csv"])

    def test_unknown_code_is_a_clear_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self.make_tree(Path(tmp))
            with self.assertRaises(ValueError) as ctx:
                sd.run("B5", "NOR", root=Path(tmp), code="ZZZZ", verbose=False)
            self.assertIn("ZZZZ", str(ctx.exception))


class B1NorRunTest(unittest.TestCase):
    """B1 NOR after its move out of NOR_NOV_tracking/: geom/ and animal/ with the clip between code and phase, R8O3
    exported as R803. The formatted files carry the code as staging spells it."""

    def test_clip_names_and_the_misspelt_code(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sub, cols in (("geom", {"tl_x": [1.0, 2.0], "tl_y": [3.0, 4.0]}),
                              ("animal", {"nose_x": [5.0, 6.0], "nose_y": [7.0, 8.0]})):
                d = root / "B1" / "NOR" / sub
                d.mkdir(parents=True)
                for i, stem in enumerate(("A1L9_GX010405_NOV", "R803_GX010395_NOV")):
                    name = f"NOR_B1_geom_predict.{i:03d}_{stem}" if sub == "geom" else stem
                    pd.DataFrame({k: [v + i for v in vs] for k, vs in cols.items()}).to_csv(
                        d / f"{name}.analysis_locs.csv", index=False)
            report = sd.run("B1", "NOR", root=root, verbose=False)
            self.assertEqual(sorted(zip(report.code, report.phase)), [("A1L9", "NOV"), ("R8O3", "NOV")])
            fmt = root / "B1" / "NOR" / "formatted"
            self.assertEqual(sorted(p.name for p in fmt.glob("*.csv")),
                             ["A1L9_NOV.merged_locs_formatted.csv", "R8O3_NOV.merged_locs_formatted.csv"])
            report = sd.run("B1", "NOR", root=root, code="R8O3", overwrite=True, verbose=False)
            self.assertEqual(list(report.out), ["R8O3_NOV.merged_locs_formatted.csv"])


class SocpRunTest(unittest.TestCase):
    """SocP is written into formatted/<PHASE>/; --keep-merged writes into merged_debug/, never into the chain's
    merged/."""

    def test_phase_folders_and_kept_merged_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sub, cols in (("geom", {"tl_x": [1.0, 2.0], "tl_y": [3.0, 4.0]}),
                              ("animal", {"nose_x": [5.0, 6.0], "nose_y": [7.0, 8.0]})):
                d = root / "B2" / "SocP" / sub
                d.mkdir(parents=True)
                for i, ph in enumerate(("S1", "S2")):
                    pd.DataFrame({k: [v + 50 * i for v in vs] for k, vs in cols.items()}).to_csv(
                        d / f"SocP-{sub}-charge4.00{i}_E9_SIS_B2_G6K2_{ph}.analysis_locs.csv", index=False)
            sd.run("B2", "SocP", root=root, keep_merged=True, verbose=False)
            test_dir = root / "B2" / "SocP"
            self.assertEqual(sorted(p.relative_to(test_dir).as_posix() for p in test_dir.rglob("*.csv")
                                    if p.parent.name not in ("geom", "animal")),
                             ["formatted/S1/G6K2_S1.merged_locs_formatted.csv",
                              "formatted/S2/G6K2_S2.merged_locs_formatted.csv",
                              "merged_debug/G6K2_S1.merged.csv", "merged_debug/G6K2_S2.merged.csv"])
            self.assertFalse((test_dir / "merged").exists())

    def test_misspelt_and_numbered_exports_take_the_staged_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sub, cols in (("geom", {"tl_x": [1.0, 2.0], "tl_y": [3.0, 4.0]}),
                              ("animal", {"nose_x": [5.0, 6.0], "nose_y": [7.0, 8.0]})):
                d = root / "B2" / "SocP" / sub
                d.mkdir(parents=True)
                for i, stem in enumerate(("SocP-{sub}-charge5.004_E9_SIS_B2_O6O6_S1",
                                          "E9_B2_SocP-{sub}_retrack_G6K2.000_E9_SIS_B2_G6K2_S1_1")):
                    pd.DataFrame({k: [v + i for v in vs] for k, vs in cols.items()}).to_csv(
                        d / f"{stem.format(sub=sub)}.analysis_locs.csv", index=False)
            report = sd.run("B2", "SocP", root=root, verbose=False)
            self.assertEqual(sorted(zip(report.code, report.phase)), [("G6K2", "S1"), ("Q6Q6", "S1")])
            fmt = root / "B2" / "SocP" / "formatted" / "S1"
            self.assertEqual(sorted(p.name for p in fmt.glob("*.csv")),
                             ["G6K2_S1.merged_locs_formatted.csv", "Q6Q6_S1.merged_locs_formatted.csv"])
            for code in ("Q6Q6", "G6K2"):
                report = sd.run("B2", "SocP", root=root, code=code, overwrite=True, verbose=False)
                self.assertEqual(list(report.out), [f"{code}_S1.merged_locs_formatted.csv"])


class FormatMergedTest(unittest.TestCase):
    """format_merged writes what process_one writes, from the merged csv, and renames bodyparts on request."""

    def test_same_file_as_process_one_and_rename(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            geom = pd.DataFrame({"tl_x": [1.0, 2.0], "tl_y": [3.0, 4.0]})
            animal = pd.DataFrame({"center_x": [5.0, 6.0], "center_y": [7.0, 8.0]})
            geom.to_csv(root / "geom.csv", index=False)
            animal.to_csv(root / "animal.csv", index=False)
            pd.concat([geom, animal], axis=1).to_csv(root / "merged.csv", index=False)

            sd.process_one(root / "geom.csv", root / "animal.csv", root / "a.csv")
            info = sd.format_merged(root / "merged.csv", root / "b.csv")
            self.assertEqual(info, {"bodyparts": 2, "cols": 7, "frames": 2})
            self.assertEqual((root / "a.csv").read_bytes(), (root / "b.csv").read_bytes())

            sd.format_merged(root / "merged.csv", root / "c.csv", {"center": "bodycentre"})
            rows = (root / "c.csv").read_bytes().decode("utf-8").splitlines()
            self.assertEqual(rows[1], "bodyparts,tl,tl,tl,bodycentre,bodycentre,bodycentre")
            self.assertEqual(rows[3:], ["0,1.0,3.0,1,5.0,7.0,1", "1,2.0,4.0,1,6.0,8.0,1"])

    def test_ten_significant_digits_as_the_open_field_b1_b2_files(self) -> None:
        # half-up on the number's shortest text: '%.10g' would give 1.23456789 for 1.2345678905 (its binary value
        # is just below the half); trailing zeros dropped; an empty value stays empty
        self.assertEqual(sd.significant([415.03132971989045, 126.51241302490234, 920.0, 1.2345678905, -0.000123456789012,
                                         float("nan"), 0.0], 10),
                         ["415.0313297", "126.512413", "920", "1.234567891", "-0.000123456789", "", "0"])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame({"tl_x": [129.52841186523438, 1.0], "tl_y": [float("nan"), 2.5]}).to_csv(root / "m.csv", index=False)
            sd.format_merged(root / "m.csv", root / "f.csv", digits=10)
            rows = (root / "f.csv").read_bytes().decode("utf-8").split("\r\n")
            self.assertEqual(rows[3:5], ["0,129.5284119,,1", "1,1,2.5,1"])


if __name__ == "__main__":
    unittest.main()
