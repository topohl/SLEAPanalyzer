"""Focused tests for the SLEAP-to-DLCAnalyzer CSV formatter."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


MODULE_DIR = Path(__file__).resolve().parents[1] / "01_SLEAPcoords"
sys.path.insert(0, str(MODULE_DIR))

import sleap_dataform as sd  # noqa: E402


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
    """--code formats one animal without parsing the legacy ``_NOR``-labelled siblings (B5 NOR)."""

    def make_tree(self, root: Path) -> None:
        for sub, cols in (("geom", {"tl_x": [1.0, 2.0], "tl_y": [3.0, 4.0]}),
                          ("animal", {"nose_x": [5.0, 6.0], "nose_y": [7.0, 8.0]})):
            d = root / "B5" / "NOR" / "SLEAP" / sub
            d.mkdir(parents=True)
            for name in ("H1E2_NOV", "H1E2_HAB", "F8O9_NOR", "F8O9_NOV"):     # F8O9_NOR is the legacy label
                pd.DataFrame(cols).to_csv(d / f"E9_B5_NOR-{sub}.000_E9_B5_{name}.analysis_locs.csv", index=False)

    def test_without_code_the_legacy_label_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self.make_tree(Path(tmp))
            with self.assertRaises(ValueError):
                sd.run("B5", "NOR", root=Path(tmp), dry_run=True, verbose=False)

    def test_with_code_only_that_animal_is_formatted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.make_tree(root)
            report = sd.run("B5", "NOR", root=root, code="H1E2", verbose=False)
            self.assertEqual(sorted(zip(report.code, report.phase)), [("H1E2", "HAB"), ("H1E2", "NOV")])
            out = sorted(p.name for p in (root / "B5" / "NOR" / "SLEAP" / "formatted").rglob("*.csv"))
            self.assertEqual(out, ["H1E2_HAB.merged_locs_formatted.csv", "H1E2_NOV.merged_locs_formatted.csv"])

    def test_unknown_code_is_a_clear_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self.make_tree(Path(tmp))
            with self.assertRaises(ValueError) as ctx:
                sd.run("B5", "NOR", root=Path(tmp), code="ZZZZ", verbose=False)
            self.assertIn("ZZZZ", str(ctx.exception))


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
