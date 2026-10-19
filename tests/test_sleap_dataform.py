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


if __name__ == "__main__":
    unittest.main()
