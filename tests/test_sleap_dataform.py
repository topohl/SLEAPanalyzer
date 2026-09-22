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


if __name__ == "__main__":
    unittest.main()
