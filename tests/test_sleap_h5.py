"""Tests for the SLEAP .analysis.h5 -> locs.csv conversion."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

MODULE_DIR = Path(__file__).resolve().parents[1] / "01_SLEAPcoords"
sys.path.insert(0, str(MODULE_DIR))

import sleap_h5 as sh  # noqa: E402

GEOM_NODES = ["tl", "tr", "br", "bl", "objL", "objR"]


def write_h5(path: Path, tracks: np.ndarray, nodes: list[str]) -> Path:
    """A minimal SLEAP analysis file: tracks is (track, xy, node, frame)."""
    with h5py.File(path, "w") as f:
        f["tracks"] = tracks
        f["node_names"] = np.array([n.encode() for n in nodes])
    return path


def geom_tracks(n: int = 10) -> np.ndarray:
    t = np.full((1, 2, 6, n), np.nan)
    for node in range(5):                                   # five landmarks tracked in every frame
        t[0, 0, node, :] = 100.0 + node
        t[0, 1, node, :] = 200.0 + node
    return t


class FillTest(unittest.TestCase):
    def test_geom_rule_uses_median_at_the_ends_and_interpolates_inside(self) -> None:
        y = np.array([np.nan, 1.0, np.nan, 3.0, 9.0, np.nan]).reshape(6, 1, 1, 1)
        out = sh.fill_missing_geom(y).ravel()
        self.assertEqual(out[0], 3.0)                       # median of 1, 3, 9
        self.assertEqual(out[5], 3.0)
        self.assertEqual(out[2], 2.0)                       # midway between 1 and 3
        self.assertTrue(np.isnan(y).ravel()[0])             # input not mutated

    def test_animal_rule_uses_nearest_value_at_the_ends(self) -> None:
        y = np.array([np.nan, 1.0, np.nan, 3.0, 9.0, np.nan]).reshape(6, 1, 1, 1)
        out = sh.fill_missing_animal(y).ravel()
        self.assertEqual(out[0], 1.0)
        self.assertEqual(out[5], 9.0)
        self.assertEqual(out[2], 2.0)

    def test_all_nan_column_stays_nan_under_both_rules(self) -> None:
        y = np.full((5, 1, 1, 1), np.nan)
        self.assertTrue(np.isnan(sh.fill_missing_geom(y)).all())
        self.assertTrue(np.isnan(sh.fill_missing_animal(y)).all())


class StaticLandmarkTest(unittest.TestCase):
    def labelled(self, positions: dict[int, tuple[float, float]], n: int = 10) -> np.ndarray:
        t = np.full((n, 6, 2, 1), np.nan)                   # (frame, node, xy, track), as read_tracks returns
        for frame, (x, y) in positions.items():
            t[frame, 5, 0, 0], t[frame, 5, 1, 0] = x, y
        return t

    def test_holds_the_first_labelled_position_for_every_frame(self) -> None:
        t = self.labelled({0: (713.3, 695.9), 1: (713.3, 695.9)})
        x, y, notes = sh.hold_static(t, 5, "objR")
        self.assertEqual((x, y), (713.3, 695.9))
        self.assertTrue((t[:, 5, 0, 0] == 713.3).all() and (t[:, 5, 1, 0] == 695.9).all())
        self.assertEqual(len(notes), 1)

    def test_a_labelled_frame_far_from_the_first_is_reported_and_not_used(self) -> None:
        t = self.labelled({0: (713.3, 695.9), 1: (713.3, 695.9), 2: (709.4, 693.8)})
        _, _, notes = sh.hold_static(t, 5, "objR")
        self.assertEqual(len(notes), 2)
        self.assertIn("frame 2", notes[1])
        self.assertTrue((t[:, 5, 0, 0] == 713.3).all())

    def test_no_labelled_frame_is_an_error(self) -> None:
        with self.assertRaises(ValueError):
            sh.hold_static(self.labelled({}), 5, "objR")


class ExtractTest(unittest.TestCase):
    def test_extract_names_columns_and_holds_the_static_landmark(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            t = geom_tracks()
            t[0, 0, 5, 0], t[0, 1, 5, 0] = 50.0, 60.0       # objR labelled in frame 0 only
            p = write_h5(Path(tmp) / "g.analysis.h5", t, GEOM_NODES)
            df, notes = sh.extract_locs(p, "geom", "NOR", static=("objR",))
            self.assertEqual(list(df.columns), [f"{n}_{a}" for n in GEOM_NODES for a in "xy"])
            self.assertEqual(len(df), 10)
            self.assertTrue((df.objR_x == 50.0).all() and (df.objR_y == 60.0).all())
            self.assertEqual(int(df.isna().sum().sum()), 0)
            self.assertEqual(len(notes), 1)

    def test_without_static_an_untracked_landmark_stays_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = write_h5(Path(tmp) / "g.analysis.h5", geom_tracks(), GEOM_NODES)
            df, _ = sh.extract_locs(p, "geom", "NOR")
            self.assertTrue(df.objR_x.isna().all())

    def test_wrong_node_count_unknown_assay_and_static_on_animal_are_errors(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = write_h5(Path(tmp) / "g.analysis.h5", geom_tracks(), GEOM_NODES)
            with self.assertRaises(ValueError):
                sh.extract_locs(p, "animal", "NOR")           # 6 nodes, 14 expected
            with self.assertRaises(ValueError):
                sh.extract_locs(p, "geom", "EPM")
            with self.assertRaises(ValueError):
                sh.extract_locs(p, "animal", "NOR", static=("nose",))
            with self.assertRaises(ValueError):
                sh.extract_locs(p, "geom", "NOR", static=("notALandmark",))


class Reproduces2024CorpusTest(unittest.TestCase):
    """The conversion must give the committed 2024 locs files from the 2024 h5 files (needs the S: share)."""

    SLEAP = Path(r"S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Raw Data\Behavior\B5\NOR\SLEAP")

    @unittest.skipUnless(SLEAP.exists(), "Exp9 raw data share not available")
    def test_b5_h1e2_geom_and_animal(self) -> None:
        cases = [
            ("geom", self.SLEAP / "geom" / "E9_B5_NOR-geom.005_E9_B5_H1E2_NOR.analysis.h5",
             self.SLEAP / "geom" / "E9_B5_NOR-geom.005_E9_B5_H1E2_NOR.analysis_locs.csv"),
            ("animal", self.SLEAP / "animal" / "E9_B5_NOR-animal.005_E9_B5_H1E2_NOR.analysis.h5",
             self.SLEAP / "animal" / "E9_B5_NOR-animal.005_E9_B5_H1E2_NOR.analysis_locs.csv"),
        ]
        for kind, h5, ref in cases:
            with self.subTest(kind=kind):
                got, _ = sh.extract_locs(h5, kind, "NOR")
                want = pd.read_csv(ref)
                self.assertEqual(list(got.columns), list(want.columns))
                self.assertTrue(np.array_equal(got.isna().values, want.isna().values))
                self.assertLess(float(np.nanmax(np.abs(got.values - want.values))), 1e-9)


if __name__ == "__main__":
    unittest.main()
