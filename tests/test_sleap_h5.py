"""Tests for the SLEAP .analysis.h5 -> locs.csv conversion."""

from __future__ import annotations

import ast
import json
import os
import re
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

MODULE_DIR = Path(__file__).resolve().parents[1] / "01_SLEAPcoords"
sys.path.insert(0, str(MODULE_DIR))

import exp9_paths  # noqa: E402
import sleap_h5 as sh  # noqa: E402

GEOM_NODES = ["tl", "tr", "br", "bl", "objL", "objR"]

#: EXP9_CORPUS_TESTS=1 asks for the tests that read the tracking on the share: a missing folder is then a failure.
CORPUS = bool(os.environ.get("EXP9_CORPUS_TESTS"))


def need_tracking(test: unittest.TestCase, folders: list[Path], always: bool = False) -> None:
    """Skip *test* unless the corpus tests were asked for (or *always*) and *folders* exist; asked for, a missing
    folder fails it (the tracking root may not be there yet)."""
    if not CORPUS and not always:
        raise unittest.SkipTest("set EXP9_CORPUS_TESTS=1 with the Exp9 tracking on the share available")
    missing = [str(f) for f in folders if not f.is_dir()]
    if missing and CORPUS:
        test.fail(f"EXP9_CORPUS_TESTS is set but the tracking is not there: {missing}")
    if missing:
        raise unittest.SkipTest(f"Exp9 tracking not available: {missing[0]}")


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
    """The conversion must give the committed 2024 locs files from the 2024 h5 files (reads the tracking on the share
    when B5/NOR is there; skipped when it is not, a failure with EXP9_CORPUS_TESTS=1)."""

    TEST_DIR = exp9_paths.tracking_dir("B5", "NOR")

    def test_b5_h1e2_geom_and_animal(self) -> None:
        need_tracking(self, [self.TEST_DIR], always=True)
        cases = [
            ("geom", self.TEST_DIR / "geom" / "E9_B5_NOR-geom.005_E9_B5_H1E2_NOR.analysis.h5",
             self.TEST_DIR / "geom" / "E9_B5_NOR-geom.005_E9_B5_H1E2_NOR.analysis_locs.csv"),
            ("animal", self.TEST_DIR / "animal" / "E9_B5_NOR-animal.005_E9_B5_H1E2_NOR.analysis.h5",
             self.TEST_DIR / "animal" / "E9_B5_NOR-animal.005_E9_B5_H1E2_NOR.analysis_locs.csv"),
        ]
        for kind, h5, ref in cases:
            with self.subTest(kind=kind):
                got, _ = sh.extract_locs(h5, kind, "NOR")
                want = pd.read_csv(ref)
                self.assertEqual(list(got.columns), list(want.columns))
                self.assertTrue(np.array_equal(got.isna().values, want.isna().values))
                self.assertLess(float(np.nanmax(np.abs(got.values - want.values))), 1e-9)


# --------------------------------------------------------------------------- #
# EPM, OFT and SocP
# --------------------------------------------------------------------------- #

def nb_fill_geom(Y, kind="linear"):
    """fill_missing of EPMBatch_geom.ipynb / SocPBatch_geom.ipynb, verbatim apart from comments."""
    from scipy.interpolate import interp1d
    if len(Y) == 0:
        return Y
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        if np.isnan(y[0]):
            y[0] = np.nanmedian(y)
        if np.isnan(y[-1]):
            y[-1] = np.nanmedian(y)
        non_nans = np.flatnonzero(~np.isnan(y))
        if non_nans.size > 1:
            f = interp1d(non_nans, y[non_nans], kind=kind, fill_value=np.nan, bounds_error=False)
            xq = np.flatnonzero(np.isnan(y))
            y[xq] = f(xq)
        Y[:, i] = y
    return Y.reshape(initial_shape)


def nb_fill_animal(Y, kind="linear"):
    """fill_missing of EPMBatch_animal.ipynb / SocPBatch_animal.ipynb, verbatim apart from comments."""
    from scipy.interpolate import interp1d
    if len(Y) == 0:
        return Y
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        if np.isnan(y).all():
            y[:] = np.nan
        else:
            if np.isnan(y[0]):
                y[0] = y[np.flatnonzero(~np.isnan(y))][0]
            if np.isnan(y[-1]):
                y[-1] = y[np.flatnonzero(~np.isnan(y))][-1]
            x = np.flatnonzero(~np.isnan(y))
            f = interp1d(x, y[x], kind=kind, fill_value=np.nan, bounds_error=False)
            xq = np.flatnonzero(np.isnan(y))
            y[xq] = f(xq)
        Y[:, i] = y
    return Y.reshape(initial_shape)


def nb_table(locations: np.ndarray, names: list[str]) -> pd.DataFrame:
    """The notebooks' table: nodes labelled by index, x then y."""
    df = pd.DataFrame()
    for index, name in enumerate(names):
        part = locations[:, index, :, :]
        df[f"{name}_x"] = part[:, 0].ravel()
        df[f"{name}_y"] = part[:, 1].ravel()
    return df


def gappy_tracks(n_nodes: int, n: int = 40, seed: int = 0) -> np.ndarray:
    """(track, xy, node, frame) with every kind of gap the rules treat differently."""
    rng = np.random.default_rng(seed)
    t = rng.normal(500, 50, (1, 2, n_nodes, n))
    t[0, :, 0, :3] = np.nan                                  # empty first frames
    t[0, :, 1 % n_nodes, -2:] = np.nan                       # empty last frames
    t[0, :, 2 % n_nodes, 10:15] = np.nan                     # an interior gap
    t[0, :, 3 % n_nodes, :] = np.nan                         # never tracked
    t[0, :, 4 % n_nodes, :] = np.nan
    t[0, :, 4 % n_nodes, 20] = 7.0                           # tracked in one frame only
    return t


def notebook_labels(notebook: str, variable: str) -> list[str]:
    """The node names of a batch notebook's index dict (``part_indices`` or ``geom_indices``), in index order."""
    nb = json.loads((MODULE_DIR / notebook).read_text(encoding="utf-8"))
    for cell in nb["cells"]:
        hit = re.search(rf"{variable}\s*=\s*(\{{.*?\}})", "".join(cell["source"]), re.S)
        if hit:
            index = ast.literal_eval(hit.group(1))
            return [index[i] for i in sorted(index)]
    raise AssertionError(f"{notebook} defines no {variable}")


#: The labels each assay's locs files carry, from the notebooks that wrote them. The open field has no such notebook
#: (OFT/OFTBatch_*.ipynb name the corners tl, tr, br, bl, before and after their June 2026 rewrite): its labels are
#: the header of its locs files, <tracking root>/B1-B6/OFT/{geom,animal}/*_locs.csv (B3-B6; B1 and B2 write
#: 'center' for bodycentre, BATCH_LABELS), checked by OpenFieldCorpusLabelsTest.
LABELS = {
    ("NOR", "geom"): notebook_labels("NOR/NORBatch_geom.ipynb", "geom_indices"),
    ("NOR", "animal"): notebook_labels("NOR/NORBatch_animal.ipynb", "part_indices"),
    ("EPM", "geom"): notebook_labels("EPM/EPMBatch_geom.ipynb", "geom_indices"),
    ("EPM", "animal"): notebook_labels("EPM/EPMBatch_animal.ipynb", "part_indices"),
    ("SocP", "geom"): notebook_labels("Social Preference/SocPBatch_geom.ipynb", "geom_indices"),
    ("SocP", "animal"): notebook_labels("Social Preference/SocPBatch_animal.ipynb", "part_indices"),
    ("OFT", "geom"): ["tl", "tr", "bl", "br"],
    ("OFT", "animal"): ["nose", "leftEar", "rightEar", "leftFlank", "rightFlank", "tailbase", "spine1", "leftArm",
                        "rightArm", "bodycentre", "spine2", "tail1", "tail2", "tailtip"],
}


class NotebookRulesTest(unittest.TestCase):
    """The EPM, OFT and SocP notebooks fill gaps exactly as the NOR rules do; only the node names differ."""

    def test_fill_rules_equal_the_notebook_functions(self) -> None:
        import warnings
        y = gappy_tracks(6).T                                # (frame, node, xy, track)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # the notebook takes the nanmedian of a never-tracked column
            self.assertTrue(np.array_equal(sh.fill_missing_geom(y), nb_fill_geom(y.copy()), equal_nan=True))
        self.assertTrue(np.array_equal(sh.fill_missing_animal(y), nb_fill_animal(y.copy()), equal_nan=True))

    def test_extract_gives_the_notebook_table_for_each_assay(self) -> None:
        import warnings
        for assay in ("EPM", "OFT", "SocP"):
            for kind in ("geom", "animal"):
                with self.subTest(assay=assay, kind=kind), tempfile.TemporaryDirectory() as tmp:
                    names = LABELS[(assay, kind)]
                    t = gappy_tracks(len(names), seed=len(names))
                    p = write_h5(Path(tmp) / "x.analysis.h5", t, sh.H5_NODES[(assay, kind)])
                    got, _ = sh.extract_locs(p, kind, assay)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        filler = nb_fill_geom if kind == "geom" else nb_fill_animal
                        want = nb_table(filler(t.T.copy()), names)
                    pd.testing.assert_frame_equal(got, want)
                    self.assertEqual(got.to_csv(index=False), want.to_csv(index=False))

    def test_node_names_are_the_labels_of_the_files(self) -> None:
        for key, want in LABELS.items():
            with self.subTest(assay=key[0], kind=key[1]):
                self.assertEqual(sh.NODES[key], want)


class StoredNodesTest(unittest.TestCase):
    def export(self, tmp: str, assay: str, kind: str, names: list[str] | None = None, tracks: int = 1) -> Path:
        names = names or sh.H5_NODES[(assay, kind)]
        t = np.repeat(gappy_tracks(len(names)), tracks, axis=0)
        return write_h5(Path(tmp) / "x.analysis.h5", t, names)

    def test_another_node_order_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            names = list(sh.H5_NODES[("OFT", "geom")])
            names[2], names[3] = names[3], names[2]          # bottom_right before bottom_left
            with self.assertRaises(ValueError) as ctx:
                sh.extract_locs(self.export(tmp, "OFT", "geom", names), "geom", "OFT")
            self.assertIn("mislabelled", str(ctx.exception))

    def test_more_than_one_track_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError) as ctx:
                sh.extract_locs(self.export(tmp, "SocP", "animal", tracks=2), "animal", "SocP")
            self.assertIn("2 tracks", str(ctx.exception))

    def test_spellings_of_the_stored_names_are_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            names = [n.replace(" ", "_") for n in sh.H5_NODES[("SocP", "animal")]]   # SocP B1: left_ear
            df, _ = sh.extract_locs(self.export(tmp, "SocP", "animal", names), "animal", "SocP")
            self.assertEqual(list(df.columns[:4]), ["nose_x", "nose_y", "leftEar_x", "leftEar_y"])

    def test_batch_labels(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = self.export(tmp, "EPM", "animal")
            self.assertEqual(sh.extract_locs(p, "animal", "EPM", batch="B3")[0].columns[12], "neck_x")
            self.assertEqual(sh.extract_locs(p, "animal", "EPM", batch="B6")[0].columns[12], "spine1_x")
            self.assertEqual(sh.extract_locs(p, "animal", "EPM")[0].columns[28], "headcentre_x")
        with tempfile.TemporaryDirectory() as tmp:
            p = self.export(tmp, "OFT", "animal")
            self.assertEqual(sh.extract_locs(p, "animal", "OFT", batch="B1")[0].columns[18], "center_x")
            self.assertEqual(sh.extract_locs(p, "animal", "OFT", batch="B4")[0].columns[18], "bodycentre_x")

    def test_nor_stored_names_are_not_checked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = write_h5(Path(tmp) / "g.analysis.h5", geom_tracks(), ["a", "b", "c", "d", "e", "f"])
            df, _ = sh.extract_locs(p, "geom", "NOR")
            self.assertEqual(list(df.columns[:2]), ["tl_x", "tl_y"])


class ReproducesCorpusOtherAssaysTest(unittest.TestCase):
    """One recording per assay and label variant, byte for byte (reads the tracking on the share; EXP9_CORPUS_TESTS=1)."""

    TRACKING = exp9_paths.TRACKING_ROOT
    CASES = [("EPM", "B3", "animal/E9_SIS_B3_EPM-animal.000_E9_SIS_B3_B3N2.analysis.h5"),
             ("EPM", "B6", "geom/E9_B6_EPM-geom.000_E9_B6_A2C8_EPM.analysis.h5"),
             ("OFT", "B4", "animal/E9_B4_OFT-animal.000_E9_B4_A0I4_OFT.analysis.h5"),
             ("OFT", "B4", "geom/E9_B4_OFT-geom.000_E9_B4_A0I4_OFT.analysis.h5"),
             ("SocP", "B2", "animal/SocP-animal-charge4.000_E9_SIS_B2_A0S5_HAB.analysis.h5"),
             ("SocP", "B2", "geom/SocP-geom-charge4.000_E9_SIS_B2_A0S5_HAB.analysis.h5")]

    def test_locs_files_are_rebuilt_byte_for_byte(self) -> None:
        import io
        need_tracking(self, sorted({exp9_paths.tracking_dir(b, a, self.TRACKING) for a, b, _ in self.CASES}))
        for assay, batch, rel in self.CASES:
            with self.subTest(assay=assay, batch=batch, file=rel):
                h5 = exp9_paths.tracking_dir(batch, assay, self.TRACKING) / rel
                df, _ = sh.extract_locs(h5, rel.split("/")[0], assay, batch=batch)
                buf = io.BytesIO()
                df.to_csv(buf, index=False, lineterminator="\r\n")
                self.assertEqual(buf.getvalue(), h5.with_name(h5.name[:-3] + "_locs.csv").read_bytes())


class OpenFieldCorpusLabelsTest(unittest.TestCase):
    """The open-field labels are those of the locs files on the share (reads one header line per batch and kind)."""

    TRACKING = ReproducesCorpusOtherAssaysTest.TRACKING
    BATCHES = ("B1", "B2", "B3", "B4", "B5", "B6")

    def test_open_field_locs_headers(self) -> None:
        need_tracking(self, [exp9_paths.tracking_dir(b, "OFT", self.TRACKING) for b in self.BATCHES])
        for batch in self.BATCHES:
            for kind in ("geom", "animal"):
                with self.subTest(batch=batch, kind=kind):
                    first = sorted((exp9_paths.tracking_dir(batch, "OFT", self.TRACKING) / kind).glob("*_locs.csv"))[0]
                    with open(first, encoding="utf-8") as fh:
                        header = fh.readline().strip().split(",")
                    rename = sh.BATCH_LABELS.get(("OFT", batch, kind), {})
                    self.assertEqual(header, [f"{rename.get(n, n)}_{a}" for n in LABELS[("OFT", kind)] for a in "xy"])


if __name__ == "__main__":
    unittest.main()
