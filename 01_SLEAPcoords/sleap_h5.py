"""Convert a SLEAP ``.analysis.h5`` export into the ``*_locs.csv`` that ``sleap_dataform.py`` merges.

This is the logic of the per-assay batch notebooks (``NOR/NORBatch_*.ipynb``, ``EPM/EPMBatch_*.ipynb`` and
``Social Preference/SocPBatch_*.ipynb``) as an importable, tested module with the same output. Two kinds of file exist:

* ``geom``   arena corners, maze points, objects or cups. Gaps are bridged linearly and an empty first or last frame
             takes the column median (the geom notebooks' rule).
* ``animal`` the body parts. Gaps are bridged linearly and an empty first or last frame takes the nearest observed
             value (the animal notebooks' rule).

The two rules are the same in every assay's notebooks; only the node names differ (``NODES``, the notebooks' index
dicts). The open field is the exception: no notebook in the repository wrote its locs files. ``OFT/OFTBatch_*.ipynb``
name the corners tl, tr, br, bl, both before and after their June 2026 rewrite (which also changed the fill rules),
while every open-field locs file has tl, tr, bl, br, the order the exports store (top left, top right, bottom left,
bottom right). The open-field names here are taken from those files, which the same two rules reproduce. Each assay
was added after its existing files were reproduced from their own h5 exports (2026-11-12, geom and animal locs
byte-identical):

* NOR   B5 H1E2 (to 1e-9).
* EPM   one recording in each of B1-B6, and B6 R8I2.
* OFT   B1 X1H0 (its h5 files sit under ``SLEAP/old/export_h5_files``), B2 A0S5 geom, B3-B6 one recording each.
* SocP  one recording in each of B1-B6, and B2 G6K2 S2, B2 W2N9 HAB, B4 I4O9 S2.

A few batches labelled one node differently in their locs files (``BATCH_LABELS``: EPM B2-B4 ``neck`` for spine1,
OFT B1/B2 ``center`` for bodycentre); ``batch=`` reproduces those labels.

One exception: the OFT B2 animal locs files of 2023-03-15 left a whole column empty when its first or last frame was
empty (4 of 20 files: A0S5 spine2, O9D3 nose, S1C1 rightEar and rightArm, Z8S9 tailtip; checked byte for byte). No
notebook in the repository does this and the open-field score uses none of those columns; a new B2 export is filled
by the animal rule, like every other batch.

Nodes are taken by position, as in the notebooks, so for EPM, OFT and SocP the node names stored in the export are
checked against ``H5_NODES`` first, and the export must hold exactly one track: a re-ordered skeleton or a second
instance would otherwise be mislabelled or misaligned without an error.

``static`` names landmarks of a geom file that do not move (an object, a corner) and were labelled by hand in a few
frames only. The notebook rule would interpolate them between the labelled frames and the end-of-video median, which
invents a slow drift; a static landmark is held at its first labelled position for the whole video instead. Labelled
frames that disagree with the first by more than ``STATIC_TOLERANCE_PX`` are reported and not used. Only NOR uses it;
the EPM, OFT and SocP corpora interpolate every geom point.

    python sleap_h5.py geom   path/to/geom.analysis.h5   out_locs.csv --assay NOR --static objR
    python sleap_h5.py animal path/to/animal.analysis.h5 out_locs.csv --assay NOR
    python sleap_h5.py animal path/to/animal.analysis.h5 out_locs.csv --assay EPM --batch B3
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

_ANIMAL_14 = ["nose", "leftEar", "rightEar", "leftFlank", "rightFlank", "tailbase", "spine1", "leftArm",
              "rightArm", "bodycentre", "spine2", "tail1", "tail2", "tailtip"]

#: Output names in h5 node order, per assay and kind. Add an assay here only after reproducing its existing
#: ``*_locs.csv`` files (see the module docstring for the files that were reproduced).
NODES: dict[tuple[str, str], list[str]] = {
    ("NOR", "geom"): ["tl", "tr", "br", "bl", "objL", "objR"],
    ("NOR", "animal"): list(_ANIMAL_14),
    ("EPM", "geom"): ["tl", "tr", "ctr", "rt", "rb", "cbr", "br", "bl", "cbl", "lb", "lt", "ctl"],
    ("EPM", "animal"): _ANIMAL_14 + ["headcentre"],
    ("OFT", "geom"): ["tl", "tr", "bl", "br"],
    ("OFT", "animal"): list(_ANIMAL_14),
    ("SocP", "geom"): ["tl", "tr", "br", "bl", "socl", "socr"],
    ("SocP", "animal"): list(_ANIMAL_14),
}

#: Labels a batch wrote in its locs files instead of the ``NODES`` name (same node, same values).
BATCH_LABELS: dict[tuple[str, str, str], dict[str, str]] = {
    ("EPM", "B2", "animal"): {"spine1": "neck"},
    ("EPM", "B3", "animal"): {"spine1": "neck"},
    ("EPM", "B4", "animal"): {"spine1": "neck"},
    ("OFT", "B1", "animal"): {"bodycentre": "center"},
    ("OFT", "B2", "animal"): {"bodycentre": "center"},
}

_H5_ANIMAL_14 = ["nose", "left ear", "right ear", "left flank", "right flank", "tailbase", "spine 1", "left arm",
                 "right arm", "center", "spine 2", "tail 1", "tail 2", "tail tip"]

#: Node names stored in the exports, in order. Compared after ``node_key`` (case, spaces and underscores ignored:
#: SocP B1 stores ``left_ear``, the other batches ``left ear``). NOR has no entry and is not checked.
H5_NODES: dict[tuple[str, str], list[str]] = {
    ("EPM", "geom"): NODES[("EPM", "geom")],
    ("EPM", "animal"): _H5_ANIMAL_14 + ["head center"],
    ("OFT", "geom"): ["top left", "top right", "bottom left", "bottom right"],
    ("OFT", "animal"): _H5_ANIMAL_14,
    ("SocP", "geom"): NODES[("SocP", "geom")],
    ("SocP", "animal"): _H5_ANIMAL_14,
}


def node_key(name: str) -> str:
    """``'left ear'``, ``'left_ear'`` and ``'Left Ear'`` compare equal."""
    return re.sub(r"[\s_]+", "", name).lower()


def labels(assay: str, kind: str, batch: str | None = None) -> list[str]:
    """Column names (without ``_x``/``_y``) of the locs file for *assay*, *kind* and, if given, *batch*."""
    if (assay, kind) not in NODES:
        raise ValueError(f"no node names for assay {assay!r} kind {kind!r}; known: {sorted(NODES)}")
    rename = BATCH_LABELS.get((assay, batch, kind), {}) if batch else {}
    return [rename.get(n, n) for n in NODES[(assay, kind)]]

#: A later labelled frame this far (px) from the first one is not trusted as the position of a static landmark.
STATIC_TOLERANCE_PX = 1.0


def fill_missing_geom(Y: np.ndarray, kind: str = "linear") -> np.ndarray:
    """Geom rule: interior gaps interpolated, empty first/last frame = column median. All-NaN columns stay NaN."""
    if len(Y) == 0:
        return Y
    shape = Y.shape
    Y = Y.reshape((shape[0], -1)).copy()
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        if np.isnan(y).all():
            continue
        if np.isnan(y[0]):
            y[0] = np.nanmedian(y)
        if np.isnan(y[-1]):
            y[-1] = np.nanmedian(y)
        seen = np.flatnonzero(~np.isnan(y))
        if seen.size > 1:
            f = interp1d(seen, y[seen], kind=kind, fill_value=np.nan, bounds_error=False)
            gap = np.flatnonzero(np.isnan(y))
            y[gap] = f(gap)
        Y[:, i] = y
    return Y.reshape(shape)


def fill_missing_animal(Y: np.ndarray, kind: str = "linear") -> np.ndarray:
    """Animal rule: interior gaps interpolated, empty first/last frame = nearest observed value."""
    if len(Y) == 0:
        return Y
    shape = Y.shape
    Y = Y.reshape((shape[0], -1)).copy()
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        if np.isnan(y).all():
            continue
        seen = np.flatnonzero(~np.isnan(y))
        if np.isnan(y[0]):
            y[0] = y[seen][0]
        if np.isnan(y[-1]):
            y[-1] = y[seen][-1]
        seen = np.flatnonzero(~np.isnan(y))
        f = interp1d(seen, y[seen], kind=kind, fill_value=np.nan, bounds_error=False)
        gap = np.flatnonzero(np.isnan(y))
        y[gap] = f(gap)
        Y[:, i] = y
    return Y.reshape(shape)


FILLERS = {"geom": fill_missing_geom, "animal": fill_missing_animal}


def read_tracks(h5_path: Path) -> tuple[np.ndarray, list[str]]:
    """``tracks`` as (frame, node, xy, track) and the node names stored in the file."""
    with h5py.File(h5_path, "r") as f:
        tracks = f["tracks"][:].T
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in f["node_names"][:]]
    return tracks, names


def hold_static(tracks: np.ndarray, node: int, label: str, tolerance: float = STATIC_TOLERANCE_PX) -> tuple[float, float, list[str]]:
    """Overwrite landmark *node* with its first labelled position in every frame; return (x, y, notes)."""
    xy = tracks[:, node, :, 0]                      # (frame, xy), first track
    seen = np.flatnonzero(np.isfinite(xy).all(axis=1))
    if seen.size == 0:
        raise ValueError(f"static landmark {label!r} has no labelled frame; nothing to hold")
    x, y = float(xy[seen[0], 0]), float(xy[seen[0], 1])
    notes = [f"{label}: held at ({x:.4f}, {y:.4f}) from frame {int(seen[0])} for all {tracks.shape[0]} frames"]
    far = [(int(i), float(np.hypot(xy[i, 0] - x, xy[i, 1] - y))) for i in seen[1:] if np.hypot(xy[i, 0] - x, xy[i, 1] - y) > tolerance]
    if far:
        notes.append(f"{label}: ignored {len(far)} labelled frame(s) more than {tolerance:g} px from the first: "
                     + ", ".join(f"frame {i} ({d:.1f} px)" for i, d in far))
    tracks[:, node, 0, :] = x
    tracks[:, node, 1, :] = y
    return x, y, notes


def check_export(h5_name: str, tracks: np.ndarray, stored: list[str], assay: str, kind: str) -> None:
    """Stored node names (in order) and a single track, for the assays listed in ``H5_NODES``."""
    want = H5_NODES.get((assay, kind))
    if want is None:
        return
    if [node_key(n) for n in stored] != [node_key(n) for n in want]:
        raise ValueError(f"{h5_name}: stored nodes {stored} are not the {assay} {kind} nodes {want} in that order; "
                         "nodes are labelled by position, so this export would be mislabelled")
    if tracks.ndim != 4 or tracks.shape[-1] != 1:
        raise ValueError(f"{h5_name}: {tracks.shape[-1] if tracks.ndim == 4 else '?'} tracks; a single-animal "
                         "export has exactly one (delete the extra instances or tracks in SLEAP and export again)")


def extract_locs(h5_path: Path, kind: str, assay: str = "NOR", static: tuple[str, ...] = (),
                 batch: str | None = None) -> tuple[pd.DataFrame, list[str]]:
    """Return the locs table of *h5_path* (columns ``<node>_x``, ``<node>_y``) and a list of notes.

    *batch* gives the labels that batch used in its locs files (``BATCH_LABELS``); without it the ``NODES`` names.
    """
    names = labels(assay, kind, batch)
    if static and kind != "geom":
        raise ValueError("static landmarks only apply to geom files")
    tracks, stored = read_tracks(Path(h5_path))
    if len(stored) != len(names):
        raise ValueError(f"{Path(h5_path).name}: {len(stored)} nodes {stored}, expected {len(names)} for {assay} {kind}")
    check_export(Path(h5_path).name, tracks, stored, assay, kind)
    notes: list[str] = []
    for label in static:
        if label not in names:
            raise ValueError(f"static landmark {label!r} is not one of {names}")
        notes += hold_static(tracks, names.index(label), label)[2]
    tracks = FILLERS[kind](tracks)
    out = pd.DataFrame()
    for i, name in enumerate(names):
        part = tracks[:, i, :, :]
        out[f"{name}_x"] = part[:, 0].ravel()
        out[f"{name}_y"] = part[:, 1].ravel()
    return out, notes


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("kind", choices=["geom", "animal"])
    ap.add_argument("h5", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--assay", default="NOR")
    ap.add_argument("--batch", default=None, help="write that batch's labels (EPM B2-B4 neck, OFT B1/B2 center)")
    ap.add_argument("--static", nargs="*", default=[], help="geom landmarks to hold at their first labelled position")
    args = ap.parse_args(argv)
    try:
        df, notes = extract_locs(args.h5, args.kind, args.assay, tuple(args.static), args.batch)
    except (ValueError, OSError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    for n in notes:
        print(n)
    print(f"{args.out}  ({len(df)} frames, {df.shape[1]} columns, {int(df.isna().sum().sum())} NaN)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
