"""Convert a SLEAP ``.analysis.h5`` export into the ``*_locs.csv`` that ``sleap_dataform.py`` merges.

This is the logic of ``NOR/NORBatch_geom.ipynb`` and ``NOR/NORBatch_animal.ipynb`` as an importable, tested module
with the same output (reproduced to 1e-9 from the 2024 NOR files). Two kinds of file exist:

* ``geom``   arena corners and objects: six static landmarks. Gaps are bridged linearly and an empty first or last
             frame takes the column median (the geom notebook's rule).
* ``animal`` the 14 body parts. Gaps are bridged linearly and an empty first or last frame takes the nearest observed
             value (the animal notebook's rule).

``static`` names landmarks of a geom file that do not move (an object, a corner) and were labelled by hand in a few
frames only. The notebook rule would interpolate them between the labelled frames and the end-of-video median, which
invents a slow drift; a static landmark is held at its first labelled position for the whole video instead. Labelled
frames that disagree with the first by more than ``STATIC_TOLERANCE_PX`` are reported and not used.

    python sleap_h5.py geom   path/to/geom.analysis.h5   out_locs.csv --assay NOR --static objR
    python sleap_h5.py animal path/to/animal.analysis.h5 out_locs.csv --assay NOR
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

#: Output names in h5 node order, per assay and kind. Only NOR is verified against the 2024 corpus; add an assay here
#: only after reproducing its existing ``*_locs.csv`` files.
NODES: dict[tuple[str, str], list[str]] = {
    ("NOR", "geom"): ["tl", "tr", "br", "bl", "objL", "objR"],
    ("NOR", "animal"): ["nose", "leftEar", "rightEar", "leftFlank", "rightFlank", "tailbase", "spine1", "leftArm",
                        "rightArm", "bodycentre", "spine2", "tail1", "tail2", "tailtip"],
}

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


def extract_locs(h5_path: Path, kind: str, assay: str = "NOR", static: tuple[str, ...] = ()) -> tuple[pd.DataFrame, list[str]]:
    """Return the locs table of *h5_path* (columns ``<node>_x``, ``<node>_y``) and a list of notes."""
    if (assay, kind) not in NODES:
        raise ValueError(f"no node names for assay {assay!r} kind {kind!r}; known: {sorted(NODES)}")
    if static and kind != "geom":
        raise ValueError("static landmarks only apply to geom files")
    names = NODES[(assay, kind)]
    tracks, stored = read_tracks(Path(h5_path))
    if len(stored) != len(names):
        raise ValueError(f"{Path(h5_path).name}: {len(stored)} nodes {stored}, expected {len(names)} for {assay} {kind}")
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
    ap.add_argument("--static", nargs="*", default=[], help="geom landmarks to hold at their first labelled position")
    args = ap.parse_args(argv)
    try:
        df, notes = extract_locs(args.h5, args.kind, args.assay, tuple(args.static))
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
