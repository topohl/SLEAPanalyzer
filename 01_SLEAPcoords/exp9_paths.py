"""Where the Exp9 tracking and the behaviour records are.

TRACKING_ROOT  ``<EXP>/Analysis/Behavior/SLEAP``: the SLEAP tracking chain, one folder per test,
               ``<batch>/<assay>/`` with projects/, predictions/, models/, geom/, animal/, merged/ (merged csv files
               and their ``_superseded_<date>/`` archives), formatted/ (SocP: formatted/<HAB|S1|S2>/), old/ and notes.
               B1/NOR has NOR_NOV_tracking/ instead of geom/ and animal/. Beside the batch folders: qc/ (QC
               material) and old/ (2022-23 training material), which are not batches (batch loops list B1-B6).
RECORDS_ROOT   ``<EXP>/Raw Data/Behavior``: the records of each test, ``<batch>/<assay>/`` (novelLoc.txt for NOR,
               novelLocHAB/S1/S2.txt for SocP, run sheets, BORIS/, TSE/).

EXP9_TRACKING_ROOT and EXP9_RECORDS_ROOT point them elsewhere.
"""

from __future__ import annotations

import os
from pathlib import Path

EXP9 = Path(r"S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress")
TRACKING_ROOT = Path(os.environ.get("EXP9_TRACKING_ROOT") or EXP9 / "Analysis" / "Behavior" / "SLEAP")
RECORDS_ROOT = Path(os.environ.get("EXP9_RECORDS_ROOT") or EXP9 / "Raw Data" / "Behavior")


def tracking_dir(batch: str, assay: str, root: Path | None = None) -> Path:
    """The tracking of one test: ``<tracking root>/<batch>/<assay>``."""
    return Path(TRACKING_ROOT if root is None else root) / batch / assay


def records_dir(batch: str, assay: str, root: Path | None = None) -> Path:
    """The records of one test: ``<records root>/<batch>/<assay>``."""
    return Path(RECORDS_ROOT if root is None else root) / batch / assay
