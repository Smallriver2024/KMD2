"""Classify and deduplicate GPS 5.0 score-based rewiring predictions.

GPS 5.0 itself must first score the reference and mutant sequence windows.
This module does not run GPS or claim to measure phosphorylation abundance; it
only applies the KinaseMD score/cutoff decision rules to GPS output.
"""

from __future__ import annotations

import pandas as pd


def classify_rewiring(reference_score: float, mutant_score: float, cutoff: float) -> str | None:
    """Return the predicted propensity class defined by GPS scores and cutoff.

    ``Gain`` and ``Loss`` cross the kinase-specific cutoff.  ``Up-regulation``
    and ``Down-regulation`` are abbreviated labels for higher/lower *predicted*
    phosphorylation propensity while both scores remain above that cutoff.
    """
    reference_score = float(reference_score)
    mutant_score = float(mutant_score)
    cutoff = float(cutoff)
    if mutant_score > cutoff > reference_score:
        return "Gain"
    if reference_score > cutoff > mutant_score:
        return "Loss"
    if mutant_score > reference_score > cutoff:
        return "Up-regulation"
    if reference_score > mutant_score > cutoff:
        return "Down-regulation"
    return None


def classify_gps_output(
    data: pd.DataFrame,
    reference_column: str = "SCORE1",
    mutant_column: str = "SCORE2",
    cutoff_column: str = "CUTOFF",
    output_column: str = "INFLUENCE",
) -> pd.DataFrame:
    """Apply the four GPS score rules and retain classifiable events only."""
    required = [reference_column, mutant_column, cutoff_column]
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(f"GPS output is missing columns: {missing}")
    result = data.copy()
    result[output_column] = [
        classify_rewiring(reference, mutant, cutoff)
        for reference, mutant, cutoff in zip(
            result[reference_column], result[mutant_column], result[cutoff_column]
        )
    ]
    return result.dropna(subset=[output_column]).reset_index(drop=True)


def deduplicate_rewiring_events(
    data: pd.DataFrame,
    event_columns: tuple[str, ...] = (
        "KIN_ACC_ID", "SUB_ACC_ID", "SUB_MOD_RSD", "SUB_MOD_MUT", "INFLUENCE"
    ),
) -> pd.DataFrame:
    """Deduplicate by kinase, substrate, p-site, mutation, and predicted effect."""
    missing = [column for column in event_columns if column not in data.columns]
    if missing:
        raise ValueError(f"Cannot deduplicate events; missing columns: {missing}")
    return data.drop_duplicates(list(event_columns), keep="first").reset_index(drop=True)
