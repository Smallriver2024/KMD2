"""Minimal Kaplan--Meier grouping and two-group log-rank test utilities.

This is a compact reference implementation reconstructed from the preserved
KinaseMD survival group tables.  It deliberately produces no figures and does
not include expression or clinical source data.
"""

from __future__ import annotations

import math
from collections.abc import Iterable

import pandas as pd


def split_lower_upper_half(
    data: pd.DataFrame, expression_column: str = "Expression"
) -> pd.Series:
    """Assign ordered samples to equally sized lower/upper expression groups.

    Samples with missing expression are left missing.  For an odd number of
    samples, the extra sample is assigned to ``High`` after ascending sorting.
    Row order breaks ties deterministically; callers may provide a stable,
    biologically appropriate order if tied values need special handling.
    """
    expression = pd.to_numeric(data[expression_column], errors="coerce")
    valid_index = expression.dropna().sort_values(kind="mergesort").index
    midpoint = len(valid_index) // 2
    groups = pd.Series(pd.NA, index=data.index, dtype="string")
    groups.loc[valid_index[:midpoint]] = "Low"
    groups.loc[valid_index[midpoint:]] = "High"
    return groups


def status_to_event(
    status: Iterable[object], death_labels: Iterable[str] = ("dead", "deceased", "1")
) -> pd.Series:
    """Convert vital status to 1=death, 0=censored, and missing=<NA>.

    Numeric status exports commonly encode death as ``1`` or ``1.0``. Missing
    values must remain missing so that downstream survival analysis excludes
    them rather than silently treating them as censored observations.
    """
    deaths = {str(label).strip().casefold() for label in death_labels}
    values = pd.Series(status, copy=False)
    normalized = values.astype("string").str.strip().str.casefold()
    missing = values.isna() | normalized.isna() | normalized.eq("")
    numeric = pd.to_numeric(values, errors="coerce")
    is_death = normalized.isin(deaths) | numeric.eq(1)
    events = pd.Series(pd.NA, index=values.index, dtype="Int64")
    events.loc[~missing] = is_death.loc[~missing].astype(int)
    return events


def logrank_test(
    durations: Iterable[object], events: Iterable[object], groups: Iterable[object]
) -> dict[str, float | int]:
    """Perform an unadjusted two-sided log-rank test for ``Low`` versus ``High``.

    The implementation uses the usual observed-minus-expected statistic with
    the Mantel--Haenszel variance.  The p-value is the chi-square survival
    probability with one degree of freedom.
    """
    frame = pd.DataFrame({"duration": durations, "event": events, "group": groups})
    frame["duration"] = pd.to_numeric(frame["duration"], errors="coerce")
    frame["event"] = pd.to_numeric(frame["event"], errors="coerce")
    frame = frame.dropna(subset=["duration", "event", "group"])
    frame = frame[frame["duration"].ge(0) & frame["group"].isin(["Low", "High"])].copy()
    frame["event"] = frame["event"].gt(0).astype(int)

    observed_low = expected_low = variance = 0.0
    event_times = sorted(frame.loc[frame["event"].eq(1), "duration"].unique())
    for time in event_times:
        at_risk = frame[frame["duration"].ge(time)]
        n_low = int(at_risk["group"].eq("Low").sum())
        n_high = int(at_risk["group"].eq("High").sum())
        n_total = n_low + n_high
        deaths = frame[(frame["duration"].eq(time)) & (frame["event"].eq(1))]
        d_low = int(deaths["group"].eq("Low").sum())
        d_total = int(len(deaths))
        if n_total == 0 or d_total == 0:
            continue
        observed_low += d_low
        expected_low += d_total * n_low / n_total
        if n_total > 1:
            variance += (
                d_total
                * (n_total - d_total)
                * n_low
                * n_high
                / (n_total**2 * (n_total - 1))
            )

    chi_square = 0.0 if variance == 0 else (observed_low - expected_low) ** 2 / variance
    return {
        "n": int(len(frame)),
        "low_n": int(frame["group"].eq("Low").sum()),
        "high_n": int(frame["group"].eq("High").sum()),
        "chi_square": float(chi_square),
        "logrank_p": float(math.erfc(math.sqrt(chi_square / 2))),
    }


def analyse_expression_survival(
    data: pd.DataFrame,
    duration_column: str = "Days",
    status_column: str = "Status",
    expression_column: str = "Expression",
) -> dict[str, float | int]:
    """Split a cohort by expression and return the unadjusted log-rank result."""
    groups = split_lower_upper_half(data, expression_column)
    events = status_to_event(data[status_column])
    return logrank_test(data[duration_column], events, groups)
