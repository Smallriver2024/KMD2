from __future__ import annotations

import math
import statistics
from collections import defaultdict

from .config import sample_groups, validate_biological_units


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the regularized incomplete beta function."""
    max_iter, eps, tiny = 250, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def regularized_beta(x: float, a: float, b: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    log_bt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    log_bt += a * math.log(x) + b * math.log1p(-x)
    bt = math.exp(log_bt)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def two_sided_t_pvalue(t: float, df: float) -> float:
    if not math.isfinite(t):
        return 0.0
    return min(1.0, max(0.0, regularized_beta(df / (df + t * t), df / 2.0, 0.5)))


def paired_test(control: list[float], treatment: list[float]) -> tuple[float, float, float, float]:
    if len(control) != len(treatment) or len(control) < 2:
        raise ValueError("Paired test requires at least two complete explicit pairs")
    delta = [t - c for c, t in zip(control, treatment)]
    mean_c, mean_t = statistics.mean(control), statistics.mean(treatment)
    sd = statistics.stdev(delta)
    t = math.copysign(math.inf, statistics.mean(delta)) if sd == 0 and statistics.mean(delta) != 0 else (
        0.0 if sd == 0 else statistics.mean(delta) / (sd / math.sqrt(len(delta)))
    )
    return mean_c, mean_t, t, two_sided_t_pvalue(t, len(delta) - 1)


def welch_test(control: list[float], treatment: list[float]) -> tuple[float, float, float, float]:
    if len(control) < 2 or len(treatment) < 2:
        raise ValueError("Welch test requires at least two independent patients per group")
    mean_c, mean_t = statistics.mean(control), statistics.mean(treatment)
    vc = statistics.variance(control) / len(control)
    vt = statistics.variance(treatment) / len(treatment)
    if vc + vt == 0:
        t = math.copysign(math.inf, mean_t - mean_c) if mean_t != mean_c else 0.0
        return mean_c, mean_t, t, 0.0 if mean_t != mean_c else 1.0
    t = (mean_t - mean_c) / math.sqrt(vc + vt)
    df = (vc + vt) ** 2 / (
        vc * vc / (len(control) - 1) + vt * vt / (len(treatment) - 1)
    )
    return mean_c, mean_t, t, two_sided_t_pvalue(t, df)


def bh_adjust(pvalues: list[float | None]) -> list[float | None]:
    valid = sorted(((p, i) for i, p in enumerate(pvalues) if p is not None), key=lambda x: x[0])
    out: list[float | None] = [None] * len(pvalues)
    best = 1.0
    for rank, (p, i) in reversed(list(enumerate(valid, 1))):
        best = min(best, p * len(valid) / rank)
        out[i] = best
    return out


def _by_biological_unit(score: dict, rows: list[dict]) -> tuple[list[float], list[float]]:
    groups: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        for group, sample in (("C", row["control_sample"]), ("T", row["treatment_sample"])):
            if sample and score.get(sample) is not None:
                groups[(group, row["biological_replicate_id"])].append(score[sample])
    controls = {bio: statistics.median(vals) for (group, bio), vals in groups.items() if group == "C"}
    treatments = {bio: statistics.median(vals) for (group, bio), vals in groups.items() if group == "T"}
    if rows[0]["statistical_design"] == "PAIRED":
        common = sorted(set(controls) & set(treatments))
        return [controls[x] for x in common], [treatments[x] for x in common]
    return list(controls.values()), list(treatments.values())


def compute_statistics(
    scores: list[dict], comparisons: list[dict], sample_ids: list[str],
    *, min_units: int = 3, technical_aggregation: str = "none",
    kinase_prefix: str = "KINASE-",
) -> list[dict]:
    """All group identity derives from explicit comparison rows, never score-column order."""
    if min_units < 2:
        raise ValueError("min_units must be at least two")
    control_ids, treatment_ids, pairs = sample_groups(comparisons, set(sample_ids))
    validate_biological_units(comparisons, technical_aggregation)
    if technical_aggregation == "site_median":
        raise ValueError("site_median must be applied before PTM-SEA with an explicit sample map")
    design = comparisons[0]["statistical_design"]
    result = []
    for score in scores:
        if technical_aggregation == "score_median":
            control, treatment = _by_biological_unit(score, comparisons)
        elif design == "PAIRED":
            complete = [(score.get(c), score.get(t)) for c, t in pairs
                        if score.get(c) is not None and score.get(t) is not None]
            control = [c for c, _ in complete]
            treatment = [t for _, t in complete]
        else:
            control = [score[c] for c in control_ids if score.get(c) is not None]
            treatment = [score[t] for t in treatment_ids if score.get(t) is not None]
        row = {
            "PXD": comparisons[0]["PXD"],
            "comparison_id": comparisons[0]["comparison_id"],
            "signature": score["signature"],
            "statistical_design": design,
            "n_control_units": len(control),
            "n_treatment_units": len(treatment),
            "n_pairs": len(control) if design == "PAIRED" else "",
            "mean_control": None,
            "mean_treatment": None,
            "mean_delta_treatment_minus_control": None,
            "t_statistic": None,
            "p_value": None,
            "BH_FDR_all": None,
            "BH_FDR_kinase_like": None,
            "signature_direction": "UNTESTED",
            "reason": "",
        }
        if control and treatment:
            mean_control = statistics.mean(control)
            mean_treatment = statistics.mean(treatment)
            descriptive_delta = mean_treatment - mean_control
            row.update(mean_control=mean_control, mean_treatment=mean_treatment,
                       mean_delta_treatment_minus_control=descriptive_delta,
                       signature_direction="UP" if descriptive_delta > 0 else "DOWN" if descriptive_delta < 0 else "NEUTRAL")
        if len(control) < min_units or len(treatment) < min_units:
            row["reason"] = "insufficient_independent_units"
        else:
            mc, mt, t, p = paired_test(control, treatment) if design == "PAIRED" else welch_test(control, treatment)
            delta = mt - mc
            row.update(mean_control=mc, mean_treatment=mt,
                       mean_delta_treatment_minus_control=delta,
                       t_statistic=t, p_value=p,
                       signature_direction="UP" if delta > 0 else "DOWN" if delta < 0 else "NEUTRAL")
        result.append(row)
    all_q = bh_adjust([r["p_value"] for r in result])
    kinase_ix = [i for i, row in enumerate(result) if row["signature"].startswith(kinase_prefix)]
    kinase_q = bh_adjust([result[i]["p_value"] for i in kinase_ix])
    for i, row in enumerate(result):
        row["BH_FDR_all"] = all_q[i]
    for i, q in zip(kinase_ix, kinase_q):
        result[i]["BH_FDR_kinase_like"] = q
    return result
