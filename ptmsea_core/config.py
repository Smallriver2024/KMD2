from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path

from .io import read_json, read_tsv

ACTIVE = {"VALIDATED"}
DESIGNS = {"PAIRED", "INDEPENDENT", "DESCRIPTIVE_ONLY", "HOLD"}
STATUSES = {"VALIDATED", "PROVENANCE_REVIEW", "REQUIRES_REANALYSIS", "HOLD"}
COMPARISON_COLUMNS = {
    "PXD", "analysis_subset", "comparison_id", "drug",
    "control_condition", "treatment_condition", "statistical_design",
    "control_sample", "treatment_sample", "pair_id",
    "biological_replicate_id", "technical_replicate_id",
    "experimental_block", "selected_seed", "comparison_status", "notes",
}


def load_comparisons(path: str | Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if not rows:
        raise ValueError("comparisons.tsv is empty")
    missing = COMPARISON_COLUMNS - set(rows[0])
    if missing:
        raise ValueError(f"comparisons.tsv missing columns: {sorted(missing)}")
    exception_path = Path(path).with_name("dataset_exceptions.tsv")
    minimums = {}
    if exception_path.exists():
        for exception in read_tsv(exception_path):
            key = (exception["PXD"], exception["comparison_id"])
            if key in minimums:
                raise ValueError("Duplicate dataset exception rows")
            minimums[key] = exception.get("min_independent_units", "")
    for row in rows:
        if row["statistical_design"] not in DESIGNS:
            raise ValueError(f"Invalid design for {row['comparison_id']}")
        if row["comparison_status"] not in STATUSES:
            raise ValueError(f"Invalid status for {row['comparison_id']}")
        minimum = minimums.get((row["PXD"], row["comparison_id"]))
        if minimum:
            row["min_independent_units"] = minimum
    return rows


def selected_comparison(rows: list[dict], pxd: str, comparison_id: str) -> list[dict]:
    selected = [r for r in rows if r["PXD"] == pxd and r["comparison_id"] == comparison_id]
    if not selected:
        raise ValueError(f"Comparison absent: {pxd}/{comparison_id}")
    status = {r["comparison_status"] for r in selected}
    design = {r["statistical_design"] for r in selected}
    if len(status) != 1 or len(design) != 1:
        raise ValueError("Comparison rows disagree about status or statistical design")
    if status != ACTIVE:
        raise ValueError(f"{pxd}/{comparison_id} is {next(iter(status))}; execution is blocked")
    samples = list(dict.fromkeys(r[field] for field in ("control_sample", "treatment_sample")
                                 for r in selected if r[field]))
    return [{**r, "site_quantitative_samples_json": json.dumps(samples)} for r in selected]


def target_aliases(path: str | Path | None = None) -> dict[str, list[str]]:
    """Load explicit, project-confirmed aliases; preserve curated target order."""
    path = Path(path) if path else Path(__file__).resolve().parents[1] / "config/drug_target_aliases.tsv"
    aliases: dict[str, list[str]] = {}
    for row in read_tsv(path):
        raw, approved = row["raw_target"].strip().upper(), row["approved_gene"].strip().upper()
        if not raw or not approved or not approved[0].isalpha() or any(c in approved for c in "/; ,"):
            raise ValueError("Target aliases require an explicit approved gene symbol")
        values = aliases.setdefault(raw, [])
        if approved not in values:
            values.append(approved)
    return aliases


def normalize_symbols(value: str, aliases: dict[str, list[str]], *, require_approved: bool = False) -> list[str]:
    """Semicolons separate curated symbols; slash families require a full alias."""
    output = []
    approved = {symbol for values in aliases.values() for symbol in values}
    for token in str(value or "").split(";"):
        token = token.strip().upper()
        if not token:
            continue
        symbols = aliases.get(token, [token] if (token in approved or not require_approved)
                              and token[0].isalpha() and all(c.isalnum() or c in "-_." for c in token) else [])
        for symbol in symbols:
            if symbol not in output:
                output.append(symbol)
    return output


def sample_groups(rows: list[dict], available: set[str]) -> tuple[list[str], list[str], list[tuple[str, str]]]:
    design = rows[0]["statistical_design"]
    controls = [r["control_sample"] for r in rows if r["control_sample"]]
    treatments = [r["treatment_sample"] for r in rows if r["treatment_sample"]]
    if not controls or not treatments:
        raise ValueError("Both sample groups must be explicit")
    if len(controls) != len(set(controls)) or len(treatments) != len(set(treatments)):
        raise ValueError("Duplicate sample assignment")
    if set(controls) & set(treatments):
        raise ValueError("A sample belongs to both groups")
    absent = (set(controls) | set(treatments)) - available
    if absent:
        raise ValueError(f"Samples missing from score/site table: {sorted(absent)}")
    pairs: list[tuple[str, str]] = []
    if design == "PAIRED":
        ids = [r["pair_id"] for r in rows]
        if any(not v for v in ids) or len(ids) != len(set(ids)):
            raise ValueError("PAIRED requires distinct nonempty pair_id per row")
        if len(controls) != len(rows) or len(treatments) != len(rows):
            raise ValueError("PAIRED requires one control and treatment per pair_id")
        pairs = [(r["control_sample"], r["treatment_sample"]) for r in rows]
    elif design == "INDEPENDENT":
        if any(r["pair_id"] for r in rows):
            raise ValueError("INDEPENDENT must not contain pair_id")
        if any(r["control_sample"] and r["treatment_sample"] for r in rows):
            raise ValueError("INDEPENDENT rows must name one sample, never an implied pair")
    else:
        raise ValueError(f"Statistical design {design} is not executable")
    return controls, treatments, pairs


def validate_biological_units(rows: list[dict], aggregation: str) -> None:
    by_bio: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        sample = row["control_sample"] or row["treatment_sample"]
        bio = row["biological_replicate_id"]
        technical = row["technical_replicate_id"]
        group = "C" if row["control_sample"] else "T"
        if technical and not bio:
            raise ValueError(f"Technical replicate {sample} lacks biological_replicate_id")
        if bio:
            by_bio[(group, bio)].add(sample)
    if aggregation == "none" and any(len(samples) > 1 for samples in by_bio.values()):
        raise ValueError("Technical replicates would be counted as independent units")
    if aggregation == "score_median" and any(not r["biological_replicate_id"] for r in rows):
        raise ValueError("score_median requires biological_replicate_id for every sample")
    if aggregation not in {"none", "site_median", "score_median"}:
        raise ValueError("Unknown technical replicate aggregation")


def release_config(path: str | Path) -> dict:
    cfg = read_json(path)
    if cfg.get("release_kind") not in {"HISTORICAL_RELEASE", "CORRECTED_WORKING_RELEASE"}:
        raise ValueError("release_kind must be explicit")
    if set(cfg.get("allowed_comparison_statuses", [])) != {"VALIDATED"}:
        raise ValueError("Release may allow only VALIDATED comparisons")
    return cfg
