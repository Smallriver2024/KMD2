"""Input/resource contracts for public entrypoints; no scientific computation."""
from __future__ import annotations

import json
import os
from pathlib import Path

from .io import sha256
from .ptmsea import read_flat_yaml


def resolve_path(value, root):
    text = str(value)
    if "${KINASEMD_INPUT_ROOT}" in text:
        base = os.environ.get("KINASEMD_INPUT_ROOT")
        if not base:
            raise ValueError("Missing source root: supply --input-root or KINASEMD_INPUT_ROOT (see docs/PORTABLE_INPUTS.md)")
        text = text.replace("${KINASEMD_INPUT_ROOT}", base)
    path = Path(text)
    return path if path.is_absolute() else Path(root) / path


def validate_execution(dataset, rows, exception):
    if dataset.get("dataset_status") in {"HOLD", "PROVENANCE_REVIEW", "REQUIRES_REANALYSIS"}:
        raise ValueError("Dataset design HOLD / provenance review: production execution is blocked")
    if exception.get("comparison_status", "VALIDATED") != "VALIDATED":
        raise ValueError("Comparison exception requires provenance review; execution is blocked")
    minimum = exception.get("min_independent_units") or rows[0].get("min_independent_units")
    if not str(minimum or "").isdigit() or int(minimum) < 2:
        raise ValueError("Missing explicit statistical unit minimum; design HOLD")
    if any(not r.get("biological_replicate_id") for r in rows):
        raise ValueError("Missing biological unit provenance; design HOLD")
    from .config import sample_groups, validate_biological_units
    available = {r.get(field, "") for r in rows for field in ("control_sample", "treatment_sample")}
    sample_groups(rows, available)
    validate_biological_units(rows, exception.get("technical_aggregation") or "none")
    seeds = {r.get("selected_seed", "") for r in rows}
    if len(seeds) != 1 or not next(iter(seeds)).isdigit() or int(next(iter(seeds))) < 1:
        raise ValueError("Missing consistent positive selected seed; provenance review required")


def validate_gmt(gmt, params, *, fixture=False, gct=None, signatures=None):
    from .evidence import read_gmt
    from .gct import read_gct
    cfg = read_flat_yaml(params)
    if not Path(gmt).is_file():
        raise ValueError(f"Missing GMT resource: {gmt}")
    if not cfg.get("gmt_version") or not cfg.get("gmt_sha256"):
        raise ValueError("GMT version and SHA256 must be explicit in --params")
    if sha256(gmt).lower() != cfg["gmt_sha256"].lower():
        raise ValueError("Incompatible GMT: SHA256 differs from selected --params contract")
    if fixture != (cfg["gmt_version"] == "SYNTHETIC"):
        raise ValueError("Synthetic GMT contract requires a DEMO fixture; real data requires its explicit resource version")
    sets = read_gmt(gmt)
    if not sets:
        raise ValueError("Invalid GMT: no signatures")
    if signatures is not None and set(signatures) - set(sets):
        raise ValueError("Incompatible GMT signature membership: statistics contain unknown signatures")
    if gct is not None:
        _, rows = read_gct(gct)
        names = {r["signature"] for r in rows}
        minimum = int(cfg["min_overlap"])
        if not any(len(names & {v[0] for v in members}) >= minimum for members in sets.values()):
            raise ValueError("Incompatible GMT membership: no signature meets configured min_overlap in this GCT")
    return {"gmt_version": cfg["gmt_version"], "gmt_sha256": sha256(gmt),
            "parameters_sha256": sha256(params)}


def validate_ppi_sources(rows, root):
    active = [r for r in rows if r.get("active", "").lower() == "true"]
    if not active:
        raise ValueError("Missing network resources: no active PPI source; configure licensed local tables or a verified --network")
    for row in active:
        path = resolve_path(row["source_file"], root)
        if not path.is_file():
            raise ValueError(f"Missing network resource: {row['source_name']}: {path}")
        if not row.get("sha256") or sha256(path).lower() != row["sha256"].lower():
            raise ValueError(f"PPI source SHA256 missing or mismatched: {row['source_name']}")
        if not row.get("source_version") or row["source_version"] == "UNRESOLVED":
            raise ValueError(f"PPI source version unresolved: {row['source_name']}")


def validate_release_rows(rows):
    for row in rows:
        if row.get("comparison_status") != "VALIDATED" or any(
            "DIAGNOSTIC_ONLY" in str(row.get(key, ""))
            for key in ("release_eligibility", "data_release_eligibility", "execution_scope", "result_status")
        ):
            raise ValueError("Blocked comparison or diagnostic-only result in release")
