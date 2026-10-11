from __future__ import annotations

from collections import Counter


def assemble_phosphoproteomics_release(chunks: list[list[dict]], allowed_statuses: set[str]) -> list[dict]:
    from .resources import validate_release_rows
    if allowed_statuses != {"VALIDATED"}:
        raise ValueError("Release may allow only VALIDATED comparisons")
    output = []
    for chunk in chunks:
        validate_release_rows(chunk)
        for row in chunk:
            if row.get("comparison_status") not in allowed_statuses:
                raise ValueError(f"Blocked comparison in release: {row.get('comparison_id')}")
            output.append(dict(row))
    return output


def frontend_projection(rows: list[dict], cfg: dict) -> list[dict]:
    """Presentation-only projection and rename; no statistical or biological filtering."""
    from .resources import validate_release_rows
    validate_release_rows(rows)
    columns = cfg["frontend_columns"]
    renames = cfg.get("frontend_rename", {})
    if len(set(renames.get(c, c) for c in columns)) != len(columns):
        raise ValueError("Frontend rename creates duplicate columns")
    output = []
    for row in rows:
        missing = set(columns) - set(row)
        if missing:
            raise ValueError(f"Required frontend columns absent: {sorted(missing)}")
        output.append({renames.get(c, c): row[c] for c in columns})
    if len(output) != len(rows):
        raise AssertionError("Frontend projection changed record count")
    return output


def release_audit(phospho: list[dict], annotated: list[dict], frontend: list[dict]) -> list[dict]:
    if len(phospho) != len(annotated) or len(annotated) != len(frontend):
        raise ValueError("Record universe changed between release stages")
    statuses = Counter(r["network_support_status"] for r in annotated)
    if set(statuses) - {"direct_target", "direct_partner", "none"}:
        raise ValueError("Unexpected network support status")
    return [
        {"metric": "phosphoproteomics_rows", "value": len(phospho)},
        {"metric": "annotated_full_rows", "value": len(annotated)},
        {"metric": "frontend_rows", "value": len(frontend)},
        *({"metric": f"network_{key}", "value": statuses[key]}
          for key in ("direct_target", "direct_partner", "none")),
    ]
