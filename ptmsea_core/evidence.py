from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from .config import normalize_symbols, target_aliases
from .io import float_or_none
from .sites import representative_site


def read_gmt(path: str | Path) -> dict[str, list[tuple[str, str]]]:
    result = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 3:
                continue
            members = []
            for item in parts[2:]:
                site, _, original_direction = item.partition(";")
                members.append((site, original_direction))
            result[parts[0]] = members
    return result


def is_kinase_signature(signature: str) -> bool:
    # Frozen PTMsigDB class namespace, never a substring match.
    return str(signature).upper().startswith("KINASE-")


def signature_center(signature: str) -> str:
    if not is_kinase_signature(signature):
        return ""
    for prefix in ("KINASE-PSP_", "KINASE-IKIP_", "KINASE-"):
        if signature.upper().startswith(prefix):
            return signature[len(prefix):].split("/")[-1]
    return ""


def tested_status(stat: dict, minimum: int) -> bool:
    """A valid test at the declared unit minimum; positive delta is insufficient."""
    p = float_or_none(stat.get("p_value"))
    if p is None or not 0 <= p <= 1 or stat.get("reason"):
        return False
    design = stat.get("statistical_design", "PAIRED")
    fields = ("n_pairs",) if design == "PAIRED" else ("n_control_units", "n_treatment_units")
    if design not in {"PAIRED", "INDEPENDENT"}:
        return False
    units = [float_or_none(stat.get(field)) for field in fields]
    return all(n is not None and n.is_integer() and n >= minimum for n in units)


def build_site_evidence(
    statistics_rows: list[dict], site_rows: list[dict],
    gmt: dict[str, list[tuple[str, str]]], metadata: dict,
    *, nominal_p_cutoff: float = 0.05, target_fallback_enabled: bool = False,
    aliases: dict[str, list[str]] | None = None,
) -> list[dict]:
    """One frozen first-source representative per motif; separate target evidence."""
    aliases = target_aliases() if aliases is None else aliases
    minimum = float_or_none(metadata.get("min_independent_units", 3))
    if minimum is None or not minimum.is_integer() or minimum < 2:
        raise ValueError("Minimum statistical unit rule is UNRESOLVED; execution is blocked")
    minimum = int(minimum)
    samples = json.loads(metadata.get("site_quantitative_samples_json") or "[]")
    by_name: dict[str, list[dict]] = defaultdict(list)
    for site in site_rows:
        if samples and not any(float_or_none(site.get(sample)) is not None for sample in samples):
            continue
        by_name[site["Name"]].append(site)
    targets = normalize_symbols(metadata.get("drug_target", ""), aliases, require_approved=True)
    if target_fallback_enabled and not targets:
        raise ValueError("Target fallback requested without an explicitly normalized drug_target")
    stats = {r["signature"]: r for r in statistics_rows if is_kinase_signature(r["signature"])}
    observed_signatures = set(stats)
    # Historical v2 target substrates may lack a scored/tested signature.
    # Their separate rule must never masquerade as directional evidence.
    if target_fallback_enabled:
        for sig in gmt:
            if is_kinase_signature(sig) and set(normalize_symbols(signature_center(sig), aliases)) & set(targets):
                stats.setdefault(sig, {"PXD": metadata["PXD"], "comparison_id": metadata["comparison_id"],
                                       "signature": sig, "signature_direction": "UNTESTED"})
    output = []
    for sig, stat in stats.items():
        center = signature_center(sig)
        delta = float_or_none(stat.get("mean_delta_treatment_minus_control"))
        p = float_or_none(stat.get("p_value"))
        tested = tested_status(stat, minimum)
        fallback = target_fallback_enabled and bool(set(normalize_symbols(center, aliases)) & set(targets))
        sensitive_rules = ["TARGET_SITE_FALLBACK"] if fallback else []
        if tested and delta is not None and delta < 0 and p < nominal_p_cutoff:
            sensitive_rules.append("SIGNATURE_DECREASE")
        rules = [";".join(sensitive_rules)] if sensitive_rules else []
        if tested and delta is not None and delta > 0:
            rules.append("SIGNATURE_INCREASE")
        seen_members = set()
        for member, original_direction in gmt.get(sig, []):
            if member in seen_members or member not in by_name:
                continue
            seen_members.add(member)
            sources = by_name[member]
            site = representative_site(sources)
            provenance = [{key: r.get(key, "") for key in
                           ("site_id", "Name", "accession", "position", "gene_raw", "accession_raw",
                            "position_raw", "sequence_window_raw", "source_file", "source_row")}
                          for r in sorted(sources, key=lambda r: (str(r.get("source_file", "")),
                                                                 str(r.get("source_row", "")), str(r["site_id"])))]
            target_status = "not_applicable"
            ambiguous = len({r["site_id"] for r in sources}) > 1 or any(
                ";" in str(r.get("accession_raw", "")) or ";" in str(r.get("gene_raw", "")) for r in sources)
            if fallback:
                target_status = ("in_selected_step3_with_stats" if tested else
                                 "in_selected_step3_but_missing_stats" if sig in observed_signatures else
                                 "not_in_selected_step3_but_GMT_Step2_overlap")
            for rule in rules:
                output.append({
                    "PXD": stat["PXD"], "comparison_id": stat["comparison_id"],
                    "comparison": stat["comparison_id"],
                    "drug": metadata["drug"], "drug_target": metadata.get("drug_target", ""),
                    "drug_target_approved": ";".join(targets),
                    "center_kinase": center, "regulated_site_id": site["site_id"],
                    "regulated_gene": site.get("gene", ""),
                    "regulated_protein": site.get("protein_name", ""),
                    "protein_accession": site.get("accession", ""),
                    "amino_acid": site.get("residue", ""), "position": site.get("position", ""),
                    "sequence_15aa": site.get("sequence_15aa", member.removesuffix("-p")),
                    "sequence_window_raw": site.get("sequence_window_raw", ""),
                    "accession_raw": site.get("accession_raw", ""), "position_raw": site.get("position_raw", ""),
                    "metadata_complete": site.get("metadata_complete", ""),
                    "site_identity_ambiguous": ambiguous,
                    **{key: metadata.get(key) or site.get(key, "") for key in
                       ("cell_line", "disease", "tissue", "model_context", "PMID", "project_title")},
                    "ptmsea_signature": sig, "signature_direction": stat["signature_direction"],
                    "site_response_direction": site.get("site_response_direction") or "UNKNOWN",
                    "signature_direction_original": original_direction,
                    "mean_delta_treatment_minus_control": delta,
                    "p_value": p, "BH_FDR_all": stat.get("BH_FDR_all"),
                    "BH_FDR_kinase_like": stat.get("BH_FDR_kinase_like"),
                    "statistical_test_status": "TESTED" if tested else "NOT_TESTED",
                    "min_independent_units": minimum, "n_pairs": stat.get("n_pairs", ""),
                    "n_control_units": stat.get("n_control_units", ""),
                    "n_treatment_units": stat.get("n_treatment_units", ""),
                    "target_signature_status": target_status, "evidence_rule": rule,
                    "source_file": site.get("source_file", ""), "source_row": site.get("source_row", ""),
                    **{key: site[key] for key in ("source_feature_id", "source_modified_sequence", "source_site_label",
                                                 "source_identity_reference", "source_identity_reconstruction", "source_sequence_window") if key in site},
                    "representative_policy": "first_source_row", "n_source_rows": len(sources),
                    "site_identity_provenance_json": json.dumps(provenance, sort_keys=True),
                    "selected_seed": metadata.get("selected_seed", ""),
                    "comparison_status": metadata.get("comparison_status", ""),
                })
    return output
