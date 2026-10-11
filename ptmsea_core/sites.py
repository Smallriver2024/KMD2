from __future__ import annotations

import math
import re
import statistics
from collections import defaultdict

from .io import float_or_none, nonmissing

AA = set("ACDEFGHIKLMNPQRSTVWY")
VALID_NAME = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY_]{15}-p$")


def ptmsea_name(sequence_window: str) -> str:
    # MaxQuant's '_' padding fixes the site at the original window center.
    # Removing it before cropping shifts terminal sites to a different residue.
    # Semicolon windows correspond to semicolon accessions: use the first of both.
    seq = str(sequence_window or "").split(";")[0].replace(" ", "").replace("*", "").upper()
    if seq.endswith("-P"):
        seq = seq[:-2]
    if len(seq) > 15:
        center = len(seq) // 2
        seq = seq[center - 7:center + 8]
    candidate = seq + "-p"
    valid_padding = re.fullmatch(r"_*[ACDEFGHIKLMNPQRSTVWY]+_*", seq)
    return candidate if VALID_NAME.fullmatch(candidate) and valid_padding and seq[7] in "STY" else ""


def _single(value: str) -> str:
    token = re.split(r"[;,]", str(value or ""))[0].strip()
    return token if nonmissing(token) else ""


def representative_site(rows: list[dict]) -> dict:
    """Frozen Step2 first-source rule, ordered by explicit input/row provenance."""
    def order(row):
        index = float_or_none(row.get("source_input_index"))
        source_row = float_or_none(row.get("source_row"))
        return (index if index is not None else 0, str(row.get("source_file", "")),
                source_row if source_row is not None else math.inf,
                str(row.get("source_row", "")), str(row.get("site_id", "")))
    return min(rows, key=order)


def prepare_sites(rows: list[dict], samples: list[str], dataset: dict) -> tuple[list[dict], list[dict]]:
    """Preserve documented MaxQuant class-I filter and report every rejection."""
    clean, rejected = [], []
    loc_min = float(dataset.get("localization_min") or 0.75)
    only_policy = dataset.get("only_identified_policy", "UNRESOLVED")
    quant_scale = dataset["quant_scale"]
    if quant_scale not in {"linear_intensity", "linear_ratio", "log2"}:
        raise ValueError(f"Unknown quantitative scale: {quant_scale}")
    for row in rows:
        reasons = []
        if row.get("adapter_rejection_reason"):
            reasons.append(row["adapter_rejection_reason"])
        if row.get("reverse", "").strip() == "+":
            reasons.append("reverse")
        if row.get("contaminant", "").strip() == "+":
            reasons.append("contaminant")
        loc = float_or_none(row.get("localization_probability"))
        if row["modality"].startswith("maxquant") or dataset.get("localization_min"):
            if loc is None or loc < loc_min:
                reasons.append("localization_below_threshold_or_missing")
        if row.get("only_identified_by_site", "").strip() == "+":
            if only_policy == "exclude":
                reasons.append("only_identified_by_site")
            elif only_policy != "retain":
                raise ValueError("Only-identified-by-site policy is UNRESOLVED for this dataset")
        gene = _single(row.get("gene", ""))
        accession = _single(row.get("accession", ""))
        if accession.startswith(("sp|", "tr|")):
            accession = accession.split("|")[1]
        residue = _single(row.get("residue", "")).upper()
        pos_match = re.search(r"\d+", _single(row.get("position", "")))
        position = pos_match.group() if pos_match else ""
        name = ptmsea_name(row.get("sequence_window", ""))
        identity_complete = bool(gene and accession and residue in {"S", "T", "Y"} and position)
        if not name:
            reasons.append("invalid_15aa_sequence_window")
        values = {}
        for sample in samples:
            val = float_or_none(row.get(sample))
            if val is not None and quant_scale != "log2":
                val = math.log2(val) if val > 0 else None
            values[sample] = val
        if all(v is None for v in values.values()):
            reasons.append("all_samples_missing_or_nonpositive")
        if reasons:
            rejected.append({"source_file": row.get("source_file", ""), "source_row": row.get("source_row", ""),
                             "reason": ";".join(reasons)})
            continue
        site_id = f"{gene}_{residue}{position}|{accession}" if identity_complete else f"UNMAPPED|{name}"
        clean.append({
            "PXD": row["PXD"], "analysis_subset": row["analysis_subset"],
            "site_id": site_id, "Name": name, "Description": site_id,
            "gene": gene, "accession": accession, "residue": residue, "position": position,
            "protein_name": row.get("protein_name", ""), "sequence_15aa": name[:-2],
            "metadata_complete": identity_complete,
            "site_response_direction": row.get("site_response_direction") or "UNKNOWN",
            **{key: row.get(key, dataset.get(key, "")) for key in
               ("cell_line", "disease", "tissue", "model_context", "PMID", "project_title")},
            "accession_raw": row.get("accession", ""), "position_raw": row.get("position", ""),
            "gene_raw": row.get("gene", ""),
            "sequence_window_raw": row.get("source_sequence_window", row.get("sequence_window", "")),
            "source_file": row["source_file"], "source_row": row["source_row"],
            "source_input_index": row.get("source_input_index", ""),
            **{key: row[key] for key in ("source_feature_id", "source_modified_sequence", "source_site_label",
                                       "source_identity_reference", "source_identity_reconstruction", "source_sequence_window") if key in row},
            **values,
        })
    return clean, rejected


def build_site_matrix(rows: list[dict], samples: list[str]) -> list[dict]:
    """Median of per-site log2 abundances, matching the frozen Step2 rule."""
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if not any(float_or_none(row.get(sample)) is not None for sample in samples):
            continue
        grouped[row["Name"]].append(row)
    out = []
    for name, source in grouped.items():
        representative = representative_site(source)
        record = {
            "Name": name,
            "Description": representative["site_id"],
            "n_source_rows": len(source),
            "site_ids": ";".join(dict.fromkeys(r["site_id"] for r in source)),
            "genes": ";".join(dict.fromkeys(r["gene"] for r in source)),
            "representative_source_file": representative.get("source_file", ""),
            "representative_source_row": representative.get("source_row", ""),
        }
        for sample in samples:
            valid = [r[sample] for r in source if r.get(sample) is not None]
            record[sample] = statistics.median(valid) if valid else None
        if any(record[s] is not None for s in samples):
            out.append(record)
    return out
