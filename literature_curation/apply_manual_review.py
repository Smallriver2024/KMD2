from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pandas as pd

from clean_llm_evidence import (
    clean_text,
    evidence_tier,
    join,
    key_text,
    load_hgnc_map,
    load_kinase_maps,
    normalize_disease,
    normalize_kinase_field,
    normalize_target,
)


HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
BASE = OUT / "llm_phosphorylation_evidence_cleaned_kinase_centered.parquet"
WORKBOOK = OUT / "manual_review_workbook.xlsx"


def is_blank(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return str(value).strip() == ""


def first_nonblank(*values: object) -> str:
    for value in values:
        if not is_blank(value):
            return clean_text(value)
    return ""


def split_manual_names(value: object) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    return [p.strip() for p in re.split(r"\s*(?:;|/|,)\s*", text) if p.strip()]


def map_manual_target(value: object, role: object, aliases: dict[str, list[str]], by_symbol: dict[str, dict], hgnc_map: dict[str, dict]) -> dict:
    manual = clean_text(value)
    if not manual:
        return {}

    kinase_symbols: list[str] = []
    kinase_uniprots: list[str] = []
    non_kinase_symbols: list[str] = []
    non_kinase_uniprots: list[str] = []
    unmatched: list[str] = []

    for part in split_manual_names(manual):
        kin = normalize_kinase_field(part, aliases, by_symbol)
        if kin["symbols"]:
            kinase_symbols.extend(kin["symbols"])
            kinase_uniprots.extend(kin["uniprot_ids"])
            continue
        target = normalize_target(part, role, aliases, by_symbol, hgnc_map)
        if target["non_kinase_symbols"]:
            non_kinase_symbols.extend(target["non_kinase_symbols"])
            non_kinase_uniprots.extend(target["non_kinase_uniprot_ids"])
        else:
            unmatched.append(part)

    kinase_symbols = list(dict.fromkeys(kinase_symbols))
    kinase_uniprots = list(dict.fromkeys(kinase_uniprots))
    non_kinase_symbols = list(dict.fromkeys(non_kinase_symbols))
    non_kinase_uniprots = list(dict.fromkeys(non_kinase_uniprots))
    unmatched = list(dict.fromkeys(unmatched))

    if kinase_symbols and not unmatched:
        status = "manual_mapped_kinase"
    elif non_kinase_symbols and not unmatched:
        status = "manual_mapped_non_kinase"
    elif kinase_symbols or non_kinase_symbols:
        status = "manual_partial_review"
    else:
        status = "manual_unmapped_review"

    return {
        "phospho_target_clean": manual,
        "phospho_target_kinase_symbols": join(kinase_symbols),
        "phospho_target_kinase_uniprot_ids": join(kinase_uniprots),
        "phospho_target_non_kinase_symbols": join(non_kinase_symbols),
        "phospho_target_non_kinase_uniprot_ids": join(non_kinase_uniprots),
        "phospho_target_unmatched": join(unmatched),
        "phospho_target_norm_status": status,
    }


def recompute_center(row: pd.Series) -> tuple[str, str]:
    symbols: list[str] = []
    uniprots: list[str] = []
    if clean_text(row.get("phospho_target_role")) == "kinase":
        symbols += split_manual_names(row.get("phospho_target_kinase_symbols"))
        uniprots += split_manual_names(row.get("phospho_target_kinase_uniprot_ids"))
    symbols += split_manual_names(row.get("upstream_kinase_symbols"))
    uniprots += split_manual_names(row.get("upstream_kinase_uniprot_ids"))
    return join(list(dict.fromkeys(symbols))), join(list(dict.fromkeys(uniprots)))


def apply_entity_level_edits(df: pd.DataFrame, xl: pd.ExcelFile, aliases: dict[str, list[str]], by_symbol: dict[str, dict], hgnc_map: dict[str, dict]) -> dict:
    stats = {}

    if "target_review" in xl.sheet_names:
        target_review = pd.read_excel(xl, "target_review", dtype=str).fillna("")
        target_map = {}
        for _, row in target_review.iterrows():
            raw = clean_text(row.get("phospho_target"))
            # ``phospho_target_clean`` is a prefilled suggestion in the review
            # template, not evidence of a human decision.  Only the explicit
            # manual field may trigger a manual-review label or an overwrite.
            manual = clean_text(row.get("manual_phospho_target_clean"))
            if raw and manual:
                target_map[raw] = manual

        changed = 0
        for raw, manual in target_map.items():
            mask = df["phospho_target"].fillna("").astype(str).map(clean_text).eq(raw)
            if not mask.any():
                continue
            for idx in df.index[mask]:
                mapped = map_manual_target(manual, df.at[idx, "phospho_target_role"], aliases, by_symbol, hgnc_map)
                for col, value in mapped.items():
                    df.at[idx, col] = value
                df.at[idx, "manual_merge_notes"] += f"target_review:{raw}->{manual};"
                df.at[idx, "manual_review_applied"] = "yes"
                changed += 1
        stats["target_review_rows_applied"] = changed

    if "agent_context_review" in xl.sheet_names:
        agent_review = pd.read_excel(xl, "agent_context_review", dtype=str).fillna("")
        agent_map = {}
        for _, row in agent_review.iterrows():
            raw = clean_text(row.get("primary_agent"))
            # Do not treat a prefilled template value as an approved edit.
            clean = clean_text(row.get("manual_primary_agent_clean"))
            if raw and clean:
                agent_map[raw] = clean
        changed = 0
        for raw, clean in agent_map.items():
            mask = df["primary_agent"].fillna("").astype(str).map(clean_text).eq(raw)
            if not mask.any():
                continue
            df.loc[mask, "primary_agent_clean"] = clean
            df.loc[mask, "primary_agent_norm_status"] = "manual_cleaned"
            df.loc[mask, "manual_review_applied"] = "yes"
            df.loc[mask, "manual_merge_notes"] = df.loc[mask, "manual_merge_notes"] + f"agent_review:{raw}->{clean};"
            changed += int(mask.sum())
        stats["agent_review_rows_applied"] = changed

    if "disease_review" in xl.sheet_names:
        disease_review = pd.read_excel(xl, "disease_review", dtype=str).fillna("")
        disease_map = {}
        for _, row in disease_review.iterrows():
            raw = clean_text(row.get("disease"))
            # Do not treat a prefilled template value as an approved edit.
            clean = clean_text(row.get("manual_disease_clean"))
            if raw and clean:
                disease_map[raw] = clean
        changed = 0
        for raw, clean in disease_map.items():
            mask = df["disease"].fillna("").astype(str).map(clean_text).eq(raw)
            if not mask.any():
                continue
            normalized = normalize_disease(clean)
            df.loc[mask, "disease_clean"] = normalized["clean"] or clean
            df.loc[mask, "disease_norm_status"] = "manual_cleaned" if normalized["status"] != "review" else "manual_review"
            df.loc[mask, "disease_review_reason"] = normalized["review_reason"]
            df.loc[mask, "manual_review_applied"] = "yes"
            df.loc[mask, "manual_merge_notes"] = df.loc[mask, "manual_merge_notes"] + f"disease_review:{raw}->{clean};"
            changed += int(mask.sum())
        stats["disease_review_rows_applied"] = changed

    return stats


def apply_row_level_edits(
    df: pd.DataFrame,
    original_df: pd.DataFrame,
    xl: pd.ExcelFile,
    aliases: dict[str, list[str]],
    by_symbol: dict[str, dict],
    hgnc_map: dict[str, dict],
) -> dict:
    if "row_review_top5000" not in xl.sheet_names:
        return {"row_level_rows_seen": 0, "row_level_rows_applied": 0}
    rows = pd.read_excel(xl, "row_review_top5000", dtype=str).fillna("")
    rows_seen = 0
    rows_applied = 0

    for _, row in rows.iterrows():
        if is_blank(row.get("source_row")):
            continue
        try:
            idx = int(float(row["source_row"]))
        except ValueError:
            continue
        if idx not in df.index:
            continue
        rows_seen += 1
        applied = False

        target_clean = clean_text(row.get("phospho_target_clean"))
        original_target_clean = clean_text(original_df.at[idx, "phospho_target_clean"])
        if target_clean and target_clean != original_target_clean:
            mapped = map_manual_target(target_clean, df.at[idx, "phospho_target_role"], aliases, by_symbol, hgnc_map)
            for col, value in mapped.items():
                df.at[idx, col] = value
            df.at[idx, "manual_merge_notes"] += f"row_target_clean->{target_clean};"
            applied = True

        agent_clean = clean_text(row.get("primary_agent_clean"))
        original_agent_clean = clean_text(original_df.at[idx, "primary_agent_clean"])
        if agent_clean and agent_clean != original_agent_clean:
            df.at[idx, "primary_agent_clean"] = agent_clean
            df.at[idx, "primary_agent_norm_status"] = "manual_cleaned"
            df.at[idx, "manual_merge_notes"] += f"row_agent_clean->{agent_clean};"
            applied = True

        disease_clean = clean_text(row.get("disease_clean"))
        original_disease_clean = clean_text(original_df.at[idx, "disease_clean"])
        if disease_clean and disease_clean != original_disease_clean:
            df.at[idx, "disease_clean"] = disease_clean
            df.at[idx, "disease_norm_status"] = "manual_cleaned"
            df.at[idx, "manual_merge_notes"] += f"row_disease_clean->{disease_clean};"
            applied = True

        if applied:
            df.at[idx, "manual_review_applied"] = "yes"
            rows_applied += 1

    return {"row_level_rows_seen": rows_seen, "row_level_rows_applied": rows_applied}


def main() -> None:
    if not BASE.exists():
        raise FileNotFoundError(BASE)
    if not WORKBOOK.exists():
        raise FileNotFoundError(WORKBOOK)

    _, by_symbol, aliases = load_kinase_maps()
    hgnc_map = load_hgnc_map()
    df = pd.read_parquet(BASE).copy()
    original_df = df.copy()
    df["manual_review_applied"] = "no"
    df["manual_merge_notes"] = ""

    xl = pd.ExcelFile(WORKBOOK)
    stats = {}
    stats.update(apply_entity_level_edits(df, xl, aliases, by_symbol, hgnc_map))
    stats.update(apply_row_level_edits(df, original_df, xl, aliases, by_symbol, hgnc_map))

    centers = df.apply(recompute_center, axis=1, result_type="expand")
    df["center_kinase_symbols"] = centers[0]
    df["center_kinase_uniprot_ids"] = centers[1]
    df["evidence_tier"] = df.apply(evidence_tier, axis=1)

    out_parquet = OUT / "llm_phosphorylation_evidence_cleaned_kinase_centered_manual_merged.parquet"
    out_tsv = OUT / "llm_phosphorylation_evidence_cleaned_kinase_centered_manual_merged.tsv"
    report_json = OUT / "manual_merge_report.json"
    report_md = OUT / "manual_merge_report.md"
    remaining_xlsx = OUT / "manual_review_remaining_after_merge.xlsx"

    df.to_parquet(out_parquet, index=False)
    df.to_csv(out_tsv, sep="\t", index=False)

    remaining = df[
        (df["evidence_tier"].astype(str).str.endswith("review"))
        | (df["phospho_target_norm_status"].astype(str).str.contains("review", na=False))
        | (df["primary_agent_norm_status"].astype(str).str.contains("review", na=False))
        | (df["disease_norm_status"].astype(str).str.contains("review", na=False))
    ].copy()

    summary = {
        "records": len(df),
        "manual_review_applied_records": int((df["manual_review_applied"] == "yes").sum()),
        "records_with_center_kinase": int((df["center_kinase_symbols"] != "").sum()),
        "records_with_center_kinase_uniprot": int((df["center_kinase_uniprot_ids"] != "").sum()),
        "records_with_clean_sites": int((df["phosphosites_clean"] != "").sum()),
        "evidence_tier_counts": df["evidence_tier"].value_counts().to_dict(),
        "phospho_target_norm_status": df["phospho_target_norm_status"].value_counts().to_dict(),
        "primary_agent_norm_status": df["primary_agent_norm_status"].value_counts().to_dict(),
        "disease_norm_status": df["disease_norm_status"].value_counts().to_dict(),
        "remaining_review_records": len(remaining),
        **stats,
    }

    target_remaining = (
        remaining.groupby(["phospho_target", "phospho_target_clean", "phospho_target_norm_status"], dropna=False)
        .agg(record_count=("pmid", "size"), pmid_examples=("pmid", lambda x: ";".join(map(str, list(x)[:5]))))
        .reset_index()
        .sort_values("record_count", ascending=False)
        .head(1200)
    )
    agent_remaining = (
        remaining.groupby(["primary_agent", "primary_agent_clean", "primary_agent_norm_status"], dropna=False)
        .agg(record_count=("pmid", "size"), pmid_examples=("pmid", lambda x: ";".join(map(str, list(x)[:5]))))
        .reset_index()
        .sort_values("record_count", ascending=False)
        .head(1200)
    )
    disease_remaining = (
        remaining.groupby(["disease", "disease_clean", "disease_norm_status"], dropna=False)
        .agg(record_count=("pmid", "size"), pmid_examples=("pmid", lambda x: ";".join(map(str, list(x)[:5]))))
        .reset_index()
        .sort_values("record_count", ascending=False)
        .head(1200)
    )
    with pd.ExcelWriter(remaining_xlsx, engine="xlsxwriter") as writer:
        pd.DataFrame([{"metric": k, "value": json.dumps(v, ensure_ascii=False) if isinstance(v, dict) else v} for k, v in summary.items()]).to_excel(writer, "summary", index=False)
        target_remaining.to_excel(writer, "target_remaining", index=False)
        agent_remaining.to_excel(writer, "agent_remaining", index=False)
        disease_remaining.to_excel(writer, "disease_remaining", index=False)
        remaining.head(5000).to_excel(writer, "row_remaining_top5000", index=False)
        for sheet in writer.sheets.values():
            sheet.freeze_panes(1, 0)

    summary["outputs"] = {
        "manual_merged_parquet": str(out_parquet),
        "manual_merged_tsv": str(out_tsv),
        "remaining_review_workbook": str(remaining_xlsx),
    }
    report_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "# Manual Review Merge Report",
        "",
        f"- Records: {summary['records']:,}",
        f"- Manual review applied records: {summary['manual_review_applied_records']:,}",
        f"- Records with center kinase: {summary['records_with_center_kinase']:,}",
        f"- Records with center kinase UniProt: {summary['records_with_center_kinase_uniprot']:,}",
        f"- Records with cleaned phosphosite: {summary['records_with_clean_sites']:,}",
        f"- Remaining review records: {summary['remaining_review_records']:,}",
        "",
        "## Evidence Tiers",
    ]
    for key, value in summary["evidence_tier_counts"].items():
        lines.append(f"- {key}: {value:,}")
    lines += [
        "",
        "## Outputs",
        f"- Manual merged parquet: `{out_parquet}`",
        f"- Manual merged TSV: `{out_tsv}`",
        f"- Remaining review workbook: `{remaining_xlsx}`",
    ]
    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
