from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from clean_llm_evidence import clean_text, join, key_text, load_hgnc_map, load_kinase_maps, normalize_kinase_field


HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
INPUT = OUT / "llm_phosphorylation_evidence_cleaned_kinase_centered_manual_merged.parquet"


TARGET_RULES = {
    # Ambiguous protein families/classes that are useful but should not be forced to one gene.
    "HISTONE H3": {"display": "Histone H3", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "HISTONE 3": {"display": "Histone H3", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "MYOSIN LIGHT CHAIN": {"display": "Myosin light chain", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "MYOSIN REGULATORY LIGHT CHAIN": {"display": "Myosin regulatory light chain", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "MYOSIN LIGHT-CHAIN": {"display": "Myosin light chain", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "I KAPPA B": {"display": "IkappaB family", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "IKAPPAB": {"display": "IkappaB family", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "SMAD": {"display": "SMAD family", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "SMADS": {"display": "SMAD family", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "MAPKS": {"display": "MAPK family", "type": "kinase_family", "status": "final_mapped_family", "keep": True},
    "FOXO": {"display": "FOXO family", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "SR PROTEINS": {"display": "SR proteins", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "PKA SUBSTRATES": {"display": "PKA substrates", "type": "substrate_set", "status": "final_mapped_family", "keep": True},
    "TRANSCRIPTION FACTORS": {"display": "Transcription factors", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    "MYOFILAMENT PROTEINS": {"display": "Myofilament proteins", "type": "protein_family", "status": "final_mapped_family", "keep": True},
    # Viral and pathogen proteins. Keep but label as non-human.
    "HCV NS5A": {"display": "HCV NS5A", "type": "viral_protein", "status": "final_mapped_viral", "keep": True},
    "NS5A": {"display": "HCV NS5A", "type": "viral_protein", "status": "final_mapped_viral", "keep": True},
    "HIV-1 TAT": {"display": "HIV-1 Tat", "type": "viral_protein", "status": "final_mapped_viral", "keep": True},
    "IE62": {"display": "VZV IE62", "type": "viral_protein", "status": "final_mapped_viral", "keep": True},
    "E7": {"display": "HPV E7", "type": "viral_protein", "status": "final_mapped_viral", "keep": True},
    "CAGA": {"display": "CagA", "type": "bacterial_protein", "status": "final_mapped_pathogen", "keep": True},
    # Non-protein or obvious extraction errors; exclude from final usable evidence.
    "SPHINGOSINE": {"display": "sphingosine", "type": "metabolite", "status": "final_exclude_non_protein", "keep": False},
    "GANCICLOVIR": {"display": "ganciclovir", "type": "drug", "status": "final_exclude_drug_as_target", "keep": False},
    "GEMCITABINE": {"display": "gemcitabine", "type": "drug", "status": "final_exclude_drug_as_target", "keep": False},
    "FTY720": {"display": "FTY720", "type": "drug", "status": "final_exclude_drug_as_target", "keep": False},
    "CHOLINE": {"display": "choline", "type": "metabolite", "status": "final_exclude_non_protein", "keep": False},
    "CERAMIDE": {"display": "ceramide", "type": "metabolite", "status": "final_exclude_non_protein", "keep": False},
    "PROTEIN": {"display": "protein", "type": "generic", "status": "final_exclude_too_generic", "keep": False},
    "TYROSINE-PHOSPHORYLATED PROTEINS": {"display": "tyrosine-phosphorylated proteins", "type": "generic", "status": "final_exclude_too_generic", "keep": False},
    "PROTEIN TYROSINE": {"display": "protein tyrosine", "type": "generic", "status": "final_exclude_too_generic", "keep": False},
    "PI3K PATHWAY": {"display": "PI3K pathway", "type": "pathway", "status": "final_exclude_pathway", "keep": False},
    "MAPK PATHWAY": {"display": "MAPK pathway", "type": "pathway", "status": "final_exclude_pathway", "keep": False},
}


GENE_RULES = {
    "RNA POLYMERASE II": "POLR2A",
    "RNA POLYMERASE-II": "POLR2A",
    "RNA POL II": "POLR2A",
    "RNA POLYMERASE II C-TERMINAL DOMAIN": "POLR2A",
    "RNA POLYMERASE II CTD": "POLR2A",
    "RNAPII": "POLR2A",
    "ALPHA-SYNUCLEIN": "SNCA",
    "TAU PROTEIN": "MAPT",
    "TAU": "MAPT",
    "VE-CADHERIN": "CDH5",
    "E-CADHERIN": "CDH1",
    "COFILIN": "CFL1",
    "BETA-CATENIN": "CTNNB1",
    "SS-CATENIN": "CTNNB1",
    "Β-CATENIN": "CTNNB1",
    "ß-CATENIN": "CTNNB1",
    "P66SHC": "SHC1",
    "ESTROGEN RECEPTOR ALPHA": "ESR1",
    "ESTROGEN RECEPTOR-ALPHA": "ESR1",
    "ESTROGEN RECEPTOR": "ESR1",
    "GLUCOCORTICOID RECEPTOR": "NR3C1",
    "CARDIAC TROPONIN I": "TNNI3",
    "TROPONIN I": "TNNI3",
    "EIF2ALPHA": "EIF2S1",
    "EUKARYOTIC INITIATION FACTOR 2ALPHA": "EIF2S1",
    "EUKARYOTIC TRANSLATION INITIATION FACTOR 2ALPHA": "EIF2S1",
    "EIF4E": "EIF4E",
    "EUKARYOTIC INITIATION FACTOR 4E": "EIF4E",
    "PLCgamma2".upper(): "PLCG2",
    "PLCgamma".upper(): "PLCG1",
    "RIBOSOMAL S6 PROTEIN": "RPS6",
    "RIBOSOMAL S6": "RPS6",
    "S6 RIBOSOMAL PROTEIN": "RPS6",
    "ENDOTHELIAL NO SYNTHASE": "NOS3",
    "ENDOTHELIAL NITRIC OXIDE SYNTHASE": "NOS3",
    "ENOS": "NOS3",
    "ACETYL COA CARBOXYLASE": "ACACA;ACACB",
    "ACETYL-COA CARBOXYLASE": "ACACA;ACACB",
    "GLYCOGEN SYNTHASE": "GYS1;GYS2",
    "PYRUVATE DEHYDROGENASE": "PDHA1",
    "PDHE1ALPHA": "PDHA1",
    "DYNAMIN-RELATED PROTEIN 1": "DNM1L",
    "DRP1": "DNM1L",
    "HEAT SHOCK PROTEIN 27": "HSPB1",
    "HSP27": "HSPB1",
    "STATHMIN": "STMN1",
    "OP18": "STMN1",
    "OP18/STATHMIN": "STMN1",
    "NEPHRIN": "NPHS1",
    "C/EBPBETA": "CEBPB",
    "C/EBPALPHA": "CEBPA",
    "CONNEXIN 43": "GJA1",
    "GJA1": "GJA1",
    "LAMIN A/C": "LMNA",
    "CYTOSOLIC PHOSPHOLIPASE A2": "PLA2G4A",
    "CPLA2": "PLA2G4A",
    "RPA32": "RPA2",
    "MYOSIN-BINDING SUBUNIT": "PPP1R12A",
    "MYOSIN PHOSPHATASE TARGET SUBUNIT 1": "PPP1R12A",
    "MBS": "PPP1R12A",
    "AMYLOID PRECURSOR PROTEIN": "APP",
    "APP": "APP",
    "P120-CATENIN": "CTNND1",
    "CALDESMON": "CALD1",
    "NUCLEOPHOSMIN": "NPM1",
    "HTERT": "TERT",
    "BAND 3": "SLC4A1",
    "NFKAPPAB": "NFKB1",
    "NUCLEAR FACTOR-KAPPAB": "NFKB1",
    "IKAPPABALPHA": "NFKBIA",
    "I KAPPA BALPHA": "NFKBIA",
    "IKBΑ": "NFKBIA",
    "H2A HISTONE FAMILY MEMBER X": "H2AX",
    "HISTONE H2AX": "H2AX",
    "HISTONE H2A.X": "H2AX",
    "H2A.X": "H2AX",
    "GAMMA-H2AX": "H2AX",
    "HISTONE 2AX": "H2AX",
    "SMAD2/3": "SMAD2;SMAD3",
    "SMAD 2/3": "SMAD2;SMAD3",
    "SMAD1/5": "SMAD1;SMAD5",
    "SMAD1/5/8": "SMAD1;SMAD5;SMAD9",
    "SMAD2/SMAD3": "SMAD2;SMAD3",
    "MEK1/2": "MAP2K1;MAP2K2",
    "EXTRACELLULAR SIGNAL-REGULATED KINASE 1/2": "MAPK3;MAPK1",
    "STAT1/2": "STAT1;STAT2",
    "SIGNAL TRANSDUCER AND ACTIVATOR OF TRANSCRIPTION 5": "STAT5A;STAT5B",
}


def split_symbols(value: str) -> list[str]:
    return [x for x in re.split(r"\s*;\s*", value) if x]


def get_hgnc(symbol: str, hgnc: dict[str, dict]) -> tuple[str, str]:
    rec = hgnc.get(key_text(symbol))
    if rec:
        return rec["symbol"], rec.get("uniprot_id", "")
    return symbol, ""


def apply_final_target(row: pd.Series, hgnc: dict[str, dict], kinase_aliases: dict[str, list[str]], kinase_by_symbol: dict[str, dict]) -> dict:
    raw = clean_text(row.get("phospho_target_clean")) or clean_text(row.get("phospho_target"))
    k = key_text(raw)
    compact = re.sub(r"[^A-Z0-9]", "", k)

    # Already mapped cleanly.
    status = str(row.get("phospho_target_norm_status", ""))
    if "mapped" in status and "review" not in status and "unmapped" not in status:
        kinase_symbols = split_symbols(clean_text(row.get("phospho_target_kinase_symbols")))
        non_kinase_symbols = split_symbols(clean_text(row.get("phospho_target_non_kinase_symbols")))
        kinase_uniprots = split_symbols(clean_text(row.get("phospho_target_kinase_uniprot_ids")))
        non_kinase_uniprots = split_symbols(clean_text(row.get("phospho_target_non_kinase_uniprot_ids")))
        symbols = join(list(dict.fromkeys([*kinase_symbols, *non_kinase_symbols])))
        uniprot = join(list(dict.fromkeys([*kinase_uniprots, *non_kinase_uniprots])))
        entity_type = (
            "mixed_kinase_and_non_kinase"
            if kinase_symbols and non_kinase_symbols
            else "kinase" if kinase_symbols else "gene"
        )
        return {
            "final_target_name": clean_text(row.get("phospho_target_clean")) or raw,
            "final_target_symbols": symbols,
            "final_target_uniprot_ids": uniprot,
            "final_target_entity_type": entity_type,
            "final_target_status": status,
            "final_target_keep": True,
        }

    rule = TARGET_RULES.get(k) or TARGET_RULES.get(compact)
    if rule:
        return {
            "final_target_name": rule["display"],
            "final_target_symbols": "",
            "final_target_uniprot_ids": "",
            "final_target_entity_type": rule["type"],
            "final_target_status": rule["status"],
            "final_target_keep": rule["keep"],
        }

    gene_rule = GENE_RULES.get(k) or GENE_RULES.get(compact)
    if gene_rule:
        final_symbols = []
        final_uniprots = []
        for symbol in split_symbols(gene_rule):
            kin = normalize_kinase_field(symbol, kinase_aliases, kinase_by_symbol)
            if kin["symbols"]:
                final_symbols += kin["symbols"]
                final_uniprots += kin["uniprot_ids"]
            else:
                h_symbol, h_uniprot = get_hgnc(symbol, hgnc)
                final_symbols.append(h_symbol)
                if h_uniprot:
                    final_uniprots.append(h_uniprot)
        entity_type = "kinase" if any(s in kinase_by_symbol for s in final_symbols) else "gene"
        return {
            "final_target_name": ";".join(final_symbols),
            "final_target_symbols": join(list(dict.fromkeys(final_symbols))),
            "final_target_uniprot_ids": join(list(dict.fromkeys(final_uniprots))),
            "final_target_entity_type": entity_type,
            "final_target_status": "final_mapped_rule",
            "final_target_keep": True,
        }

    # Last chance: HGNC exact/alias match on current cleaned text.
    rec = hgnc.get(k) or hgnc.get(compact)
    if rec:
        return {
            "final_target_name": rec["symbol"],
            "final_target_symbols": rec["symbol"],
            "final_target_uniprot_ids": rec.get("uniprot_id", ""),
            "final_target_entity_type": "gene",
            "final_target_status": "final_mapped_hgnc",
            "final_target_keep": True,
        }

    if not raw:
        return {
            "final_target_name": "",
            "final_target_symbols": "",
            "final_target_uniprot_ids": "",
            "final_target_entity_type": "blank",
            "final_target_status": "final_exclude_blank_target",
            "final_target_keep": False,
        }

    # Keep remaining named targets only if evidence has a center kinase and target is not generic.
    too_dirty = len(raw) > 90 or re.search(r"\b(pathway|protein phosphorylation|substrates?|activity|phosphorylation)\b", raw, flags=re.I)
    return {
        "final_target_name": raw,
        "final_target_symbols": "",
        "final_target_uniprot_ids": "",
        "final_target_entity_type": "named_target_unmapped",
        "final_target_status": "final_keep_named_unmapped" if not too_dirty else "final_exclude_dirty_target",
        "final_target_keep": not too_dirty,
    }


def classify_final(row: pd.Series) -> tuple[bool, str]:
    if str(row.get("final_target_keep")) != "True" and row.get("final_target_keep") is not True:
        return False, "exclude_bad_or_blank_target"
    if not clean_text(row.get("center_kinase_symbols")):
        return False, "exclude_no_center_kinase"
    if not clean_text(row.get("center_kinase_uniprot_ids")):
        return False, "exclude_center_kinase_no_uniprot"
    if not clean_text(row.get("disease_clean")):
        return False, "exclude_no_disease"
    if not clean_text(row.get("evidence_text")):
        return False, "exclude_no_evidence_text"
    if int(row.get("confidence", 0)) < 8 and str(row.get("needs_manual_review")) == "yes":
        return False, "exclude_low_confidence_review"
    return True, "usable"


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    _, kinase_by_symbol, kinase_aliases = load_kinase_maps()
    hgnc = load_hgnc_map()
    df = pd.read_parquet(INPUT).copy()

    target_df = df.apply(lambda row: pd.Series(apply_final_target(row, hgnc, kinase_aliases, kinase_by_symbol)), axis=1)
    df = pd.concat([df, target_df], axis=1)
    final = df.apply(lambda row: pd.Series(classify_final(row), index=["final_usable", "final_exclusion_reason"]), axis=1)
    df = pd.concat([df, final], axis=1)

    usable = df[df["final_usable"]].copy()
    full_out = OUT / "llm_phosphorylation_evidence_final_with_flags.parquet"
    usable_out = OUT / "llm_phosphorylation_evidence_final_usable.parquet"
    usable_tsv = OUT / "llm_phosphorylation_evidence_final_usable.tsv"
    review_out = OUT / "final_excluded_or_review_targets.xlsx"
    report_md = OUT / "final_cleaning_report.md"
    report_json = OUT / "final_cleaning_report.json"

    df.to_parquet(full_out, index=False)
    usable.to_parquet(usable_out, index=False)
    usable.to_csv(usable_tsv, sep="\t", index=False)

    excluded_summary = (
        df[~df["final_usable"]]
        .groupby(["final_exclusion_reason", "final_target_status", "final_target_name"], dropna=False)
        .agg(record_count=("pmid", "size"), pmid_examples=("pmid", lambda x: ";".join(map(str, list(x)[:5]))))
        .reset_index()
        .sort_values("record_count", ascending=False)
        .head(2000)
    )
    unmapped_keep = (
        usable[usable["final_target_entity_type"].eq("named_target_unmapped")]
        .groupby(["final_target_name", "phospho_target"], dropna=False)
        .agg(record_count=("pmid", "size"), pmid_examples=("pmid", lambda x: ";".join(map(str, list(x)[:5]))))
        .reset_index()
        .sort_values("record_count", ascending=False)
        .head(1000)
    )
    with pd.ExcelWriter(review_out, engine="xlsxwriter") as writer:
        excluded_summary.to_excel(writer, "excluded_summary", index=False)
        unmapped_keep.to_excel(writer, "kept_named_unmapped", index=False)
        for sheet in writer.sheets.values():
            sheet.freeze_panes(1, 0)

    report = {
        "input_records": len(df),
        "final_usable_records": len(usable),
        "final_excluded_records": int((~df["final_usable"]).sum()),
        "usable_fraction": round(len(usable) / len(df), 4),
        "usable_center_kinases": int(usable["center_kinase_symbols"].str.split(";").explode().nunique()),
        "usable_pmids": int(usable["pmid"].nunique()),
        "target_entity_type_counts": usable["final_target_entity_type"].value_counts().to_dict(),
        "target_status_counts": usable["final_target_status"].value_counts().to_dict(),
        "exclusion_reason_counts": df.loc[~df["final_usable"], "final_exclusion_reason"].value_counts().to_dict(),
        "evidence_tier_counts_usable": usable["evidence_tier"].value_counts().to_dict(),
        "outputs": {
            "final_usable_parquet": str(usable_out),
            "final_usable_tsv": str(usable_tsv),
            "final_with_flags_parquet": str(full_out),
            "excluded_or_review_workbook": str(review_out),
        },
    }
    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "# Final Cleaning Report",
        "",
        f"- Input records: {len(df):,}",
        f"- Final usable records: {len(usable):,}",
        f"- Excluded records: {(~df['final_usable']).sum():,}",
        f"- Usable fraction: {len(usable) / len(df):.1%}",
        f"- Usable PMIDs: {usable['pmid'].nunique():,}",
        f"- Center kinases in usable set: {report['usable_center_kinases']:,}",
        "",
        "## Usable Target Entity Types",
    ]
    for key, value in report["target_entity_type_counts"].items():
        lines.append(f"- {key}: {value:,}")
    lines += ["", "## Exclusion Reasons"]
    for key, value in report["exclusion_reason_counts"].items():
        lines.append(f"- {key}: {value:,}")
    lines += [
        "",
        "## Outputs",
        f"- Final usable parquet: `{usable_out}`",
        f"- Final usable TSV: `{usable_tsv}`",
        f"- Full table with flags: `{full_out}`",
        f"- Excluded/review workbook: `{review_out}`",
    ]
    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
