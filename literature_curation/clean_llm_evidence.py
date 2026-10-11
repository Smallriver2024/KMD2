from __future__ import annotations

import ast
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

SOURCE = ROOT / "events_from2000_human_known_disease_cleaned_v2.parquet"
KINASE_FILE = ROOT / "kinase_base_fixed.xlsx"
HGNC_FILE = Path(__file__).resolve().parent / "reference" / "hgnc_complete_set.txt"


GREEK = {
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "κ": "kappa",
    "μ": "mu",
    "–": "-",
    "—": "-",
    "−": "-",
    "’": "'",
    "“": '"',
    "”": '"',
}


KINASE_ALIAS = {
    "AKT": ["AKT1", "AKT2", "AKT3"],
    "AKT-1": ["AKT1"],
    "AKT-2": ["AKT2"],
    "AKT-3": ["AKT3"],
    "PKB": ["AKT1", "AKT2", "AKT3"],
    "ERK": ["MAPK1", "MAPK3"],
    "ERK1": ["MAPK3"],
    "ERK2": ["MAPK1"],
    "ERK1/2": ["MAPK3", "MAPK1"],
    "ERK1 AND ERK2": ["MAPK3", "MAPK1"],
    "P38": ["MAPK14"],
    "P38 MAPK": ["MAPK14"],
    "P38 MAP KINASE": ["MAPK14"],
    "JNK": ["MAPK8", "MAPK9", "MAPK10"],
    "JNK1": ["MAPK8"],
    "JNK2": ["MAPK9"],
    "JNK3": ["MAPK10"],
    "MTORC1": ["MTOR"],
    "MTORC2": ["MTOR"],
    "PDK1": ["PDPK1"],
    "AMPK": ["PRKAA1", "PRKAA2"],
    "AMPKALPHA": ["PRKAA1", "PRKAA2"],
    "AMPK-ALPHA": ["PRKAA1", "PRKAA2"],
    "PKA": ["PRKACA", "PRKACB", "PRKACG"],
    "PKC": ["PRKCA", "PRKCB", "PRKCG"],
    "PKC-ALPHA": ["PRKCA"],
    "PKC-BETA": ["PRKCB"],
    "PKC-GAMMA": ["PRKCG"],
    "MEK": ["MAP2K1", "MAP2K2"],
    "MEK1": ["MAP2K1"],
    "MEK2": ["MAP2K2"],
    "S6K": ["RPS6KB1", "RPS6KB2"],
    "P70S6K": ["RPS6KB1"],
    "IKK": ["CHUK", "IKBKB", "IKBKE"],
    "IKKALPHA": ["CHUK"],
    "IKKBETA": ["IKBKB"],
    "CK2": ["CSNK2A1", "CSNK2A2"],
    "CDK": ["CDK1", "CDK2", "CDK4", "CDK6"],
    "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
    "PI 3-KINASE": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
}


NON_KINASE_ALIAS = {
    "P53": ("TP53", "P04637", "tumor protein p53"),
    "TAU": ("MAPT", "P10636", "microtubule-associated protein tau"),
    "C-JUN": ("JUN", "P05412", "Jun proto-oncogene"),
    "JUN": ("JUN", "P05412", "Jun proto-oncogene"),
    "CREB": ("CREB1", "P16220", "cAMP responsive element binding protein 1"),
    "BETA-CATENIN": ("CTNNB1", "P35222", "beta-catenin"),
    "Β-CATENIN": ("CTNNB1", "P35222", "beta-catenin"),
    "CTNNB1": ("CTNNB1", "P35222", "beta-catenin"),
    "P65": ("RELA", "Q04206", "NF-kappa-B p65 subunit"),
    "NF-KB P65": ("RELA", "Q04206", "NF-kappa-B p65 subunit"),
    "IKAPPABALPHA": ("NFKBIA", "P25963", "NF-kappa-B inhibitor alpha"),
    "IKBΑ": ("NFKBIA", "P25963", "NF-kappa-B inhibitor alpha"),
    "IKBALPHA": ("NFKBIA", "P25963", "NF-kappa-B inhibitor alpha"),
    "RB": ("RB1", "P06400", "retinoblastoma protein"),
    "RETINOBLASTOMA PROTEIN": ("RB1", "P06400", "retinoblastoma protein"),
    "HSP27": ("HSPB1", "P04792", "heat shock protein beta-1"),
    "BAD": ("BAD", "Q92934", "BCL2 associated agonist of cell death"),
    "FOXO3A": ("FOXO3", "O43524", "forkhead box O3"),
    "FOXO1": ("FOXO1", "Q12778", "forkhead box O1"),
    "GSK3": ("GSK3B", "P49841", "glycogen synthase kinase 3 beta"),
}


AGENT_CANONICAL = {
    "TNF-ALPHA": "TNF-alpha",
    "TNF-Α": "TNF-alpha",
    "TGF-BETA1": "TGF-beta1",
    "TGF-B1": "TGF-beta1",
    "TGF-Β1": "TGF-beta1",
    "IL-6": "IL-6",
    "INTERLEUKIN-6": "IL-6",
    "IL-8": "IL-8",
    "INTERLEUKIN-8": "IL-8",
    "EGF": "EGF",
    "HGF": "HGF",
    "VEGF": "VEGF",
    "INSULIN": "insulin",
    "HIGH GLUCOSE": "high glucose",
    "HYPOXIA": "hypoxia",
    "IONIZING RADIATION": "ionizing radiation",
}


DISEASE_CANONICAL = {
    "NON-SMALL CELL LUNG CARCINOMA": "Non-Small Cell Lung Cancer",
    "NSCLC": "Non-Small Cell Lung Cancer",
    "HEPATOCELLULAR CARCINOMA": "Hepatocellular Carcinoma",
    "HCC": "Hepatocellular Carcinoma",
    "COLORECTAL CARCINOMA": "Colorectal Cancer",
    "COLON CARCINOMA": "Colon Cancer",
    "ALZHEIMER DISEASE": "Alzheimer's Disease",
    "ALZHEIMER'S DISEASE": "Alzheimer's Disease",
    "TYPE II DIABETES": "Type 2 Diabetes",
    "TYPE 2 DIABETES MELLITUS": "Type 2 Diabetes",
}


def clean_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null", "n/a"}:
        return ""
    for old, new in GREEK.items():
        text = text.replace(old, new)
    text = re.sub(r"\s+", " ", text)
    text = text.strip(" \t\r\n;")
    return text


def key_text(value: object) -> str:
    text = clean_text(value)
    text = text.replace("β", "beta")
    text = text.replace("κ", "kappa")
    text = re.sub(r"\s+", " ", text)
    return text.upper()


def split_entities(value: object, split_and: bool = False) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    if len(text) > 120 and not re.search(r"[/;]", text):
        return [text]
    sep = r"\s*(?:/|;)\s*"
    if "," in text and re.search(r"\b[A-Z0-9]{2,}[,;]", text):
        sep = r"\s*(?:/|;|,)\s*"
    if split_and:
        sep = r"\s*(?:/|;|,|\band\b)\s*"
    parts = [p.strip() for p in re.split(sep, text) if p.strip()]
    return list(dict.fromkeys(parts))


def parse_listish_sites(value: object) -> list[str]:
    text = clean_text(value)
    if not text or text == "[]":
        return []
    return list(dict.fromkeys(re.findall(r"[STY]\d+", text.upper())))


def split_pipe_list(value: object) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    return [x.strip() for x in text.split("|") if x.strip()]


def load_hgnc_map() -> dict[str, dict]:
    if not HGNC_FILE.exists():
        return {}
    hgnc = pd.read_csv(HGNC_FILE, sep="\t", dtype=str).fillna("")
    mapping: dict[str, dict] = {}
    for _, row in hgnc.iterrows():
        symbol = clean_text(row.get("symbol"))
        if not symbol:
            continue
        rec = {
            "symbol": symbol,
            "name": clean_text(row.get("name")),
            "hgnc_id": clean_text(row.get("hgnc_id")),
            "uniprot_id": split_pipe_list(row.get("uniprot_ids"))[0] if split_pipe_list(row.get("uniprot_ids")) else "",
        }
        keys = [symbol]
        keys += split_pipe_list(row.get("alias_symbol"))
        keys += split_pipe_list(row.get("prev_symbol"))
        keys += [rec["name"]]
        for key in keys:
            k = key_text(key)
            if not k:
                continue
            # Prefer exact approved symbols over aliases if collisions occur.
            if k not in mapping or k == key_text(symbol):
                mapping[k] = rec
            compact = re.sub(r"[^A-Z0-9]", "", k)
            if compact and compact not in mapping:
                mapping[compact] = rec
    return mapping


def load_known_drugs() -> dict[str, str]:
    names: dict[str, str] = {}
    drug_file = ROOT / "all_drug_response.txt"
    resistance_file = ROOT / "drug_resistance.xlsx"
    if drug_file.exists():
        for chunk in pd.read_csv(drug_file, sep="\t", usecols=["drug"], dtype=str, chunksize=200000):
            for drug in chunk["drug"].dropna().astype(str).unique():
                names.setdefault(key_text(drug), clean_text(drug))
    if resistance_file.exists():
        try:
            dr = pd.read_excel(resistance_file, usecols=["Drug"], dtype=str)
            for drug in dr["Drug"].dropna().astype(str).unique():
                names.setdefault(key_text(drug), clean_text(drug))
        except Exception:
            pass
    return names


def load_kinase_maps() -> tuple[pd.DataFrame, dict[str, dict], dict[str, list[str]]]:
    km = pd.read_excel(KINASE_FILE, sheet_name="kinase_master")
    by_symbol = {}
    aliases = defaultdict(list)
    for _, row in km.iterrows():
        symbol = str(row["gene_symbol"]).strip()
        rec = {
            "symbol": symbol,
            "uniprot_id": str(row["uniprot_id"]).strip(),
            "kinase_id": str(row["kinase_id"]).strip(),
            "group_name": clean_text(row.get("group_name")),
            "family_name": clean_text(row.get("family_name")),
        }
        by_symbol[symbol] = rec
        aliases[key_text(symbol)].append(symbol)
        protein = clean_text(row.get("protein_name"))
        for alias in re.findall(r"\(([^()]{2,45})\)", protein):
            if alias.startswith("EC ") or alias.startswith("CD "):
                continue
            aliases[key_text(alias)].append(symbol)
    for alias, symbols in KINASE_ALIAS.items():
        aliases[key_text(alias)].extend(symbols)
    aliases = {k: list(dict.fromkeys(v)) for k, v in aliases.items()}
    return km, by_symbol, aliases


def normalize_kinase_field(value: object, aliases: dict[str, list[str]], by_symbol: dict[str, dict]) -> dict:
    raw = clean_text(value)
    parts = split_entities(raw, split_and=True)
    symbols: list[str] = []
    unmatched: list[str] = []
    for part in parts:
        k = key_text(part)
        k_compact = re.sub(r"[^A-Z0-9]", "", k)
        found = aliases.get(k) or aliases.get(k_compact) or []
        if not found and k in by_symbol:
            found = [k]
        if found:
            symbols.extend(found)
        else:
            unmatched.append(part)
    symbols = list(dict.fromkeys([s for s in symbols if s in by_symbol]))
    ids = [by_symbol[s]["uniprot_id"] for s in symbols]
    if not raw:
        status = "blank"
    elif symbols and not unmatched:
        status = "mapped"
    elif symbols and unmatched:
        status = "partial_review"
    else:
        status = "unmapped_review"
    return {
        "raw": raw,
        "symbols": symbols,
        "uniprot_ids": ids,
        "unmatched": unmatched,
        "status": status,
    }


def normalize_target(value: object, role: object, aliases: dict[str, list[str]], by_symbol: dict[str, dict], hgnc_map: dict[str, dict]) -> dict:
    raw = clean_text(value)
    kinase = normalize_kinase_field(raw, aliases, by_symbol)
    non_kinase_symbols: list[str] = []
    non_kinase_uniprot: list[str] = []
    non_kinase_unmatched: list[str] = []
    if not kinase["symbols"]:
        for part in split_entities(raw, split_and=False):
            key = key_text(part)
            key_compact = re.sub(r"[^A-Z0-9]", "", key)
            hit = NON_KINASE_ALIAS.get(key) or NON_KINASE_ALIAS.get(key_compact)
            if hit:
                non_kinase_symbols.append(hit[0])
                non_kinase_uniprot.append(hit[1])
            elif key in hgnc_map or key_compact in hgnc_map:
                rec = hgnc_map.get(key) or hgnc_map.get(key_compact)
                non_kinase_symbols.append(rec["symbol"])
                if rec.get("uniprot_id"):
                    non_kinase_uniprot.append(rec["uniprot_id"])
            elif part:
                non_kinase_unmatched.append(part)
    if kinase["symbols"]:
        status = kinase["status"]
    elif non_kinase_symbols and not non_kinase_unmatched:
        status = "mapped_non_kinase"
    elif non_kinase_symbols:
        status = "partial_non_kinase_review"
    elif raw:
        status = "unmapped_non_kinase_review"
    else:
        status = "blank"
    return {
        "raw": raw,
        "role": clean_text(role),
        "kinase_symbols": kinase["symbols"],
        "kinase_uniprot_ids": kinase["uniprot_ids"],
        "non_kinase_symbols": list(dict.fromkeys(non_kinase_symbols)),
        "non_kinase_uniprot_ids": list(dict.fromkeys(non_kinase_uniprot)),
        "unmatched": kinase["unmatched"] or non_kinase_unmatched,
        "status": status,
    }


def normalize_agent(value: object, known_drugs: dict[str, str]) -> dict:
    raw = clean_text(value)
    if not raw:
        return {"clean": "", "status": "blank", "review_reason": ""}
    cleaned = re.sub(r"^(treatment with|stimulation with|exposure to|treated with)\s+", "", raw, flags=re.I)
    key = key_text(cleaned)
    canonical = AGENT_CANONICAL.get(key) or known_drugs.get(key) or cleaned
    review = []
    if len(cleaned) > 90:
        review.append("long_descriptive_agent")
    if re.search(r"\s(?:and|or)\s|[,;/]", cleaned, flags=re.I) and len(split_entities(cleaned, split_and=True)) > 1:
        review.append("multi_agent_or_context")
    if re.search(r"\b(inhibitors?|activators?|conditioned medium|extract|mixture|combination)\b", cleaned, flags=re.I):
        review.append("class_or_mixture")
    return {
        "clean": canonical,
        "status": "review" if review else "cleaned",
        "review_reason": ";".join(review),
    }


def normalize_disease(value: object) -> dict:
    raw = clean_text(value)
    if not raw:
        return {"clean": "", "status": "blank", "review_reason": ""}
    key = key_text(raw)
    canonical = DISEASE_CANONICAL.get(key)
    if not canonical:
        canonical = " ".join(w if w.isupper() and len(w) <= 5 else w.capitalize() for w in raw.split())
        canonical = canonical.replace("And", "and").replace("Of", "of").replace("With", "with")
        canonical = canonical.replace("Hiv", "HIV").replace("Dna", "DNA").replace("Rna", "RNA")
    review = []
    if len(raw) > 95:
        review.append("long_multi_disease")
    if re.search(r";|,", raw, flags=re.I):
        review.append("possible_multi_disease")
    return {
        "clean": canonical,
        "status": "review" if review else "cleaned",
        "review_reason": ";".join(review),
    }


def parse_review_reasons(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, float) and math.isnan(value):
        return []
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    text = str(value).strip()
    if not text or text == "[]":
        return []
    try:
        parsed = ast.literal_eval(text)
        if isinstance(parsed, list):
            return [str(v) for v in parsed]
    except Exception:
        pass
    return [x.strip().strip("\"'") for x in re.split(r"[,;]", text.strip("[]")) if x.strip()]


def evidence_tier(row: pd.Series) -> str:
    has_center = bool(row["center_kinase_symbols"])
    if not has_center:
        return "no_center_kinase_review"
    if row["needs_manual_review"] == "no" and int(row["confidence"]) >= 8:
        return "curated_candidate"
    if int(row["confidence"]) >= 8:
        return "literature_mined_review"
    return "low_confidence_review"


def join(values: list[str]) -> str:
    return ";".join(values)


def main() -> None:
    _, by_symbol, aliases = load_kinase_maps()
    hgnc_map = load_hgnc_map()
    known_drugs = load_known_drugs()
    df = pd.read_parquet(SOURCE).copy()
    cleaned_rows = []
    review_rows = []

    for idx, row in df.iterrows():
        upstream = normalize_kinase_field(row.get("upstream_kinase"), aliases, by_symbol)
        target = normalize_target(row.get("phospho_target"), row.get("phospho_target_role"), aliases, by_symbol, hgnc_map)
        agent = normalize_agent(row.get("primary_agent"), known_drugs)
        disease = normalize_disease(row.get("disease"))
        sites = parse_listish_sites(row.get("phosphosites"))

        center_symbols: list[str] = []
        center_uniprots: list[str] = []
        if clean_text(row.get("phospho_target_role")) == "kinase":
            center_symbols.extend(target["kinase_symbols"])
            center_uniprots.extend(target["kinase_uniprot_ids"])
        center_symbols.extend(upstream["symbols"])
        center_uniprots.extend(upstream["uniprot_ids"])
        center_symbols = list(dict.fromkeys(center_symbols))
        center_uniprots = list(dict.fromkeys(center_uniprots))

        out = row.to_dict()
        out.update(
            {
                "phospho_target_clean": target["raw"],
                "phospho_target_kinase_symbols": join(target["kinase_symbols"]),
                "phospho_target_kinase_uniprot_ids": join(target["kinase_uniprot_ids"]),
                "phospho_target_non_kinase_symbols": join(target["non_kinase_symbols"]),
                "phospho_target_non_kinase_uniprot_ids": join(target["non_kinase_uniprot_ids"]),
                "phospho_target_norm_status": target["status"],
                "phospho_target_unmatched": join(target["unmatched"]),
                "phosphosites_clean": join(sites),
                "upstream_kinase_clean": upstream["raw"],
                "upstream_kinase_symbols": join(upstream["symbols"]),
                "upstream_kinase_uniprot_ids": join(upstream["uniprot_ids"]),
                "upstream_kinase_norm_status": upstream["status"],
                "upstream_kinase_unmatched": join(upstream["unmatched"]),
                "primary_agent_clean": agent["clean"],
                "primary_agent_norm_status": agent["status"],
                "primary_agent_review_reason": agent["review_reason"],
                "disease_clean": disease["clean"],
                "disease_norm_status": disease["status"],
                "disease_review_reason": disease["review_reason"],
                "center_kinase_symbols": join(center_symbols),
                "center_kinase_uniprot_ids": join(center_uniprots),
            }
        )
        out["evidence_tier"] = evidence_tier(pd.Series(out))
        cleaned_rows.append(out)

        reasons = set(parse_review_reasons(row.get("review_reasons")))
        for status_col, reason in [
            ("phospho_target_norm_status", "target_" + out["phospho_target_norm_status"]),
            ("upstream_kinase_norm_status", "upstream_" + out["upstream_kinase_norm_status"]),
            ("primary_agent_norm_status", "agent_" + out["primary_agent_norm_status"]),
            ("disease_norm_status", "disease_" + out["disease_norm_status"]),
        ]:
            if "review" in out[status_col]:
                reasons.add(reason)
        if out["evidence_tier"].endswith("review") or reasons:
            review_rows.append(
                {
                    "source_row": idx,
                    "pmid": out["pmid"],
                    "pub_year": out["pub_year"],
                    "center_kinase_symbols": out["center_kinase_symbols"],
                    "center_kinase_uniprot_ids": out["center_kinase_uniprot_ids"],
                    "phospho_target": row.get("phospho_target"),
                    "phospho_target_clean": out["phospho_target_clean"],
                    "phospho_target_norm_status": out["phospho_target_norm_status"],
                    "upstream_kinase": row.get("upstream_kinase"),
                    "upstream_kinase_symbols": out["upstream_kinase_symbols"],
                    "upstream_kinase_norm_status": out["upstream_kinase_norm_status"],
                    "primary_agent": row.get("primary_agent"),
                    "primary_agent_clean": out["primary_agent_clean"],
                    "disease": row.get("disease"),
                    "disease_clean": out["disease_clean"],
                    "phosphosites_clean": out["phosphosites_clean"],
                    "regulation_mode": out["regulation_mode"],
                    "confidence": out["confidence"],
                    "needs_manual_review": out["needs_manual_review"],
                    "evidence_tier": out["evidence_tier"],
                    "review_reasons_combined": ";".join(sorted(reasons)),
                    "evidence_text": out["evidence_text"],
                    "relationship": out["relationship"],
                }
            )

    cleaned = pd.DataFrame(cleaned_rows)
    review = pd.DataFrame(review_rows)

    cleaned_parquet = OUT / "llm_phosphorylation_evidence_cleaned_kinase_centered.parquet"
    cleaned_tsv = OUT / "llm_phosphorylation_evidence_cleaned_kinase_centered.tsv"
    review_csv = OUT / "manual_review_rows_full.csv"
    review_xlsx = OUT / "manual_review_workbook.xlsx"
    report_json = OUT / "normalization_report.json"
    report_md = OUT / "normalization_report.md"

    cleaned.to_parquet(cleaned_parquet, index=False)
    cleaned.to_csv(cleaned_tsv, sep="\t", index=False)
    review.to_csv(review_csv, index=False)

    def entity_review(field_raw: str, field_clean: str, status_col: str, extra_cols: list[str], top_n: int = 1200) -> pd.DataFrame:
        cols = [field_raw, field_clean, status_col] + extra_cols
        tmp = cleaned[cols + ["pmid"]].copy()
        tmp[field_raw] = tmp[field_raw].fillna("").astype(str)
        grouped = (
            tmp.groupby(cols, dropna=False)
            .agg(record_count=("pmid", "size"), pmid_examples=("pmid", lambda x: ";".join(map(str, list(x)[:5]))))
            .reset_index()
            .sort_values(["record_count", field_raw], ascending=[False, True])
        )
        return grouped[grouped[status_col].astype(str).str.contains("review|unmapped|partial", regex=True)].head(top_n)

    upstream_review = entity_review("upstream_kinase", "upstream_kinase_clean", "upstream_kinase_norm_status", ["upstream_kinase_unmatched"], 900)
    target_review = entity_review("phospho_target", "phospho_target_clean", "phospho_target_norm_status", ["phospho_target_unmatched"], 1200)
    agent_review = entity_review("primary_agent", "primary_agent_clean", "primary_agent_norm_status", ["primary_agent_review_reason"], 1500)
    disease_review = entity_review("disease", "disease_clean", "disease_norm_status", ["disease_review_reason"], 1200)

    with pd.ExcelWriter(review_xlsx, engine="xlsxwriter") as writer:
        summary = pd.DataFrame(
            [
                {"metric": "source_records", "value": len(df)},
                {"metric": "cleaned_records", "value": len(cleaned)},
                {"metric": "manual_review_rows", "value": len(review)},
                {"metric": "curated_candidate", "value": int((cleaned["evidence_tier"] == "curated_candidate").sum())},
                {"metric": "records_with_center_kinase", "value": int((cleaned["center_kinase_symbols"] != "").sum())},
                {"metric": "records_with_clean_sites", "value": int((cleaned["phosphosites_clean"] != "").sum())},
            ]
        )
        summary.to_excel(writer, sheet_name="summary", index=False)
        upstream_review.to_excel(writer, sheet_name="upstream_kinase_review", index=False)
        target_review.to_excel(writer, sheet_name="target_review", index=False)
        agent_review.to_excel(writer, sheet_name="agent_context_review", index=False)
        disease_review.to_excel(writer, sheet_name="disease_review", index=False)
        review.head(5000).to_excel(writer, sheet_name="row_review_top5000", index=False)
        for sheet in writer.sheets.values():
            sheet.freeze_panes(1, 0)

    report = {
        "source_records": len(df),
        "cleaned_records": len(cleaned),
        "manual_review_rows": len(review),
        "evidence_tier_counts": cleaned["evidence_tier"].value_counts().to_dict(),
        "records_with_center_kinase": int((cleaned["center_kinase_symbols"] != "").sum()),
        "records_with_clean_sites": int((cleaned["phosphosites_clean"] != "").sum()),
        "hgnc_reference_loaded": bool(hgnc_map),
        "known_drug_names_loaded": len(known_drugs),
        "upstream_kinase_norm_status": cleaned["upstream_kinase_norm_status"].value_counts().to_dict(),
        "phospho_target_norm_status": cleaned["phospho_target_norm_status"].value_counts().to_dict(),
        "primary_agent_norm_status": cleaned["primary_agent_norm_status"].value_counts().to_dict(),
        "disease_norm_status": cleaned["disease_norm_status"].value_counts().to_dict(),
        "top_center_kinases": cleaned.loc[cleaned["center_kinase_symbols"] != "", "center_kinase_symbols"].str.split(";").explode().value_counts().head(30).to_dict(),
        "outputs": {
            "cleaned_parquet": str(cleaned_parquet),
            "cleaned_tsv": str(cleaned_tsv),
            "manual_review_csv": str(review_csv),
            "manual_review_xlsx": str(review_xlsx),
        },
    }
    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "# Kinase-Centered LLM Evidence Normalization Report",
        "",
        f"- Source records: {len(df):,}",
        f"- Cleaned records: {len(cleaned):,}",
        f"- Records with center kinase: {report['records_with_center_kinase']:,}",
        f"- Records with cleaned phosphosite: {report['records_with_clean_sites']:,}",
        f"- Manual review rows: {len(review):,}",
        "",
        "## Evidence Tiers",
    ]
    for k, v in report["evidence_tier_counts"].items():
        lines.append(f"- {k}: {v:,}")
    lines += [
        "",
        "## Output Files",
        f"- Cleaned parquet: `{cleaned_parquet}`",
        f"- Cleaned TSV: `{cleaned_tsv}`",
        f"- Full manual review CSV: `{review_csv}`",
        f"- Manual review workbook: `{review_xlsx}`",
    ]
    report_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
