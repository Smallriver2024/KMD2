"""
Step 0: 合并激酶结构域注释

Supplementary_Table_S2_latest.csv 作为人工/已有注释优先；
drkinase_results_new.csv 仅补充同一 UniProt + 同一结构域类型下
不与已有注释重叠的预测区间。

输出:
  output/kinase_domain_annotation_merged.txt
  output/kinase_domain_annotation_merged.xlsx
"""

import os
import re

import pandas as pd


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

SUPPLEMENTARY_PATH = os.path.join(BASE_DIR, "Supplementary_Table_S2_latest.csv")
PREDICTED_PATH = os.path.join(BASE_DIR, "drkinase_results_new.csv")
KINASE_BASE_PATH = os.path.join(BASE_DIR, "kinase_base_fixed.xlsx")

DOMAIN_COLUMNS = {
    "Gatekeeper": "Gatekeeper",
    "A-loop": "A-loop",
    "G-loop": "G-loop",
    "αC-helix": "aC-Helix",
}


def parse_domain_ranges(value):
    """Parse position/range strings such as '315', '380..404', '245...256'."""
    if pd.isna(value):
        return []

    text = str(value).strip()
    if not text or text == "[N/A]":
        return []

    ranges = []
    for part in re.split(r"\s*,\s*", text):
        part = part.strip()
        if not part or part == "[N/A]":
            continue

        match = re.match(r"^(\d+)\.{2,3}(\d+)$", part)
        if match:
            start, end = int(match.group(1)), int(match.group(2))
            ranges.append((min(start, end), max(start, end)))
            continue

        if part.isdigit():
            pos = int(part)
            ranges.append((pos, pos))
            continue

        raise ValueError(f"Cannot parse domain range: {part!r}")

    return ranges


def intervals_overlap(start_a, end_a, start_b, end_b):
    return start_a <= end_b and start_b <= end_a


def load_kinase_lookup():
    df_kin = pd.read_excel(
        KINASE_BASE_PATH,
        sheet_name="kinase_master",
        usecols=["gene_symbol", "uniprot_id"],
    ).drop_duplicates(subset=["uniprot_id"])

    return df_kin.set_index("uniprot_id")["gene_symbol"].to_dict()


def build_known_domains():
    df_sup = pd.read_csv(SUPPLEMENTARY_PATH)

    rows = []
    for _, row in df_sup.iterrows():
        entry = str(row["unprot_id"]).split("-")[0]
        for source_col, domain_type in DOMAIN_COLUMNS.items():
            for start, end in parse_domain_ranges(row[source_col]):
                rows.append(
                    {
                        "Entry": entry,
                        "Kinase": row["Kinase"],
                        "gene_id": row.get("gene_id"),
                        "hgnc_id": row.get("hgnc_id"),
                        "domain_type": domain_type,
                        "domain_start": start,
                        "domain_end": end,
                        "domain_score": 1.0,
                        "source": "Supplementary_Table_S2_latest",
                        "source_priority": 1,
                        "source_domain_type": source_col,
                    }
                )

    df_known = pd.DataFrame(rows)
    return df_known.drop_duplicates(
        subset=["Entry", "domain_type", "domain_start", "domain_end"]
    )


def build_predicted_domains(kinase_lookup):
    df_pred = pd.read_csv(PREDICTED_PATH)
    df_pred = df_pred.rename(
        columns={
            "DR type": "domain_type",
            "Start": "domain_start",
            "End": "domain_end",
            "Score": "domain_score",
        }
    )
    df_pred["Entry"] = df_pred["Entry"].astype(str).str.split("-").str[0]
    df_pred["domain_type"] = df_pred["domain_type"].replace({"αC-helix": "aC-Helix"})
    df_pred["domain_start"] = df_pred["domain_start"].astype(int)
    df_pred["domain_end"] = df_pred["domain_end"].astype(int)
    df_pred["Kinase"] = df_pred["Entry"].map(kinase_lookup)
    df_pred["gene_id"] = pd.NA
    df_pred["hgnc_id"] = pd.NA
    df_pred["source"] = "drkinase_results_new"
    df_pred["source_priority"] = 2
    df_pred["source_domain_type"] = df_pred["domain_type"]

    keep_cols = [
        "Entry",
        "Kinase",
        "gene_id",
        "hgnc_id",
        "domain_type",
        "domain_start",
        "domain_end",
        "domain_score",
        "source",
        "source_priority",
        "source_domain_type",
    ]
    df_pred = df_pred[keep_cols]
    df_pred = df_pred.dropna(subset=["Kinase"])
    return df_pred.drop_duplicates(
        subset=["Entry", "domain_type", "domain_start", "domain_end"],
        keep="first",
    )


def filter_non_overlapping_predictions(df_known, df_pred):
    known_by_key = {
        key: grp[["domain_start", "domain_end"]].to_records(index=False)
        for key, grp in df_known.groupby(["Entry", "domain_type"])
    }

    keep = []
    for _, row in df_pred.iterrows():
        intervals = known_by_key.get((row["Entry"], row["domain_type"]), [])
        overlaps_known = any(
            intervals_overlap(
                row["domain_start"],
                row["domain_end"],
                int(start),
                int(end),
            )
            for start, end in intervals
        )
        keep.append(not overlaps_known)

    return df_pred.loc[keep].copy()


def main():
    print("=" * 60)
    print("Step 0: 合并激酶结构域注释")
    print("=" * 60)

    kinase_lookup = load_kinase_lookup()
    df_known = build_known_domains()
    df_pred = build_predicted_domains(kinase_lookup)
    df_pred_keep = filter_non_overlapping_predictions(df_known, df_pred)

    df_merged = pd.concat([df_known, df_pred_keep], ignore_index=True)
    df_merged = df_merged.sort_values(
        ["Entry", "source_priority", "domain_type", "domain_start", "domain_end"]
    ).reset_index(drop=True)

    out_cols = [
        "Entry",
        "Kinase",
        "gene_id",
        "hgnc_id",
        "domain_type",
        "domain_start",
        "domain_end",
        "domain_score",
        "source",
        "source_priority",
        "source_domain_type",
    ]
    df_merged = df_merged[out_cols]

    txt_path = os.path.join(OUT_DIR, "kinase_domain_annotation_merged.txt")
    xlsx_path = os.path.join(OUT_DIR, "kinase_domain_annotation_merged.xlsx")

    df_merged.to_csv(txt_path, sep="\t", index=False)
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df_merged.to_excel(writer, sheet_name="merged_domains", index=False)

        summary = pd.DataFrame(
            {
                "metric": [
                    "known_domains",
                    "predicted_domains_input",
                    "predicted_domains_kept",
                    "merged_domains",
                    "merged_uniprot_entries",
                    "merged_kinases",
                ],
                "value": [
                    len(df_known),
                    len(df_pred),
                    len(df_pred_keep),
                    len(df_merged),
                    df_merged["Entry"].nunique(),
                    df_merged["Kinase"].nunique(),
                ],
            }
        )
        summary.to_excel(writer, sheet_name="summary", index=False)

        by_type_source = (
            df_merged.groupby(["domain_type", "source"], dropna=False)
            .size()
            .reset_index(name="n")
        )
        by_type_source.to_excel(writer, sheet_name="by_type_source", index=False)

    print(f"  已有注释结构域: {len(df_known):,}")
    print(f"  预测结构域输入: {len(df_pred):,}")
    print(f"  预测结构域保留: {len(df_pred_keep):,}")
    print(f"  合并后结构域:   {len(df_merged):,}")
    print("\n  合并后结构域分布:")
    print(df_merged["domain_type"].value_counts().to_string())
    print("\n  来源分布:")
    print(df_merged["source"].value_counts().to_string())
    print(f"\n  保存: {txt_path}")
    print(f"  保存: {xlsx_path}")


if __name__ == "__main__":
    main()
