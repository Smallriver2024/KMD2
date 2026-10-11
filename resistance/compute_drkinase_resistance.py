from __future__ import annotations

import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


# Reference-only code snapshot: override this at runtime instead of exposing a
# contributor's local path. The surrounding source data are intentionally absent.
ROOT = Path(os.environ.get("KINASEMD_DATA_ROOT", Path(__file__).resolve().parents[2]))
OUT_DIR = ROOT / "outputs"

MUTATION_PATH = ROOT / "mutation_data.tsv"
RESPONSE_PATH = ROOT / "drug_response_data.tsv"
PREDICTION_PATH = ROOT / "predicted_hotspots.csv"
MODEL_PATH = ROOT / "model_metadata.csv"
OLD_RESISTANCE_PATH = ROOT / "reference_resistance.xlsx"

KEYS = ["Kinase", "Drug", "Cancer code", "Cancer name", "Sub-tissue"]


def clean_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if not np.isfinite(value):
            return None
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def records_for_json(df: pd.DataFrame) -> list[dict]:
    return [
        {column: clean_value(value) for column, value in row.items()}
        for row in df.to_dict(orient="records")
    ]


def normalize_name(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.upper()
        .str.replace(r"[^A-Z0-9]", "", regex=True)
    )


def treatment_component(value: object) -> str:
    text = "Unknown" if pd.isna(value) or str(value).strip() == "" else str(value).strip()
    text = re.sub(r"\s+", "_", text)
    return text.replace(":", "-")


def fallback_code(name: str) -> str:
    words = re.findall(r"[A-Za-z0-9]+", str(name))
    if not words:
        return "UNK"
    if len(words) == 1:
        return words[0][:5].upper()
    return "".join(word[0] for word in words[:5]).upper()


def load_cancer_code_map() -> dict[str, str]:
    old = pd.read_excel(OLD_RESISTANCE_PATH, usecols=["Cancer code", "Cancer name"])
    old = old.dropna().drop_duplicates("Cancer name")
    return dict(zip(old["Cancer name"].astype(str), old["Cancer code"].astype(str)))


def load_model_metadata() -> pd.DataFrame:
    columns = [
        "model_id",
        "model_name",
        "tissue",
        "cancer_type",
        "cancer_type_ncit_id",
        "cancer_type_detail",
        "BROAD_ID",
        "COSMIC_ID",
    ]
    model = pd.read_csv(MODEL_PATH, usecols=columns)
    model["_has_broad"] = model["BROAD_ID"].notna()
    model["_has_model_id"] = model["model_id"].notna()
    model = (
        model.sort_values(
            ["model_name", "_has_broad", "_has_model_id"],
            ascending=[True, False, False],
        )
        .drop_duplicates("model_name", keep="first")
        .reset_index(drop=True)
    )
    return model


def map_mutations_to_predicted_hotspots(model: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    usecols = [
        "dataset",
        "kinase",
        "gene_id",
        "hgnc_id",
        "cancer",
        "sample_id",
        "AAchange",
        "mut_type",
        "AApos",
        "AAx",
        "uniprot_id",
    ]
    mutation = pd.read_csv(MUTATION_PATH, sep="\t", usecols=usecols)
    mutation = mutation[mutation["dataset"].eq("GDSC")].copy()
    mutation["Entry"] = mutation["uniprot_id"].astype(str).str.split("-").str[0]
    mutation["AApos"] = pd.to_numeric(mutation["AApos"], errors="coerce")
    mutation = mutation.dropna(subset=["AApos"])
    mutation["AApos"] = mutation["AApos"].astype(int)

    predicted = pd.read_csv(PREDICTION_PATH).rename(
        columns={
            "DR type": "Hotspot",
            "Start": "Hotspot start",
            "End": "Hotspot end",
            "Score": "Hotspot score",
        }
    )
    predicted["Entry"] = predicted["Entry"].astype(str).str.split("-").str[0]

    frames = []
    for entry, mutation_group in mutation.groupby("Entry", sort=False):
        intervals = predicted[predicted["Entry"].eq(entry)]
        if intervals.empty:
            continue
        for _, interval in intervals.iterrows():
            start = int(interval["Hotspot start"])
            end = int(interval["Hotspot end"])
            hit = mutation_group[mutation_group["AApos"].between(start, end)].copy()
            if hit.empty:
                continue
            hit["Hotspot"] = interval["Hotspot"]
            hit["Hotspot start"] = start
            hit["Hotspot end"] = end
            hit["Hotspot score"] = float(interval["Hotspot score"])
            hit["Hotspot source"] = "DrKinase predicted hotspot"
            frames.append(hit)

    hotspot = pd.concat(frames, ignore_index=True)
    hotspot = (
        hotspot.sort_values("Hotspot score", ascending=False)
        .drop_duplicates(["sample_id", "kinase", "AAchange", "Hotspot"], keep="first")
        .reset_index(drop=True)
    )

    hotspot = hotspot.merge(
        model.drop(columns=["_has_broad", "_has_model_id"]),
        left_on="sample_id",
        right_on="model_name",
        how="left",
    )

    response_ids = pd.read_csv(
        RESPONSE_PATH,
        sep="\t",
        usecols=["dataset", "cell_id", "cell_name"],
    )
    response_ids = response_ids[response_ids["dataset"].isin(["GDSC1", "GDSC2"])]
    response_ids = response_ids[["cell_id", "cell_name"]].drop_duplicates()
    response_ids["cell_name_norm"] = normalize_name(response_ids["cell_name"])
    name_counts = response_ids.groupby("cell_name_norm")["cell_id"].nunique()
    unique_names = set(name_counts[name_counts.eq(1)].index)
    name_to_id = (
        response_ids[response_ids["cell_name_norm"].isin(unique_names)]
        .drop_duplicates("cell_name_norm")
        .set_index("cell_name_norm")["cell_id"]
        .to_dict()
    )
    valid_response_ids = set(response_ids["cell_id"].astype(str))

    hotspot["response_cell_id"] = hotspot["BROAD_ID"]
    direct_mask = hotspot["response_cell_id"].isna() & hotspot["sample_id"].astype(str).isin(valid_response_ids)
    hotspot.loc[direct_mask, "response_cell_id"] = hotspot.loc[direct_mask, "sample_id"]
    fallback_mask = hotspot["response_cell_id"].isna()
    hotspot.loc[fallback_mask, "response_cell_id"] = normalize_name(
        hotspot.loc[fallback_mask, "sample_id"]
    ).map(name_to_id)

    qc = {
        "Raw GDSC mutation rows with AA position": int(len(mutation)),
        "Predicted hotspot intervals": int(len(predicted)),
        "Mapped mutation-hotspot rows": int(len(hotspot)),
        "Mapped mutation-hotspot rows with response cell ID": int(hotspot["response_cell_id"].notna().sum()),
        "Unmatched mutation-hotspot rows": int(hotspot["response_cell_id"].isna().sum()),
    }
    hotspot = hotspot.dropna(subset=["response_cell_id"]).copy()
    hotspot = hotspot.drop_duplicates(
        ["response_cell_id", "kinase", "AAchange", "Hotspot"], keep="first"
    )
    return hotspot, qc


def load_response(model: pd.DataFrame, cancer_codes: dict[str, str]) -> tuple[pd.DataFrame, dict]:
    response = pd.read_csv(RESPONSE_PATH, sep="\t")
    response = response[response["dataset"].isin(["GDSC1", "GDSC2"])].copy()
    response["ic50"] = pd.to_numeric(response["ic50"], errors="coerce")
    before_filter = len(response)
    response = response[response["ic50"].gt(0) & response["ic50"].le(10000)].copy()

    model_by_broad = (
        model.dropna(subset=["BROAD_ID"])
        .sort_values(["BROAD_ID", "_has_model_id"], ascending=[True, False])
        .drop_duplicates("BROAD_ID")
        .drop(columns=["_has_broad", "_has_model_id"])
    )
    response = response.merge(
        model_by_broad,
        left_on="cell_id",
        right_on="BROAD_ID",
        how="left",
    )
    response["Cancer name"] = response["cancer_type"].fillna(response["cancer"])
    response["Sub-tissue"] = response["cancer_type_detail"].fillna(response["Cancer name"])
    response["Tissue"] = response["tissue"].fillna("Unknown")
    response["Cancer code"] = response["Cancer name"].map(cancer_codes)
    response["Cancer code"] = response["Cancer code"].fillna(
        response["Cancer name"].map(fallback_code)
    )
    response["Cell line model ID"] = response["model_id"]
    response["Cell line name"] = response["cell_name"]
    response["Cell line ID"] = response["cell_id"]
    response["IC50 (uM)"] = response["ic50"]
    response["Kinase"] = response["kinase"]
    response["Drug"] = response["drug"]

    qc = {
        "Raw GDSC drug-response rows": int(before_filter),
        "Valid GDSC drug-response rows": int(len(response)),
        "Response rows with model metadata": int(response["model_id"].notna().sum()),
        "Response rows without model metadata": int(response["model_id"].isna().sum()),
        "GDSC1 response rows": int(response["dataset"].eq("GDSC1").sum()),
        "GDSC2 response rows": int(response["dataset"].eq("GDSC2").sum()),
    }
    return response, qc


def build_cell_candidates(response: pd.DataFrame, hotspot: pd.DataFrame) -> pd.DataFrame:
    group_keys = KEYS
    cell_stats = (
        response.groupby(
            group_keys + ["Cell line ID", "Cell line name", "Cell line model ID", "Tissue"],
            dropna=False,
        )
        .agg(
            **{
                "Median(IC50, G1)": ("IC50 (uM)", "median"),
                "G1 response records": ("IC50 (uM)", "size"),
                "Source datasets": ("dataset", lambda values: "; ".join(sorted(set(values)))),
            }
        )
        .reset_index()
    )

    hotspot_cells = hotspot[["kinase", "response_cell_id"]].drop_duplicates().rename(
        columns={"kinase": "Kinase", "response_cell_id": "Cell line ID"}
    )
    hotspot_cells["Has predicted hotspot mutation"] = True
    cell_stats = cell_stats.merge(hotspot_cells, on=["Kinase", "Cell line ID"], how="left")
    has_hotspot = cell_stats["Has predicted hotspot mutation"]
    cell_stats["Has predicted hotspot mutation"] = has_hotspot.where(
        has_hotspot.notna(), False
    ).astype(bool)

    background = (
        response.groupby(group_keys, dropna=False)
        .agg(
            **{
                "Median(IC50, treatment background)": ("IC50 (uM)", "median"),
                "Treatment response records": ("IC50 (uM)", "size"),
                "Treatment cell lines": ("Cell line ID", "nunique"),
            }
        )
        .reset_index()
    )

    strict_source = response.merge(
        hotspot_cells,
        on=["Kinase", "Cell line ID"],
        how="left",
    )
    strict_source = strict_source[strict_source["Has predicted hotspot mutation"].isna()]
    strict = (
        strict_source.groupby(group_keys, dropna=False)
        .agg(
            **{
                "Median(IC50, non-hotspot G2)": ("IC50 (uM)", "median"),
                "Non-hotspot G2 response records": ("IC50 (uM)", "size"),
                "Non-hotspot G2 cell lines": ("Cell line ID", "nunique"),
            }
        )
        .reset_index()
    )

    candidates = cell_stats[cell_stats["Has predicted hotspot mutation"]].copy()
    candidates = candidates.merge(background, on=group_keys, how="left")
    candidates = candidates.merge(strict, on=group_keys, how="left")
    candidates["FC(G1/treatment background)"] = (
        candidates["Median(IC50, G1)"] / candidates["Median(IC50, treatment background)"]
    )
    candidates["FC(G1/non-hotspot G2)"] = (
        candidates["Median(IC50, G1)"] / candidates["Median(IC50, non-hotspot G2)"]
    )
    candidates["Pass S1-like FC >= 1.5"] = candidates[
        "FC(G1/treatment background)"
    ].ge(1.5)
    candidates["Pass strict-method FC >= 1.5"] = candidates[
        "FC(G1/non-hotspot G2)"
    ].ge(1.5)
    return candidates


def add_treatment_fields(df: pd.DataFrame, prefix: str, baseline: str) -> pd.DataFrame:
    out = df.copy()
    treatment_values = []
    for _, row in out.iterrows():
        treatment_values.append(
            ":".join(
                [
                    treatment_component(row["Kinase"]),
                    treatment_component(row["Drug"]),
                    treatment_component(row["Cancer code"]),
                    treatment_component(row["Sub-tissue"]),
                ]
            )
        )
    out["Treatment"] = treatment_values
    treatments = sorted(out["Treatment"].drop_duplicates())
    event_map = {value: f"{prefix}_{index:06d}" for index, value in enumerate(treatments, 1)}
    out["Treatment ID"] = out["Treatment"].map(event_map)
    out["Comparison baseline"] = baseline
    return out


def attach_mutation_evidence(events: pd.DataFrame, hotspot: pd.DataFrame) -> pd.DataFrame:
    evidence_columns = [
        "kinase",
        "response_cell_id",
        "AAchange",
        "mut_type",
        "AApos",
        "AAx",
        "uniprot_id",
        "Entry",
        "Hotspot",
        "Hotspot start",
        "Hotspot end",
        "Hotspot score",
        "Hotspot source",
    ]
    evidence = events.merge(
        hotspot[evidence_columns],
        left_on=["Kinase", "Cell line ID"],
        right_on=["kinase", "response_cell_id"],
        how="left",
    )
    evidence = evidence.drop_duplicates(
        KEYS + ["Cell line ID", "AAchange", "Hotspot"], keep="first"
    )
    evidence = evidence.rename(
        columns={
            "AAchange": "Mutation",
            "mut_type": "Mutation type",
            "AApos": "AA position",
            "AAx": "Reference AA",
            "uniprot_id": "Mutation UniProt ID",
            "Entry": "UniProt ID",
        }
    )
    evidence["Evidence ID"] = (
        evidence.groupby(["Treatment ID", "Cell line ID"]).cumcount().add(1)
    )
    evidence["Evidence ID"] = evidence["Treatment ID"] + "_E" + evidence["Evidence ID"].map(
        lambda value: f"{value:03d}"
    )
    return evidence


def build_primary_secondary(response: pd.DataFrame, hotspot: pd.DataFrame) -> pd.DataFrame:
    mutation_sets = (
        hotspot[hotspot["kinase"].eq("EGFR")]
        .dropna(subset=["response_cell_id", "AAchange"])
        .groupby("response_cell_id")["AAchange"]
        .agg(set)
    )
    double_cells = set(
        mutation_sets[mutation_sets.map(lambda values: {"p.L858R", "p.T790M"}.issubset(values))].index
    )
    primary_only_cells = set(
        mutation_sets[
            mutation_sets.map(lambda values: "p.L858R" in values and "p.T790M" not in values)
        ].index
    )
    if not double_cells or not primary_only_cells:
        return pd.DataFrame()

    rows = []
    egfr = response[response["Kinase"].eq("EGFR")]
    for key, group in egfr.groupby(KEYS, dropna=False):
        g1 = group[group["Cell line ID"].isin(double_cells)]
        g2 = group[group["Cell line ID"].isin(primary_only_cells)]
        if g1.empty or g2.empty:
            continue
        g2_median = g2["IC50 (uM)"].median()
        g2_records = len(g2)
        g2_cells = g2["Cell line ID"].nunique()
        for cell_id, cell_group in g1.groupby("Cell line ID"):
            g1_median = cell_group["IC50 (uM)"].median()
            fc = g1_median / g2_median
            if not np.isfinite(fc) or fc < 1.5:
                continue
            first = cell_group.iloc[0]
            row = dict(zip(KEYS, key))
            row.update(
                {
                    "Cell line ID": cell_id,
                    "Cell line name": first["Cell line name"],
                    "Cell line model ID": first["Cell line model ID"],
                    "Tissue": first["Tissue"],
                    "Source datasets": "; ".join(sorted(set(cell_group["dataset"]))),
                    "G1 response records": len(cell_group),
                    "Median(IC50, G1)": g1_median,
                    "Median(IC50, G2)": g2_median,
                    "FC(G1/G2)": fc,
                    "Primary mutation": "p.L858R",
                    "Secondary mutation": "p.T790M",
                    "Primary-only G2 response records": g2_records,
                    "Primary-only G2 cell lines": g2_cells,
                    "Comparison baseline": "EGFR L858R+T790M versus L858R-only",
                }
            )
            rows.append(row)

    if not rows:
        return pd.DataFrame()
    events = pd.DataFrame(rows)
    events["Treatment"] = events.apply(
        lambda row: ":".join(
            [
                treatment_component(row["Kinase"]),
                treatment_component(row["Drug"]),
                treatment_component(row["Cancer code"]),
                treatment_component(row["Cell line name"]),
            ]
        ),
        axis=1,
    )
    treatment_ids = {
        treatment: f"KMD2_DRK_PS_{index:06d}"
        for index, treatment in enumerate(sorted(events["Treatment"].unique()), 1)
    }
    events["Treatment ID"] = events["Treatment"].map(treatment_ids)

    secondary = hotspot[
        hotspot["kinase"].eq("EGFR") & hotspot["AAchange"].eq("p.T790M")
    ][
        [
            "response_cell_id",
            "AAchange",
            "mut_type",
            "AApos",
            "AAx",
            "Entry",
            "Hotspot",
            "Hotspot start",
            "Hotspot end",
            "Hotspot score",
            "Hotspot source",
        ]
    ].drop_duplicates(["response_cell_id", "AAchange", "Hotspot"])
    events = events.merge(secondary, left_on="Cell line ID", right_on="response_cell_id", how="left")
    events = events.rename(
        columns={
            "AAchange": "Mutation",
            "mut_type": "Mutation type",
            "AApos": "AA position",
            "AAx": "Reference AA",
            "Entry": "UniProt ID",
        }
    )
    events["Evidence ID"] = events["Treatment ID"] + "_E001"
    return events


def summarize_s1_events(evidence: pd.DataFrame) -> pd.DataFrame:
    summary = (
        evidence.groupby(["Treatment ID", "Treatment"] + KEYS, dropna=False)
        .agg(
            **{
                "Resistant cell-line events": ("Cell line ID", "nunique"),
                "Evidence rows": ("Evidence ID", "nunique"),
                "Mutations": ("Mutation", "nunique"),
                "Hotspot types": ("Hotspot", "nunique"),
                "Hotspots": ("Hotspot", lambda values: "; ".join(sorted(set(values.dropna().astype(str))))),
                "Median FC": ("FC(G1/G2)", "median"),
                "Maximum FC": ("FC(G1/G2)", "max"),
                "Median(IC50, G2)": ("Median(IC50, G2)", "first"),
                "Treatment cell lines": ("Treatment cell lines", "first"),
                "Strict-method passing cell lines": (
                    "Pass strict-method FC >= 1.5",
                    lambda values: int(values.groupby(evidence.loc[values.index, "Cell line ID"]).max().sum()),
                ),
            }
        )
        .reset_index()
    )
    return summary.sort_values(["Kinase", "Drug", "Cancer name", "Sub-tissue"]).reset_index(drop=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cancer_codes = load_cancer_code_map()
    model = load_model_metadata()
    hotspot, mutation_qc = map_mutations_to_predicted_hotspots(model)
    response, response_qc = load_response(model, cancer_codes)
    candidates = build_cell_candidates(response, hotspot)

    s1_events = candidates[candidates["Pass S1-like FC >= 1.5"]].copy()
    s1_events["Median(IC50, G2)"] = s1_events["Median(IC50, treatment background)"]
    s1_events["FC(G1/G2)"] = s1_events["FC(G1/treatment background)"]
    s1_events = add_treatment_fields(
        s1_events,
        "KMD2_DRK",
        "Treatment-background median, matching DrKinase Supplementary Table S1 output",
    )
    s1_evidence = attach_mutation_evidence(s1_events, hotspot)

    strict_events = candidates[candidates["Pass strict-method FC >= 1.5"]].copy()
    strict_events["Median(IC50, G2)"] = strict_events["Median(IC50, non-hotspot G2)"]
    strict_events["FC(G1/G2)"] = strict_events["FC(G1/non-hotspot G2)"]
    strict_events = add_treatment_fields(
        strict_events,
        "KMD2_DRK_STRICT",
        "Median of cell lines without predicted hotspot mutations",
    )
    strict_evidence = attach_mutation_evidence(strict_events, hotspot)

    primary_secondary = build_primary_secondary(response, hotspot)
    event_summary = summarize_s1_events(s1_evidence)

    evidence_columns = [
        "Treatment ID",
        "Treatment",
        "Evidence ID",
        "Kinase",
        "Drug",
        "Cancer code",
        "Cancer name",
        "Sub-tissue",
        "Tissue",
        "Cell line ID",
        "Cell line model ID",
        "Cell line name",
        "Hotspot",
        "Mutation",
        "Mutation type",
        "AA position",
        "Reference AA",
        "UniProt ID",
        "Mutation UniProt ID",
        "Hotspot start",
        "Hotspot end",
        "Hotspot score",
        "Hotspot source",
        "Median(IC50, G1)",
        "Median(IC50, G2)",
        "FC(G1/G2)",
        "Median(IC50, non-hotspot G2)",
        "FC(G1/non-hotspot G2)",
        "Pass strict-method FC >= 1.5",
        "Source datasets",
        "G1 response records",
        "Treatment response records",
        "Treatment cell lines",
        "Non-hotspot G2 response records",
        "Non-hotspot G2 cell lines",
        "Comparison baseline",
    ]
    strict_columns = [
        "Treatment ID",
        "Treatment",
        "Evidence ID",
        "Kinase",
        "Drug",
        "Cancer code",
        "Cancer name",
        "Sub-tissue",
        "Tissue",
        "Cell line ID",
        "Cell line model ID",
        "Cell line name",
        "Hotspot",
        "Mutation",
        "Mutation type",
        "AA position",
        "Reference AA",
        "UniProt ID",
        "Hotspot start",
        "Hotspot end",
        "Hotspot score",
        "Median(IC50, G1)",
        "Median(IC50, G2)",
        "FC(G1/G2)",
        "Source datasets",
        "G1 response records",
        "Non-hotspot G2 response records",
        "Non-hotspot G2 cell lines",
        "Comparison baseline",
    ]
    ps_columns = [
        "Treatment ID",
        "Treatment",
        "Evidence ID",
        "Kinase",
        "Drug",
        "Cancer code",
        "Cancer name",
        "Sub-tissue",
        "Tissue",
        "Cell line ID",
        "Cell line model ID",
        "Cell line name",
        "Primary mutation",
        "Secondary mutation",
        "Hotspot",
        "Mutation",
        "Mutation type",
        "AA position",
        "Reference AA",
        "UniProt ID",
        "Hotspot start",
        "Hotspot end",
        "Hotspot score",
        "Median(IC50, G1)",
        "Median(IC50, G2)",
        "FC(G1/G2)",
        "Source datasets",
        "G1 response records",
        "Primary-only G2 response records",
        "Primary-only G2 cell lines",
        "Comparison baseline",
    ]

    s1_evidence = s1_evidence[evidence_columns].sort_values(
        ["Kinase", "Drug", "Cancer name", "Sub-tissue", "Cell line name", "Mutation", "Hotspot"]
    )
    strict_evidence = strict_evidence[strict_columns].sort_values(
        ["Kinase", "Drug", "Cancer name", "Sub-tissue", "Cell line name", "Mutation", "Hotspot"]
    )
    if not primary_secondary.empty:
        primary_secondary = primary_secondary[ps_columns].sort_values(
            ["Drug", "Cell line name", "Mutation", "Hotspot"]
        )
    else:
        primary_secondary = pd.DataFrame(columns=ps_columns)

    qc_metrics = {
        **mutation_qc,
        **response_qc,
        "Hotspot-mutated cell-line treatment candidates": int(len(candidates)),
        "S1-like resistant cell-line events": int(
            s1_evidence[["Treatment ID", "Cell line ID"]].drop_duplicates().shape[0]
        ),
        "S1-like unique treatment groups": int(s1_evidence["Treatment ID"].nunique()),
        "S1-like evidence rows": int(len(s1_evidence)),
        "S1-like exact duplicate evidence rows": int(
            s1_evidence.duplicated(["Treatment", "Cell line ID", "Mutation", "Hotspot"]).sum()
        ),
        "Strict-method resistant cell-line events": int(
            strict_evidence[["Treatment ID", "Cell line ID"]].drop_duplicates().shape[0]
        ),
        "Strict-method unique treatment groups": int(strict_evidence["Treatment ID"].nunique()),
        "Strict-method evidence rows": int(len(strict_evidence)),
        "Primary-secondary resistant cell-line events": int(
            primary_secondary[["Treatment ID", "Cell line ID"]].drop_duplicates().shape[0]
        ),
        "Primary-secondary evidence rows": int(len(primary_secondary)),
        "Missing mutation in S1-like evidence": int(s1_evidence["Mutation"].isna().sum()),
        "Missing G2 median in S1-like evidence": int(s1_evidence["Median(IC50, G2)"].isna().sum()),
        "S1-like rows below FC threshold": int((s1_evidence["FC(G1/G2)"] < 1.5).sum()),
        "Strict-method rows below FC threshold": int((strict_evidence["FC(G1/G2)"] < 1.5).sum()),
    }

    payload = {
        "event_summary": records_for_json(event_summary),
        "s1_like_evidence": records_for_json(s1_evidence),
        "strict_method_evidence": records_for_json(strict_evidence),
        "primary_secondary": records_for_json(primary_secondary),
        "qc": [{"Metric": key, "Value": value} for key, value in qc_metrics.items()],
        "notes": [
            {
                "Section": "Primary result",
                "Description": "S1-like evidence uses GDSC1/GDSC2 only. Each row is a qualifying cell-line mutation/hotspot record within a kinase-drug-cancer-sub-tissue treatment.",
            },
            {
                "Section": "S1-like G1",
                "Description": "Median IC50 for one hotspot-mutated cell line across its available GDSC1/GDSC2 response records.",
            },
            {
                "Section": "S1-like G2",
                "Description": "Median of all drug-response records in the same treatment subgroup. This operational definition reproduces the pattern in the supplied DrKinase Supplementary Table S1.",
            },
            {
                "Section": "S1-like resistance rule",
                "Description": "Median(IC50,G1) / Median(IC50,G2) >= 1.5.",
            },
            {
                "Section": "Strict-method comparison",
                "Description": "A separate sheet follows the paper wording literally: G2 contains only cell lines without mutations in any predicted hotspot for the selected kinase.",
            },
            {
                "Section": "Primary-secondary comparison",
                "Description": "A separate sheet applies the paper's second setting to the curated EGFR example: L858R+T790M versus L858R-only. Other mutation pairs are not inferred automatically because primary/secondary status is not encoded in the raw mutation table.",
            },
            {
                "Section": "Hotspot definition",
                "Description": "Mutations were remapped directly from the original KMD2 mutation table to DrKinase-predicted Gatekeeper, G-loop, aC-Helix and A-loop intervals. Curated KMD1 substructure intervals were not mixed into this primary result.",
            },
            {
                "Section": "Treatment definition",
                "Description": "Kinase:Drug:Cancer code:Sub-tissue. Sub-tissue uses GDSC cancer_type_detail and falls back to the cancer name when detail is unavailable.",
            },
            {
                "Section": "Source: mutation",
                "Description": "all_mutation_filtered_kmd2_0420.txt",
            },
            {
                "Section": "Source: drug response",
                "Description": "all_drug_response.txt; the pre-existing label column was ignored and labels were recalculated.",
            },
            {
                "Section": "Source: hotspot prediction",
                "Description": "drkinase_results_new.csv",
            },
            {
                "Section": "Source: cell-line metadata",
                "Description": "model_metadata.csv",
            },
            {
                "Section": "Reference method",
                "Description": "Dr. Kinase, Nucleic Acids Research 2025, DOI: 10.1093/nar/gkaf366.",
            },
        ],
    }
    with (OUT_DIR / "drkinase_resistance_data.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False)
    with (OUT_DIR / "build_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(qc_metrics, handle, ensure_ascii=False, indent=2)
    print(json.dumps(qc_metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
