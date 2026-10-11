"""
Step 3: IC50 数据整合
整合 CCLE、GDSC1、GDSC2 的药物敏感数据，映射药物→激酶靶点
输出: ic50_integrated.txt
格式: kinase, dataset, cell_id, cell_name, cancer, drug, ic50, uniprot_id
      (与 V1 all_drug_response.txt 格式一致)
"""

import pandas as pd
import numpy as np
import os
import math

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR  = os.path.join(BASE_DIR, "output")

print("=" * 60)
print("Step 3: IC50 数据整合")
print("=" * 60)

# ── A. 读取激酶参考表（uniprot_id 映射）────────────────────────────────
print("\n[A] 读取激酶序列库...")
xl_kin = pd.ExcelFile(os.path.join(BASE_DIR, "kinase_base_fixed.xlsx"))
df_kin = xl_kin.parse("kinase_master", usecols=["gene_symbol","uniprot_id"])
df_kin = df_kin.drop_duplicates(subset="gene_symbol")
# gene_symbol → uniprot_id 字典
kin2uni = dict(zip(df_kin["gene_symbol"].str.upper(),
                   df_kin["uniprot_id"]))
kin_set_upper = set(kin2uni.keys())  # 快速查找用的 uppercase 集合
print(f"  激酶 gene_symbol 数: {len(kin2uni)}")

# ── B. 读取 CCLE 细胞系 ID 映射表 ────────────────────────────────────
print("\n[B] 读取 CCLE Model.csv...")
df_model = pd.read_csv(
    os.path.join(BASE_DIR, "Depmap", "Model.csv"),
    usecols=["ModelID","CellLineName","StrippedCellLineName",
             "OncotreePrimaryDisease","CCLEName"]
)
# 建立 CCLEName (如 1321N1_CENTRAL_NERVOUS_SYSTEM) → ACH-ID 的映射
ccle_id2ach = {}
for _, row in df_model.iterrows():
    if pd.notna(row["CCLEName"]):
        ccle_id2ach[row["CCLEName"]] = row["ModelID"]
# 另建 stripped_name → info
strip2info = {}
for _, row in df_model.iterrows():
    strip2info[str(row["StrippedCellLineName"]).upper()] = {
        "cell_id":   row["ModelID"],
        "cell_name": row["CellLineName"],
        "cancer":    row["OncotreePrimaryDisease"]
    }
print(f"  CCLEName 映射: {len(ccle_id2ach)} 条  StrippedName 映射: {len(strip2info)} 条")

# ── C. CCLE 药物数据 ───────────────────────────────────────────────────
print("\n[C] 处理 CCLE 药物数据...")
df_ccle = pd.read_csv(
    os.path.join(BASE_DIR, "Depmap", "CCLE_NP24.2009_Drug_data_2015.02.24.csv"),
    usecols=["CCLE Cell Line Name","Primary Cell Line Name",
             "Compound","Target","IC50 (uM)"]
)
df_ccle = df_ccle.rename(columns={
    "CCLE Cell Line Name": "ccle_full_name",
    "Primary Cell Line Name": "cell_name_raw",
    "Compound":  "drug",
    "Target":    "target_raw",
    "IC50 (uM)": "ic50"
})
# 去掉 IC50 无效值
df_ccle = df_ccle.dropna(subset=["ic50"])

# 映射 cell_id (ACH-)
df_ccle["cell_id"] = df_ccle["ccle_full_name"].map(ccle_id2ach)
# 对于没映射到的，用 stripped name 再试
def try_strip_map(row):
    if pd.notna(row["cell_id"]):
        return row["cell_id"]
    stripped = row["ccle_full_name"].split("_")[0].upper()
    info = strip2info.get(stripped, {})
    return info.get("cell_id", np.nan)

df_ccle["cell_id"] = df_ccle.apply(try_strip_map, axis=1)

# 获取 cancer / cell_name
def get_info(cell_id, field, model_df):
    row = model_df[model_df["ModelID"]==cell_id]
    if row.empty:
        return np.nan
    return row.iloc[0]["OncotreePrimaryDisease"] if field=="cancer" else row.iloc[0]["CellLineName"]

model_lookup = df_model.set_index("ModelID")
df_ccle["cancer"]    = df_ccle["cell_id"].map(lambda x: model_lookup.loc[x,"OncotreePrimaryDisease"] if pd.notna(x) and x in model_lookup.index else np.nan)
df_ccle["cell_name"] = df_ccle["cell_id"].map(lambda x: model_lookup.loc[x,"CellLineName"] if pd.notna(x) and x in model_lookup.index else np.nan)
df_ccle["dataset"]   = "CCLE"

# CCLE 靶点别名字典：非标准名 → [标准基因符号列表]
# 多靶点药物会生成多行记录
CCLE_TARGET_ALIAS = {
    "C-MET":   ["MET"],
    "MET":     ["MET"],
    "ABL":     ["ABL1"],
    "ABL1":    ["ABL1"],
    "IGF1R":   ["IGF1R"],
    "EGFR":    ["EGFR"],
    "ALK":     ["ALK"],
    "FGFR":    ["FGFR1","FGFR2","FGFR3","FGFR4"],
    "RAF":     ["RAF1","BRAF"],
    "RTK":     ["KIT","FLT3","PDGFRA","PDGFRB","VEGFR2"],  # sorafenib-like RTK
    "MEK":     ["MAP2K1","MAP2K2"],
    "CDK4":    ["CDK4"],
    "CDK6":    ["CDK6"],
    # 非激酶靶点，跳过
    "HSP90":   [],
    "MDM2":    [],
    "GS":      [],
    "TOP1":    [],
    "XIAP":    [],
    "HDAC":    [],
    "TUBB1":   [],
}

def resolve_ccle_targets(target_raw: str):
    """将 CCLE Target 字段解析为标准激酶基因符号列表"""
    result = []
    parts = [t.strip() for t in str(target_raw).replace(";",",").split(",")]
    for p in parts:
        p_up = p.upper()
        # 先查别名字典
        if p_up in CCLE_TARGET_ALIAS:
            result.extend(CCLE_TARGET_ALIAS[p_up])
        else:
            # 直接尝试匹配 kinase 库
            result.append(p)
    return list(dict.fromkeys(result))  # 去重保序

ccle_rows = []
for _, row in df_ccle.iterrows():
    targets = resolve_ccle_targets(row["target_raw"])
    for tgt in targets:
        if not tgt:  # 跳过空字符串
            continue
        tgt_up = tgt.upper()
        uni = kin2uni.get(tgt_up, np.nan)
        if pd.isna(uni):
            continue  # 非激酶，跳过
        ccle_rows.append({
            "kinase":     tgt,
            "dataset":    "CCLE",
            "cell_id":    row["cell_id"],
            "cell_name":  row["cell_name"],
            "cancer":     row["cancer"],
            "drug":       row["drug"],
            "ic50":       row["ic50"],
            "uniprot_id": uni
        })

df_ccle_out = pd.DataFrame(ccle_rows)
# 只保留能映射到激酶且有cell_id的记录
df_ccle_out = df_ccle_out.dropna(subset=["uniprot_id","cell_id"])
print(f"  CCLE 记录: {len(df_ccle_out):,}  细胞系: {df_ccle_out['cell_id'].nunique()}  "
      f"药物: {df_ccle_out['drug'].nunique()}  激酶: {df_ccle_out['kinase'].nunique()}")
print(f"  激酶列表: {sorted(df_ccle_out['kinase'].unique())}")

# ── D. GDSC 细胞系 ID 映射表 ─────────────────────────────────────────
print("\n[D] 读取 GDSC model_list...")
df_gdsc_model = pd.read_csv(
    os.path.join(BASE_DIR, "model_metadata.csv"),
    usecols=["model_id","model_name","tissue","cancer_type",
             "BROAD_ID","COSMIC_ID"]
)
# cell_line_name (model_name) → info
gdsc_name2info = {}
for _, row in df_gdsc_model.iterrows():
    gdsc_name2info[row["model_name"]] = {
        "cell_id":        row["BROAD_ID"],    # ACH-xxx
        "sanger_id":      row["model_id"],     # SIDMxxxxx
        "cancer":         row["cancer_type"],
        "cell_name":      row["model_name"],
        "COSMIC_ID":      row["COSMIC_ID"]
    }
# SANGER_MODEL_ID → info
sidm2info = {row["model_id"]: {"cell_id": row["BROAD_ID"],
                                 "cell_name": row["model_name"],
                                 "cancer": row["cancer_type"]}
             for _, row in df_gdsc_model.iterrows()}
print(f"  GDSC name 映射: {len(gdsc_name2info)} 条")

# ── E. GDSC1 + GDSC2 Excel ────────────────────────────────────────────

# GDSC PUTATIVE_TARGET 别名字典（处理非标准格式）
GDSC_TARGET_ALIAS = {
    # 大小写/格式不同
    "EEF2K":  ["EEF2K"],
    "EEF2K":  ["EEF2K"],
    "eEF2K":  ["EEF2K"],
    # 缩写映射
    "ABL":    ["ABL1"],
    "BCR-ABL":["ABL1"],
    "AKT":    ["AKT1","AKT2","AKT3"],
    "CAMK2":  ["CAMK2B"],
    "ALK5":   ["TGFBR1"],
    "DNAPK":  ["PRKDC"],
    "TAK1":   ["MAP3K7"],
    "PKCB":   ["PRKCB"],
    "PKC":    ["PRKCB","PRKCD","PRKCQ","PRKD1"],
    "FGRF1":  ["FGFR1"],  # typo in GDSC
    # 括号形式 "PDK1 (PDPK1)" → PDPK1
    "PDK1 (PDPK1)": ["PDPK1"],
    "PDK1":   ["PDK1"],
    # IKK = CHUK(IKKalpha) + IKBKB(IKKbeta)
    "IKK":    ["CHUK","IKBKB"],
    "TIE2":   ["TEK"],
    "c-FGR":  ["FGR"],
    "VEGFR2": ["KDR"],
    "CRAF":   ["RAF1"],
    "ARAF, BRAF, CRAF": ["ARAF","BRAF","RAF1"],
}

def parse_gdsc_targets(target_raw: str, kin_set: set):
    """解析 GDSC PUTATIVE_TARGET 字段，返回有效激酶基因符号列表"""
    # 先整体检查别名字典
    stripped = str(target_raw).strip()
    if stripped in GDSC_TARGET_ALIAS:
        return [k for k in GDSC_TARGET_ALIAS[stripped] if k.upper() in kin_set]

    result = []
    # 分割（逗号或分号）
    parts = [p.strip() for p in stripped.replace(";",",").split(",")]
    for p in parts:
        p_stripped = p.strip()
        # 检查别名字典
        if p_stripped in GDSC_TARGET_ALIAS:
            result.extend([k for k in GDSC_TARGET_ALIAS[p_stripped]
                           if k.upper() in kin_set])
        elif p_stripped.upper() in kin_set:
            result.append(p_stripped)
        else:
            # 去括号内容再试 e.g. "PDK1 (PDPK1)" → "PDPK1"
            import re
            m = re.search(r'\(([^)]+)\)', p_stripped)
            if m:
                inner = m.group(1).strip()
                if inner.upper() in kin_set:
                    result.append(inner)
    return list(dict.fromkeys(result))  # 去重

def process_gdsc_excel(filepath, dataset_label):
    print(f"\n  处理 {dataset_label} ({os.path.basename(filepath)})...")
    xl = pd.ExcelFile(filepath)
    df = xl.parse(xl.sheet_names[0],
                  usecols=["CELL_LINE_NAME","SANGER_MODEL_ID",
                           "CANCER_TYPE","DRUG_NAME","PUTATIVE_TARGET","LN_IC50"])

    df = df.dropna(subset=["LN_IC50"])
    # IC50(μM) = exp(LN_IC50)
    df["ic50"] = np.exp(df["LN_IC50"])

    rows = []
    for _, row in df.iterrows():
        # cell line ID → ACH
        sidm = row["SANGER_MODEL_ID"]
        cname = row["CELL_LINE_NAME"]
        info = sidm2info.get(sidm, gdsc_name2info.get(cname, {}))
        cell_id  = info.get("cell_id", np.nan)
        cancer   = info.get("cancer", row.get("CANCER_TYPE", np.nan))
        cn       = info.get("cell_name", cname)

        # 多靶点拆分（使用带别名字典的解析函数）
        targets = parse_gdsc_targets(str(row["PUTATIVE_TARGET"]), kin_set_upper)
        for tgt in targets:
            tgt_up = tgt.upper()
            uni = kin2uni.get(tgt_up, np.nan)
            if pd.isna(uni):
                continue  # 仍不在激酶库中，跳过
            rows.append({
                "kinase":     tgt,
                "dataset":    dataset_label,
                "cell_id":    cell_id,
                "cell_name":  cn,
                "cancer":     cancer,
                "drug":       row["DRUG_NAME"],
                "ic50":       row["ic50"],
                "uniprot_id": uni
            })

    df_out = pd.DataFrame(rows)
    df_out = df_out.dropna(subset=["cell_id"])
    print(f"    记录: {len(df_out):,}  细胞系: {df_out['cell_id'].nunique()}  "
          f"药物: {df_out['drug'].nunique()}  激酶: {df_out['kinase'].nunique()}")
    return df_out

print("\n[E] 处理 GDSC1 和 GDSC2...")
df_gdsc1_out = process_gdsc_excel(
    os.path.join(BASE_DIR, "GDSC", "GDSC1_fitted_dose_response_27Oct23.xlsx"),
    "GDSC1"
)
df_gdsc2_out = process_gdsc_excel(
    os.path.join(BASE_DIR, "GDSC", "GDSC2_fitted_dose_response_27Oct23.xlsx"),
    "GDSC2"
)

# ── F. 合并所有数据集 ─────────────────────────────────────────────────
print("\n[F] 合并所有数据集...")
df_all = pd.concat([df_ccle_out, df_gdsc1_out, df_gdsc2_out], ignore_index=True)

# 确保列顺序与 V1 一致
col_order = ["kinase","dataset","cell_id","cell_name","cancer","drug","ic50","uniprot_id"]
df_all = df_all[[c for c in col_order if c in df_all.columns]]

# 去重（同一 dataset/cell_id/kinase/drug 保留一条，按 IC50 无偏好取第一条）
df_all = df_all.drop_duplicates(subset=["dataset","cell_id","kinase","drug"])

out_path = os.path.join(OUT_DIR, "ic50_integrated.txt")
df_all.to_csv(out_path, sep="\t", index=False)

print(f"\n  ✅ 保存: {out_path}")
print(f"\n  ━━━ 总体统计 ━━━")
print(f"  总记录数: {len(df_all):,}")
print(f"  各数据集:")
print(df_all["dataset"].value_counts().to_string())
print(f"\n  唯一细胞系 (cell_id): {df_all['cell_id'].nunique():,}")
print(f"  唯一药物:   {df_all['drug'].nunique():,}")
print(f"  唯一激酶:   {df_all['kinase'].nunique():,}")
print(f"  IC50 范围:  {df_all['ic50'].min():.4f} ~ {df_all['ic50'].max():.2f} μM")

print("\n" + "=" * 60)
print("Step 3 完成！")
print("=" * 60)
