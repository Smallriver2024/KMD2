"""
Step 1: 激酶结构域突变映射
将 all_mutation_filtered_kmd2_0420.txt 中的突变按UniProt坐标映射到
output/kinase_domain_annotation_merged.txt 定义的四个结构域
输出: domain_mutations.txt
"""

import pandas as pd
import numpy as np
import os
import time

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR  = os.path.join(BASE_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("Step 1: 结构域突变映射")
print("=" * 60)

# ── 1. 读取突变数据（只保留 CCLE 和 GDSC，其余无药物数据）──────────
print("\n[1/4] 读取突变文件...")
t0 = time.time()
df_mut = pd.read_csv(
    os.path.join(BASE_DIR, "all_mutation_filtered_kmd2_0420.txt"),
    sep="\t",
    dtype={"AApos": float}
)
print(f"  原始记录: {len(df_mut):,}  用时: {time.time()-t0:.1f}s")

# 只保留 CCLE 和 GDSC（有药物敏感数据的数据集）
df_mut = df_mut[df_mut["dataset"].isin(["CCLE", "GDSC"])].copy()
print(f"  CCLE+GDSC 记录: {len(df_mut):,}")

# 处理 uniprot_id 亚型后缀：P21127-2 → P21127
df_mut["uniprot_main"] = df_mut["uniprot_id"].str.split("-").str[0]

# 去掉 AApos 为空的行
df_mut = df_mut.dropna(subset=["AApos"])
df_mut["AApos"] = df_mut["AApos"].astype(int)
print(f"  去除 AApos 空值后: {len(df_mut):,}")

# ── 2. 读取结构域数据 ────────────────────────────────────────────────
print("\n[2/4] 读取结构域文件...")
domain_path = os.path.join(OUT_DIR, "kinase_domain_annotation_merged.txt")
if not os.path.exists(domain_path):
    raise FileNotFoundError(
        f"未找到合并结构域注释表: {domain_path}\n"
        "请先运行: python step0_merge_domain_annotations.py"
    )

df_dom = pd.read_csv(domain_path, sep="\t")
print(f"  结构域记录: {len(df_dom):,}  UniProt条目: {df_dom['Entry'].nunique()}")
print(f"  结构域分布:\n{df_dom['domain_type'].value_counts().to_string()}")
print(f"  来源分布:\n{df_dom['source'].value_counts().to_string()}")

# ── 3. 区间包含判断（向量化）─────────────────────────────────────────
print("\n[3/4] 执行区间映射...")

results = []
for entry, grp_dom in df_dom.groupby("Entry"):
    # 取该UniProt对应的所有突变
    sub_mut = df_mut[df_mut["uniprot_main"] == entry]
    if sub_mut.empty:
        continue
    # 对每个突变位点，检查是否落在任一结构域区间内
    for _, dom_row in grp_dom.iterrows():
        start, end = int(dom_row["domain_start"]), int(dom_row["domain_end"])
        dr_type = dom_row["domain_type"]
        score = dom_row["domain_score"]
        # 向量化判断
        hit = sub_mut[(sub_mut["AApos"] >= start) & (sub_mut["AApos"] <= end)].copy()
        if hit.empty:
            continue
        hit["domain_type"]  = dr_type
        hit["domain_start"] = start
        hit["domain_end"]   = end
        hit["domain_score"] = score
        hit["domain_source"] = dom_row["source"]
        hit["domain_source_priority"] = dom_row["source_priority"]
        results.append(hit)

if results:
    df_domain_mut = pd.concat(results, ignore_index=True)
else:
    df_domain_mut = pd.DataFrame()

print(f"  命中结构域的突变记录: {len(df_domain_mut):,}")

# ── 4. 整理输出列 ─────────────────────────────────────────────────────
print("\n[4/4] 整理并保存结果...")

if not df_domain_mut.empty:
    # 列排序
    out_cols = [
        "dataset", "kinase", "gene_id", "hgnc_id", "cancer",
        "sample_id", "AAchange", "mut_type", "AApos", "AAx",
        "uniprot_id", "uniprot_main",
        "domain_type", "domain_start", "domain_end", "domain_score",
        "domain_source", "domain_source_priority"
    ]
    df_out = df_domain_mut[[c for c in out_cols if c in df_domain_mut.columns]]

    # 去重：同一突变命中多个同类型结构域时，先保留已有注释，再保留最高 score
    df_out = (df_out.sort_values(
                        ["domain_source_priority", "domain_score"],
                        ascending=[True, False]
                    )
                    .drop_duplicates(
                        subset=["dataset", "sample_id", "kinase", "AApos", "domain_type"],
                        keep="first"
                    ))

    out_path = os.path.join(OUT_DIR, "domain_mutations.txt")
    df_out.to_csv(out_path, sep="\t", index=False)
    print(f"  ✅ 保存: {out_path}")
    print(f"  记录数: {len(df_out):,}")
    print(f"\n  结构域分布:")
    print(df_out["domain_type"].value_counts().to_string())
    print(f"\n  数据集分布:")
    print(df_out["dataset"].value_counts().to_string())
    print(f"\n  涉及细胞系数: {df_out['sample_id'].nunique():,}")
    print(f"  涉及激酶数: {df_out['kinase'].nunique():,}")
else:
    print("  ⚠️ 未找到任何结构域突变，请检查数据！")

print("\n" + "=" * 60)
print("Step 1 完成！")
print("=" * 60)
