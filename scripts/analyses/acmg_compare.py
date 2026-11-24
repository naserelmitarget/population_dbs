#!/usr/bin/env python3
"""Compare BIAS-2015 classifications: A (gnomAD only) vs B1 / B2 (population datasets added), per sample and pooled.

usage: acmg_compare.py <acmg_impact dir> <sample> [<sample> ...]
reads  <dir>/<sample>/bias_{A,B1,B2}.tsv
writes <dir>/tables/{summary.tsv, class_transitions_A_to_B{1,2}[_<sample>].tsv, variants_with_class_change.tsv},
       <dir>/figures/{class_counts_by_scenario, vus_resolved_per_sample}.{png,pdf}
Frequency codes: BA1, BS1 (any strength), PM2 (any strength).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D, SAMPLES = sys.argv[1], sys.argv[2:]
T, F = os.path.join(D, "tables"), os.path.join(D, "figures")
os.makedirs(T, exist_ok=True)
os.makedirs(F, exist_ok=True)
ORDER = ["pathogenic", "likely pathogenic", "uncertain", "likely benign", "benign"]   # BIAS-2015 labels
PLP, LBB = ORDER[:2], ORDER[3:]


def load(sample, name):
    x = pd.read_csv(os.path.join(D, sample, f"bias_{name}.tsv"), sep="\t", dtype=str, low_memory=False)
    x["key"] = x[list(x.columns[:4])].agg(":".join, axis=1)
    codes = x.acmgCodes.fillna("").str.split(",")
    x["freq_codes"] = codes.map(lambda c: ",".join(sorted(k for k in c if k.split("_")[0] in ("BA1", "BS1", "PM2"))))
    return x.drop_duplicates("key").set_index("key")


def summarise(m, b):
    a, n = m.acmgClassification, m[f"acmgClassification_{b}"]
    fa, fb = m.freq_codes, m[f"freq_codes_{b}"]
    return {"variants": len(m), "class changed": int((a != n).sum()), "frequency codes changed": int((fa != fb).sum()),
            "PM2 lost": int((fa.str.contains("PM2") & ~fb.str.contains("PM2")).sum()),
            "BS1 gained": int((~fa.str.contains("BS1") & fb.str.contains("BS1")).sum()),
            "BA1 gained": int((~fa.str.contains("BA1") & fb.str.contains("BA1")).sum()),
            "P/LP in A": int(a.isin(PLP).sum()), "P/LP downgraded": int((a.isin(PLP) & ~n.isin(PLP)).sum()),
            "towards pathogenic": int((~a.isin(PLP) & n.isin(PLP)).sum()),
            "VUS in A": int((a == "uncertain").sum()), "VUS -> LB/B": int(((a == "uncertain") & n.isin(LBB)).sum()),
            "LB -> B": int(((a == "likely benign") & (n == "benign")).sum())}


rows, changed, merged, counts = [], [], {"B1": [], "B2": []}, []
for s in SAMPLES:
    A = load(s, "A")
    counts.append(A.acmgClassification.value_counts().rename((s, "A")))
    for b in ("B1", "B2"):
        B = load(s, b)
        counts.append(B.acmgClassification.value_counts().rename((s, b)))
        cols = ["acmgClassification", "acmgCodes", "freq_codes"]
        gene = [c for c in A.columns if c.lower() in ("genename", "gene", "symbol")][:1]
        m = A[cols + gene].join(B[cols], rsuffix=f"_{b}", how="inner")
        merged[b].append(m)
        rows.append({"sample": s, "scenario": b, **summarise(m, b)})
        pd.crosstab(m.acmgClassification, m[f"acmgClassification_{b}"]).reindex(index=ORDER, columns=ORDER, fill_value=0) \
            .to_csv(os.path.join(T, f"class_transitions_A_to_{b}_{s}.tsv"), sep="\t")
        ch = m[m.acmgClassification != m[f"acmgClassification_{b}"]]
        changed.append(ch.assign(sample=s, scenario=b).rename(columns={f"acmgClassification_{b}": "class_new",
                                                                        f"acmgCodes_{b}": "codes_new",
                                                                        f"freq_codes_{b}": "freq_codes_new"}))
for b in ("B1", "B2"):
    m = pd.concat(merged[b])
    rows.append({"sample": "all samples", "scenario": b, **summarise(m, b)})
    pd.crosstab(m.acmgClassification, m[f"acmgClassification_{b}"]).reindex(index=ORDER, columns=ORDER, fill_value=0) \
        .to_csv(os.path.join(T, f"class_transitions_A_to_{b}.tsv"), sep="\t")
S = pd.DataFrame(rows)
S["% VUS resolved"] = (100 * S["VUS -> LB/B"] / S["VUS in A"]).round(1)
S.to_csv(os.path.join(T, "summary.tsv"), sep="\t", index=False)
pd.concat(changed).rename_axis("variant").reset_index().to_csv(os.path.join(T, "variants_with_class_change.tsv"), sep="\t", index=False)

C = pd.concat(counts, axis=1).reindex(ORDER).fillna(0)
pooled = C.T.groupby(level=1).sum().T
fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
xs = np.arange(len(ORDER))
labels = {"A": "A: gnomAD only", "B1": "B1: + highest population AF", "B2": "B2: + AF agreed by 2 datasets"}
for i, (sc, col) in enumerate(zip(["A", "B1", "B2"], ["#7f7f7f", "#1f77b4", "#2ca02c"])):
    ax.bar(xs + (i - 1) * 0.27, pooled[sc], width=0.27, color=col, label=labels[sc])
ax.set_xticks(xs)
ax.set_xticklabels(ORDER)
ax.set_yscale("log")
ax.set_ylabel("variants (log scale)")
ax.set_title(f"ACMG classes of exonic variants in {len(SAMPLES)} genome(s), BIAS-2015, with and without population datasets")
ax.legend()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(F, f"class_counts_by_scenario.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
per = S[S["sample"] != "all samples"].pivot(index="sample", columns="scenario", values="% VUS resolved")
per.plot.bar(ax=ax, color=["#1f77b4", "#2ca02c"], rot=0)
ax.set_ylabel("% of VUS moved to likely benign / benign")
ax.set_title("VUS resolved by adding population datasets, per genome")
ax.legend(["B1: highest population AF", "B2: two datasets agree"])
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(F, f"vus_resolved_per_sample.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
print(S.to_string(index=False))
