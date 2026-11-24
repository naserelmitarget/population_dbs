#!/usr/bin/env python3
"""Gene-level variant burden, exonic regions (gene_burden.sh output).

usage: gene_burden_report.py <gene_burden dir>

A. Across populations, gnomAD v4.1 as reference. For gene g in dataset d (per variant class and AF class):
     share_d(g) = n_d(g) / N_d        share_gnomAD(g) = n_gn(g) / N_gn
     ratio = share_d / share_gnomAD;  expected = N_d * share_gnomAD;  two-sided binomial test, BH-FDR per dataset.
   Gene length cancels (same gene in both); dataset size cancels through the shares.
   Called "high" / "low" when FDR < 0.05, |log2 ratio| >= 1 and max(observed, expected) >= 5.
B. Within a dataset: variant density per kb of the gene's exonic footprint, robust z-score of log10 density
   (median / MAD over genes with >= 300 bp footprint). Outlier when |z| >= 3.
C. Uniqueness: in how many datasets (and population groups) each gene is high in A / outlier in B.
Protein-coding genes only.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

D = sys.argv[1]
T, F = os.path.join(D, "tables"), os.path.join(D, "figures")
os.makedirs(T, exist_ok=True)
os.makedirs(F, exist_ok=True)

GROUP = {"tommo": "East Asia", "kova": "East Asia", "wbbc": "East Asia", "sg10k_npm": "East Asia",
         "mcps": "Latin America", "finngen": "Europe", "swegen": "Europe", "gonl": "Europe", "ukbb": "Europe",
         "saudi": "Middle East", "gme": "Middle East"}
REF = "gnomad41"
MINLEN = 300


def save(fig, name):
    fig.savefig(os.path.join(F, name + ".png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(F, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


gc = os.path.join(D, "gene_counts.tsv")
c = pd.read_csv(gc if os.path.exists(gc) else gc + ".gz", sep="\t", dtype={"symbol": str})
c = c[c.biotype == "protein_coding"]
L = pd.read_csv(os.path.join(D, "gene_length.tsv"), sep="\t", dtype={"symbol": str})
L = L.groupby("gene", as_index=False).exonic_positions.max()
sym = c.drop_duplicates("gene").set_index("gene").symbol

# ---------- A: each dataset vs gnomAD ----------
rows = []
for (cl, af), g in c.groupby(["class", "af"]):
    w = g.pivot_table(index="gene", columns="dataset", values="variants", aggfunc="sum", fill_value=0)
    if REF not in w:
        continue
    ref = w[REF]
    Ngn = ref.sum()
    for d in w.columns:
        if d == REF:
            continue
        obs, Nd = w[d], w[d].sum()
        keep = (ref > 0) | (obs > 0)
        sh = (ref[keep] + 0.5) / (Ngn + 0.5 * keep.sum())           # smoothed reference share
        exp = Nd * sh
        o = obs[keep]
        p = [binomtest(int(k), int(Nd), float(s)).pvalue if (k >= 5 or e >= 5) else 1.0
             for k, e, s in zip(o.values, exp.values, sh.values)]
        r = pd.DataFrame({"gene": o.index, "dataset": d, "class": cl, "af": af, "observed": o.values,
                          "expected": exp.values, "gnomad": ref[keep].values,
                          "log2_ratio": np.log2((o.values + 0.5) / (exp.values + 0.5)), "p": p})
        r["fdr"] = bh(r.p)
        rows.append(r)
A = pd.concat(rows, ignore_index=True)
A["symbol"] = A.gene.map(sym)
A["call"] = np.where((A.fdr < 0.05) & (A.log2_ratio >= 1) & (np.maximum(A.observed, A.expected) >= 5), "high",
                     np.where((A.fdr < 0.05) & (A.log2_ratio <= -1) & (np.maximum(A.observed, A.expected) >= 5), "low", ""))
A["group"] = A.dataset.map(GROUP).fillna("multi-population / other")
A.to_csv(os.path.join(T, "A_gene_vs_gnomad_all.tsv.gz"), sep="\t", index=False, float_format="%.4g", compression="gzip")
A[A.call != ""].sort_values(["class", "af", "dataset", "log2_ratio"], ascending=[True, True, True, False]).to_csv(
    os.path.join(T, "A_gene_vs_gnomad_significant.tsv"), sep="\t", index=False, float_format="%.4g")

# uniqueness for A: number of datasets / groups where the gene is high (or low)
for call in ["high", "low"]:
    s = A[A.call == call]
    u = s.groupby(["class", "af", "gene"]).agg(symbol=("symbol", "first"), n_datasets=("dataset", "nunique"),
                                                datasets=("dataset", lambda x: ",".join(sorted(x))),
                                                groups=("group", lambda x: ",".join(sorted(set(x)))),
                                                max_log2_ratio=("log2_ratio", "max" if call == "high" else "min")).reset_index()
    u["pattern"] = np.where(u.n_datasets == 1, "unique to one dataset",
                            np.where(~u.groups.str.contains(","), "one population group", "several groups"))
    u.sort_values(["class", "af", "n_datasets"], ascending=[True, True, False]).to_csv(
        os.path.join(T, f"C_{call}_genes_vs_gnomad_by_uniqueness.tsv"), sep="\t", index=False, float_format="%.3g")

# ---------- B: within-dataset density outliers ----------
B = c.merge(L, on="gene", how="inner")
B = B[B.exonic_positions >= MINLEN].copy()
B["per_kb"] = B.variants / (B.exonic_positions / 1000)
out = []
for (d, cl, af), g in B.groupby(["dataset", "class", "af"]):
    x = np.log10(g.per_kb + 0.01)
    med, mad = x.median(), (x - x.median()).abs().median() * 1.4826
    g = g.assign(z=(x - med) / mad if mad > 0 else 0.0)
    out.append(g)
B = pd.concat(out, ignore_index=True)
B["outlier"] = np.where(B.z >= 3, "high", np.where(B.z <= -3, "low", ""))
B.to_csv(os.path.join(T, "B_gene_density_within_dataset_all.tsv.gz"), sep="\t", index=False, float_format="%.4g", compression="gzip")
Bo = B[B.outlier != ""]
Bo.sort_values(["dataset", "class", "af", "z"], ascending=[True, True, True, False]).to_csv(
    os.path.join(T, "B_gene_density_outliers.tsv"), sep="\t", index=False, float_format="%.4g")
nds = B.groupby(["class", "af"]).dataset.nunique()
ub = (Bo[Bo.outlier == "high"].groupby(["class", "af", "gene"])
      .agg(symbol=("symbol", "first"), n_datasets=("dataset", "nunique"), datasets=("dataset", lambda x: ",".join(sorted(x))),
           median_z=("z", "median")).reset_index())
ub["of_datasets"] = [nds.get((a, b), 0) for a, b in zip(ub["class"], ub["af"])]
ub["pattern"] = np.where(ub.n_datasets == 1, "unique to one dataset",
                         np.where(ub.n_datasets >= 0.5 * ub.of_datasets, "shared (most datasets)", "some datasets"))
ub.sort_values(["class", "af", "n_datasets"], ascending=[True, True, False]).to_csv(
    os.path.join(T, "C_high_density_genes_by_uniqueness.tsv"), sep="\t", index=False, float_format="%.3g")

# ---------- figures ----------
# 1. number of significant genes per dataset (A), protein-altering, all and common
for af in ["all", "common"]:
    s = A[(A["class"] == "protein_altering") & (A.af == af) & (A.call != "")]
    if s.empty:
        continue
    hi = A[(A["class"] == "protein_altering") & (A.af == af) & (A.call == "high")]
    uniq = hi.groupby("gene").dataset.nunique()
    hi = hi.assign(kind=np.where(hi.gene.map(uniq) == 1, "high, only this dataset", "high, also elsewhere"))
    tab = pd.concat([hi.groupby(["dataset", "kind"]).size().unstack(fill_value=0),
                     s[s.call == "low"].groupby("dataset").size().rename("low")], axis=1).fillna(0)
    tab = tab.loc[tab.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(13, 0.42 * len(tab) + 2.5), constrained_layout=True)
    left = np.zeros(len(tab))
    for k, col in [("high, only this dataset", "#d62728"), ("high, also elsewhere", "#ff9f40"), ("low", "#1f77b4")]:
        if k in tab:
            ax.barh(tab.index, tab[k], left=left, color=col, label=k)
            left += tab[k].values
    ax.set_xlabel("protein-coding genes with significantly more / fewer variants than gnomAD predicts (FDR < 5%, >= 2-fold)")
    ax.set_title(f"Protein-altering variants ({'all' if af == 'all' else 'AF >= 1%'}): genes deviating from gnomAD v4.1, per dataset")
    ax.legend()
    save(fig, f"A_significant_genes_per_dataset_protein_altering_{af}")

# 2. heatmap: log2 ratio of the most often high genes, protein-altering common, country datasets
for af in ["common", "all"]:
    s = A[(A["class"] == "protein_altering") & (A.af == af)]
    top = s[s.call == "high"].groupby("gene").size().sort_values(ascending=False).head(40).index
    if len(top) == 0:
        continue
    cols = [d for d in list(GROUP) if d in s.dataset.unique()] + sorted(set(s.dataset) - set(GROUP))
    h = s[s.gene.isin(top)].pivot_table(index="gene", columns="dataset", values="log2_ratio").reindex(index=top, columns=cols)
    sig = s[s.gene.isin(top)].pivot_table(index="gene", columns="dataset", values="call", aggfunc="first").reindex(index=top, columns=cols)
    fig, ax = plt.subplots(figsize=(0.55 * len(cols) + 5, 0.32 * len(top) + 3), constrained_layout=True)
    im = ax.imshow(h.values, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            if sig.iat[i, j] in ("high", "low"):
                ax.text(j, i, "*", ha="center", va="center", fontsize=9)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=9)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([f"{sym.get(g, g)}" for g in top], fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.03, label="log2 (observed / expected from gnomAD)")
    ax.set_title(f"Genes most often high vs gnomAD, protein-altering ({'AF >= 1%' if af == 'common' else 'all'}); * = significant")
    save(fig, f"A_heatmap_top_genes_protein_altering_{af}")

# 3. within-dataset density: top high-density genes per country dataset (protein-altering, all)
s = B[(B["class"] == "protein_altering") & (B.af == "all") & B.dataset.isin(list(GROUP) + [REF])]
if not s.empty:
    top = (s[s.outlier == "high"].groupby("gene").size().sort_values(ascending=False).head(40).index)
    cols = [d for d in [REF] + list(GROUP) if d in s.dataset.unique()]
    h = s[s.gene.isin(top)].pivot_table(index="gene", columns="dataset", values="z").reindex(index=top, columns=cols)
    fig, ax = plt.subplots(figsize=(0.6 * len(cols) + 5, 0.32 * len(top) + 3), constrained_layout=True)
    im = ax.imshow(h.values, cmap="magma_r", vmin=0, vmax=8, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90, fontsize=9)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([sym.get(g, g) for g in top], fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.03, label="robust z of variants per kb (within dataset)")
    ax.set_title("Highest protein-altering variant density within each dataset (gnomAD + country datasets)")
    save(fig, "B_heatmap_high_density_genes_protein_altering")
print("GENE_REPORT_DONE")
