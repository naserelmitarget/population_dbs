#!/usr/bin/env python3
"""Flag genes and variants whose population frequencies should not be trusted; validate against GIAB regions.

usage: unreliable_regions.py <extra_analyses dir>   (called by unreliable_regions.sh)

Gene rules (exonic, protein-coding):
  R1  high for common (AF >= 1%) protein-altering variants vs gnomAD in >= 5 datasets across >= 2 population groups
  R2  high for common synonymous variants vs gnomAD in >= 5 datasets (a technical signal: synonymous variants are
      not under the same selection, so excess points to calling / mapping problems)
  R3  within-dataset protein-altering density outlier (robust z >= 3) in at least half of the datasets
Variant list: ClinVar P/LP, VUS and conflicting variants whose high frequency (>= 5%) rests on one dataset only.
Validation: share of each gene's exonic footprint inside GIAB segmental duplications / low-mappability regions,
flagged genes against all other protein-coding genes (Mann-Whitney U).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

W = sys.argv[1]
K = os.environ.get("K", "/home/NASER/annotation/population_study/data/archive")
O = os.path.join(W, "unreliable_regions")
T, F, TMP = (os.path.join(O, x) for x in ("tables", "figures", "tmp"))

u = pd.read_csv(f"{K}/gene_level/tables/C_high_genes_vs_gnomad_by_uniqueness.tsv", sep="\t")
d = pd.read_csv(f"{K}/gene_level/tables/C_high_density_genes_by_uniqueness.tsv", sep="\t")
r1 = u[(u["class"] == "protein_altering") & (u.af == "common") & (u.n_datasets >= 5) & u.groups.str.contains(",")]
r2 = u[(u["class"] == "synonymous") & (u.af == "common") & (u.n_datasets >= 5)]
r3 = d[(d["class"] == "protein_altering") & (d.af == "all") & (d.n_datasets >= 0.5 * d.of_datasets)]
genes = pd.DataFrame(index=sorted(set(r1.gene) | set(r2.gene) | set(r3.gene)))
genes["R1_pa_high_vs_gnomad_datasets"] = r1.set_index("gene").n_datasets
genes["R2_syn_high_vs_gnomad_datasets"] = r2.set_index("gene").n_datasets
genes["R3_density_outlier_datasets"] = r3.set_index("gene").n_datasets
sym = pd.concat([u[["gene", "symbol"]], d[["gene", "symbol"]]]).drop_duplicates("gene").set_index("gene").symbol
genes["symbol"] = genes.index.map(sym)
genes["n_rules"] = genes[[c for c in genes.columns if c.startswith("R")]].notna().sum(axis=1)

# GIAB overlap per gene (key "ENSG|SYMBOL")
cov = None
for b in ("segdups", "lowmappabilityall"):
    x = pd.read_csv(f"{TMP}/gene_{b}.tsv", sep="\t", header=None, names=["key", "positions", f"in_{b}"])
    cov = x if cov is None else cov.merge(x[["key", f"in_{b}"]], on="key")
cov["gene"] = cov.key.str.split("|").str[0]
cov["symbol"] = cov.key.str.split("|").str[1]
cov["frac_segdups"] = cov.in_segdups / cov.positions
cov["frac_lowmap"] = cov.in_lowmappabilityall / cov.positions
cov = cov.groupby("gene").agg(symbol=("symbol", "first"), positions=("positions", "sum"),
                              frac_segdups=("frac_segdups", "mean"), frac_lowmap=("frac_lowmap", "mean"))
genes = genes.join(cov[["positions", "frac_segdups", "frac_lowmap"]], how="left")
genes["symbol"] = genes.symbol.fillna(cov.symbol)
genes["tier"] = np.where(genes.n_rules >= 2, "high confidence (>= 2 rules)", "one rule")
genes = genes.sort_values(["n_rules", "R1_pa_high_vs_gnomad_datasets"], ascending=[False, False])
genes.rename_axis("gene").to_csv(os.path.join(T, "unreliable_genes.tsv"), sep="\t", float_format="%.3f")
span = pd.read_csv(f"{TMP}/gene_span.bed", sep="\t", header=None, names=["chr", "start", "end", "key"])
span["gene"] = span.key.str.split("|").str[0]
bed = span[span.gene.isin(genes.index)].copy()
bed["name"] = bed.key.str.split("|").str[1] + "|" + bed.gene.map(genes.tier).str.split(" ").str[0]
bed.assign(chr="chr" + bed.chr.astype(str))[["chr", "start", "end", "name"]].to_csv(
    os.path.join(T, "unreliable_genes.hg38.bed"), sep="\t", index=False, header=False)

# validation
rows = []
flag = cov.index.isin(genes.index)
for col, lab in (("frac_segdups", "segmental duplications"), ("frac_lowmap", "low mappability")):
    a, b = cov.loc[flag, col].dropna(), cov.loc[~flag, col].dropna()
    rows.append({"GIAB region": lab, "flagged genes": len(a), "median share (flagged)": a.median(),
                 "mean share (flagged)": a.mean(), "flagged with > 50% inside": (a > 0.5).mean(),
                 "other genes": len(b), "median share (other)": b.median(), "mean share (other)": b.mean(),
                 "other with > 50% inside": (b > 0.5).mean(), "Mann-Whitney p": mannwhitneyu(a, b, alternative="greater").pvalue})
pd.DataFrame(rows).to_csv(os.path.join(T, "giab_enrichment.tsv"), sep="\t", index=False, float_format="%.4g")

# variant list: ClinVar frequencies resting on one dataset
vl = []
for reg in ("exonic", "non_exonic"):
    x = pd.read_csv(f"{K}/study/{reg}/tables/step4_clinvar_PLP_VUS_common_somewhere_rare_in_gnomad.tsv", sep="\t", dtype={"chr": str})
    x = x[x.support <= 1].assign(region=reg)
    vl.append(x)
vl = pd.concat(vl)
vl[["region", "chr", "start", "end", "ref", "alt", "gene", "clnsig_group", "clnsig", "review_stars", "max_af",
    "datasets_af_ge_5pct", "gnomad_popmax", "support"]].to_csv(
    os.path.join(T, "unreliable_variants_clinvar_single_dataset.tsv"), sep="\t", index=False)
sd = vl.datasets_af_ge_5pct.str.split(";").explode().str.split(":").str[0].value_counts()
sd.rename_axis("dataset").rename("single-dataset ClinVar frequency flags").to_csv(
    os.path.join(T, "unreliable_variants_by_dataset.tsv"), sep="\t")

# figures
fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
for a_, col, lab in ((ax[0], "frac_segdups", "segmental duplications"), (ax[1], "frac_lowmap", "low-mappability regions")):
    data = [cov.loc[~flag, col].dropna() * 100, cov.loc[flag, col].dropna() * 100]
    a_.boxplot(data, showfliers=False, widths=0.5)
    a_.set_xticks([1, 2])
    a_.set_xticklabels([f"other genes (n={len(data[0]):,})", f"flagged genes (n={len(data[1]):,})"])
    a_.set_ylabel(f"% of exonic footprint in GIAB {lab}")
    a_.set_title(lab)
fig.suptitle("Genes flagged as unreliable lie in GIAB difficult regions far more often")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(F, f"giab_overlap_flagged_vs_other.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
fig, ax = plt.subplots(figsize=(10, 0.4 * len(sd) + 2), constrained_layout=True)
s = sd.iloc[::-1]
ax.barh(s.index, s.values, color="#d62728")
ax.set_xlabel("ClinVar P/LP, VUS or conflicting variants whose >= 5% frequency comes from this dataset alone")
ax.set_title("Single-dataset frequency flags per dataset (probable data errors)")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(F, f"single_dataset_flags_per_dataset.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"unreliable genes: {len(genes)} ({int((genes.n_rules >= 2).sum())} with >= 2 rules); "
      f"single-dataset ClinVar variants: {len(vl)}")
