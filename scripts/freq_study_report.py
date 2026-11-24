#!/usr/bin/env python3
"""Tables and figures for the frequency study (freq_study.sh + freq_clinvar.sh outputs).

usage: freq_study_report.py <freq_study region dir> [region label]
writes <dir>/tables/*.tsv and <dir>/figures/*.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = sys.argv[1]
REG = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(os.path.normpath(D))
REGT = {"exonic": "exonic regions (gnomAD v4 exome calling intervals)",
        "non_exonic": "non-exonic regions (gnomAD genomes only)"}.get(REG, REG)
T, F = os.path.join(D, "tables"), os.path.join(D, "figures")
os.makedirs(T, exist_ok=True)
os.makedirs(F, exist_ok=True)

COUNTRY = {"tommo": "Japan (ToMMo)", "kova": "Korea (KOVA)", "wbbc": "China (WBBC)", "sg10k_npm": "Singapore (SG10K)",
           "mcps": "Mexico (MCPS)", "finngen": "Finland (FinnGen)", "swegen": "Sweden (SweGen)",
           "gonl": "Netherlands (GoNL)", "ukbb": "UK (UKBB)", "saudi": "Saudi Arabia", "gme": "Middle East (GME)"}
# kept in the tables, left out of the figures: UCSC GenomeAsia frequencies are unreliable (common variants
# with REF matching hg38 but absent from both TOPMed and gnomAD, ~25% even at AF 10-90%)
FLAG = {"genomeasia"}
GCLS = ["absent", "qcfail", "rare", "low", "poly"]
GCOL = {"absent": "#d62728", "qcfail": "#bbbbbb", "rare": "#ff9f40", "low": "#f2d15c", "poly": "#2ca02c"}
GLAB = {"absent": "not in gnomAD", "qcfail": "in gnomAD, no AF (QC fail)", "rare": "gnomAD <5%",
        "low": "gnomAD 5-10%", "poly": "gnomAD >10%"}

c = pd.read_csv(os.path.join(D, "counts.tsv"), sep="\t")
ov = c[c.gnomad_view == "overall"]


def save(fig, name):
    fig.savefig(os.path.join(F, name + ".png"), dpi=150, bbox_inches="tight")
    fig.savefig(os.path.join(F, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------- Step 1: share of each dataset not in gnomAD v4.1 (exome + genome) ----------
rows = []
for (ds, t), g in ov.groupby(["dataset", "type"]):
    tot = g.variants.sum()
    ab = g[g.gnomad_class == "absent"].variants.sum()
    common = g[g.ds_class.isin(["r5", "low", "poly"])]
    ctot, cab = common.variants.sum(), common[common.gnomad_class == "absent"].variants.sum()
    rows.append(dict(dataset=ds, type=t, variants=tot, not_in_gnomad=ab, pct_not_in_gnomad=100 * ab / tot if tot else np.nan,
                     variants_af_ge_1pct=ctot, not_in_gnomad_af_ge_1pct=cab,
                     pct_not_in_gnomad_af_ge_1pct=100 * cab / ctot if ctot else np.nan))
s1 = pd.DataFrame(rows).sort_values(["type", "pct_not_in_gnomad"], ascending=[True, False])
s1.to_csv(os.path.join(T, "step1_unique_vs_gnomad.tsv"), sep="\t", index=False, float_format="%.2f")

for t in ["snv", "indel"]:
    x = s1[(s1.type == t) & ~s1.dataset.isin(FLAG)].sort_values("not_in_gnomad_af_ge_1pct")
    fig, ax = plt.subplots(1, 3, figsize=(22, 0.42 * len(x) + 2.5), sharey=True, constrained_layout=True)
    y = np.arange(len(x))
    ax[0].barh(y, x.pct_not_in_gnomad, color="#7f7f7f")
    ax[0].set_xlabel("% of all variants\nnot in gnomAD v4.1")
    ax[1].barh(y, x.pct_not_in_gnomad_af_ge_1pct, color="#1f77b4")
    ax[1].set_xlabel("% of variants with AF >= 1%\nnot in gnomAD")
    ax[2].barh(y, x.not_in_gnomad_af_ge_1pct.clip(lower=1), color="#d62728")
    ax[2].set_xscale("log")
    ax[2].set_xlabel("variants with AF >= 1% not in gnomAD\n(count, log scale)")
    ax[0].set_yticks(y)
    ax[0].set_yticklabels(x.dataset, fontsize=10)
    fig.suptitle(f"{t.upper()}, {REGT}: what each dataset adds beyond gnomAD v4.1")
    save(fig, f"step1_{t}_unique_vs_gnomad")

# ---------- Step 2: dataset polymorphisms (>10%) and rare variants (<5%) by gnomAD class ----------
c["ds_group"] = c.ds_class.map({"poly": "polymorphism (>10%)", "r1": "rare (<5%)", "r5": "rare (<5%)",
                                "low": "low (5-10%)", "noaf": "no AF"})
s2 = (c[c.ds_group.isin(["polymorphism (>10%)", "rare (<5%)"])]
      .groupby(["dataset", "type", "ds_group", "gnomad_view", "gnomad_class"]).variants.sum().unstack(fill_value=0))
for k in GCLS:
    if k not in s2:
        s2[k] = 0
s2 = s2[GCLS]
pct = s2.div(s2.sum(axis=1), axis=0) * 100
out = s2.join(pct.add_prefix("pct_")).reset_index()
out.to_csv(os.path.join(T, "step2_poly_rare_vs_gnomad.tsv"), sep="\t", index=False, float_format="%.2f")

for view in ["overall", "popmax", "matched"]:
    for grp in ["polymorphism (>10%)", "rare (<5%)"]:
        n = len(pct.xs((grp, view), level=["ds_group", "gnomad_view"]).index.get_level_values(0).unique())
        fig, axes = plt.subplots(1, 2, figsize=(22, 0.42 * n + 3), sharey=False, constrained_layout=True)
        for ax, t in zip(axes, ["snv", "indel"]):
            p = pct.xs((t, grp, view), level=["type", "ds_group", "gnomad_view"], drop_level=True)
            p = p[~p.index.isin(FLAG)]
            if p.empty:
                ax.set_visible(False)
                continue
            p = p.sort_values("absent")
            left = np.zeros(len(p))
            for k in GCLS:
                ax.barh(p.index, p[k], left=left, color=GCOL[k], label=GLAB[k])
                left += p[k].values
            ax.set_xlim(0, 100)
            ax.set_xlabel("% of the dataset's variants")
            ax.tick_params(axis="y", labelsize=10)
            ax.set_title(t.upper())
        h, l = axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="outside lower center", ncol=5, fontsize=10)
        fig.suptitle(f"{REGT}: dataset {grp} variants by gnomAD v4.1 class ({view} AF)")
        save(fig, f"step2_{'poly' if grp.startswith('poly') else 'rare'}_{view}")

# ---------- Step 3: markers + population structure ----------
m = pd.read_csv(os.path.join(D, "markers.tsv.gz"), sep="\t", dtype={"chr": str},
                usecols=["chr", "ref", "alt", "dataset", "ds_af", "gnomad_popmax"])
m["type"] = np.where((m.ref == "-") | (m.alt == "-"), "indel",
                     np.where((m.ref.str.len() == 1) & (m.alt.str.len() == 1), "snv", "mnv"))
m["ds_af"] = pd.to_numeric(m.ds_af, errors="coerce")
m["gnomad"] = np.where(m.gnomad_popmax == ".", "absent/QC-fail", "rare (<5%)")
mk = (m.assign(cls=np.where(m.ds_af > 0.10, ">10%", "5-10%"))
      .groupby(["dataset", "type", "cls", "gnomad"]).size().unstack(["cls", "gnomad"], fill_value=0))
mk.to_csv(os.path.join(T, "step3_markers_common_in_dataset_rare_or_absent_in_gnomad.tsv"), sep="\t")
snv = m[(m.type == "snv") & (m.ds_af > 0.10) & m.dataset.isin(COUNTRY)]
cnt = snv.groupby(["dataset", "gnomad"]).size().unstack(fill_value=0)
cnt = cnt.loc[cnt.sum(axis=1).sort_values().index]
fig, ax = plt.subplots(figsize=(13, 6), constrained_layout=True)
left = np.zeros(len(cnt))
for k, col in [("absent/QC-fail", "#d62728"), ("rare (<5%)", "#ff9f40")]:
    if k in cnt:
        ax.barh([COUNTRY[d] for d in cnt.index], cnt[k], left=left, color=col, label=f"gnomAD popmax: {k}")
        left += cnt[k].values
ax.set_xlabel("SNVs with AF > 10% in the country dataset")
ax.set_title(f"Marker SNVs, {REGT}:\ncommon locally, rare or absent in every gnomAD v4.1 group")
ax.legend()
save(fig, "step3_country_marker_snvs")

pc = pd.read_csv(os.path.join(D, "pca_snvs.tsv.gz"), sep="\t", dtype={"chr": str}, na_values=["NA", "."])
popcols = [x for x in pc.columns if x.startswith("gnomad_")]
labels = {**COUNTRY, **{p: "gnomAD " + p.split("_")[1].upper() for p in popcols}}


def pca_plot(cols, name, title):
    x = pc.dropna(subset=cols)
    x = x[(x[cols] >= 0.01).any(axis=1)]
    x = x.assign(win=x.chr + ":" + (x.start // 100000).astype(str)).drop_duplicates("win")   # 1 SNV / 100 kb
    a = x[cols].to_numpy().T                                                                   # populations x SNVs
    z = a - a.mean(axis=0)
    u, sv, _ = np.linalg.svd(z, full_matrices=False)
    var = sv ** 2 / (sv ** 2).sum() * 100
    pcs = u * sv
    pd.DataFrame(pcs[:, :4], index=cols, columns=["PC1", "PC2", "PC3", "PC4"]).to_csv(os.path.join(T, f"{name}_coords.tsv"), sep="\t")
    # pairwise distance: mean |AF difference|
    dist = np.abs(a[:, None, :] - a[None, :, :]).mean(axis=2)
    pd.DataFrame(dist, index=cols, columns=cols).to_csv(os.path.join(T, f"{name}_mean_abs_af_diff.tsv"), sep="\t", float_format="%.4f")
    fig, ax = plt.subplots(1, 2, figsize=(24, 10), constrained_layout=True, gridspec_kw={"width_ratios": [1.15, 1]})
    cmap = plt.get_cmap("tab20")
    for i, col in enumerate(cols):
        g = col.startswith("gnomad_")
        ax[0].scatter(pcs[i, 0], pcs[i, 1], s=170, marker="s" if g else "o", color=cmap(i % 20),
                      edgecolor="black", linewidth=0.6, label=f"{i + 1}. {labels[col]}", zorder=3)
        ax[0].annotate(str(i + 1), (pcs[i, 0], pcs[i, 1]), fontsize=7, ha="center", va="center", zorder=4)
    ax[0].legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=9, frameon=False,
                 title="circle = country dataset\nsquare = gnomAD group")
    ax[0].grid(alpha=0.3)
    ax[0].set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax[0].set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax[0].set_title(f"{title}\n{a.shape[1]:,} SNVs (AF >= 1% somewhere, 1 per 100 kb)")
    im = ax[1].imshow(dist, cmap="viridis")
    ax[1].set_xticks(range(len(cols)))
    ax[1].set_xticklabels([labels[c] for c in cols], rotation=90, fontsize=9)
    ax[1].set_yticks(range(len(cols)))
    ax[1].set_yticklabels([labels[c] for c in cols], fontsize=9)
    ax[1].set_title("mean |AF difference|")
    fig.colorbar(im, ax=ax[1], fraction=0.046)
    save(fig, name)


# datasets with AF for at least 20% of the matrix SNVs in this region (GME has none outside exons)
have = [k for k in COUNTRY if k in pc.columns and pc[k].notna().mean() >= 0.2]
pca_plot(have + popcols, "step3_pca", f"Population PCA on allele frequencies, {REGT}")

# ---------- Step 4: ClinVar ----------
cvf = os.path.join(D, "clinvar_hits.tsv.gz")
if os.path.exists(cvf):
    cv = pd.read_csv(cvf, sep="\t", dtype={"chr": str})
    cv["popmax"] = pd.to_numeric(cv.gnomad_popmax, errors="coerce")
    cv["gnomad"] = np.where(cv.popmax.isna(), "absent/QC-fail",
                            np.where(cv.popmax < 0.05, "rare (<5%)", np.where(cv.popmax <= 0.10, "5-10%", ">10%")))
    cv["max_af_class"] = np.where(cv.max_af > 0.10, ">10%", "5-10%")
    t4 = cv.groupby(["clnsig_group", "max_af_class", "gnomad"]).size().unstack(fill_value=0)
    t4.to_csv(os.path.join(T, "step4_clinvar_common_somewhere_by_gnomad.tsv"), sep="\t")
    flag = cv[cv.clnsig_group.isin(["P/LP", "VUS", "Conflicting"]) & cv.gnomad.isin(["absent/QC-fail", "rare (<5%)"])].copy()
    # independent support: datasets with AF >= 5%, gnomAD v2 / ExAC counted as gnomAD (not independent).
    # A single dataset is usually a data error (ALFA MYH7 74%, SVatalog CFTR F508del 77%, GME BRCA2 del 47%).
    DUP = {"gnomad_exome_v2", "gnomad_genome_v2", "exac03"}
    flag["support"] = flag.datasets_af_ge_5pct.map(lambda x: len({i.split(":")[0] for i in str(x).split(";")} - DUP))
    flag.sort_values(["clnsig_group", "support", "max_af"], ascending=[True, False, False]).to_csv(
        os.path.join(T, "step4_clinvar_PLP_VUS_common_somewhere_rare_in_gnomad.tsv"), sep="\t", index=False)
    rep = flag[flag.support >= 2]
    rep.sort_values(["clnsig_group", "support", "max_af"], ascending=[True, False, False]).to_csv(
        os.path.join(T, "step4_clinvar_flagged_supported_by_2plus_datasets.tsv"), sep="\t", index=False)
    pd.crosstab(flag.clnsig_group, flag.support.clip(upper=4).map(lambda v: "4+" if v == 4 else str(v))).to_csv(
        os.path.join(T, "step4_clinvar_flagged_by_support.tsv"), sep="\t")
    fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
    sup = pd.crosstab(flag.clnsig_group, flag.support.clip(upper=4))
    sup = sup.reindex([k for k in ["P/LP", "Conflicting", "VUS"] if k in sup.index])
    left = np.zeros(len(sup))
    for v, col in zip(sup.columns, ["#bbbbbb", "#d62728", "#ff9f40", "#2ca02c", "#1f77b4"]):
        lab = {0: "only gnomAD v2/ExAC", 1: "1 dataset (likely data error)", 4: "4+ datasets"}.get(v, f"{v} datasets")
        ax.barh(sup.index, sup[v], left=left, color=col, label=lab)
        left += sup[v].values
    ax.set_xlabel("ClinVar variants with AF >= 5% somewhere and gnomAD popmax < 5% (or absent)")
    ax.set_title(f"{REGT}: how many independent datasets support the high frequency")
    ax.legend(fontsize=9)
    save(fig, "step4_clinvar_flagged_by_support")
    per = []
    for _, r in rep.iterrows():
        for item in str(r.datasets_af_ge_5pct).split(";"):
            ds = item.split(":")[0]
            per.append((ds, r.clnsig_group))
    per = pd.DataFrame(per, columns=["dataset", "clnsig_group"]).value_counts().unstack(fill_value=0)
    per.to_csv(os.path.join(T, "step4_clinvar_flagged_supported_per_dataset.tsv"), sep="\t")
    per = per[~per.index.isin(FLAG)]
    cols = [k for k in ["P/LP", "Conflicting", "VUS"] if k in per]
    per = per.loc[per[cols].sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(1, len(cols), figsize=(8 * len(cols), 0.42 * len(per) + 2.5), sharey=True, constrained_layout=True)
    for a_, k in zip(np.atleast_1d(ax), cols):
        a_.barh(per.index, per[k], color={"P/LP": "#d62728", "Conflicting": "#9467bd", "VUS": "#ff9f40"}[k])
        a_.set_title(k)
        a_.tick_params(axis="y", labelsize=10)
    fig.supxlabel("ClinVar variants with AF >= 5% in the dataset and in at least one other dataset, gnomAD popmax < 5% (or absent)")
    fig.suptitle(f"ClinVar 2026-11-15, {REGT}: common in a population, rare in gnomAD v4.1")
    save(fig, "step4_clinvar_flagged_supported_per_dataset")
print("REPORT_DONE")
