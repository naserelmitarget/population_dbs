#!/usr/bin/env python3
"""Candidate ancestry-informative SNVs for the Middle East (selection only, not validated on genotypes).

usage: ancestry_markers.py <extra_analyses dir>
input : data/archive/study/{exonic,non_exonic}/pca_snvs.tsv.gz (SNV AFs of country datasets and gnomAD genome groups)
output: <dir>/ancestry_markers/{tables,figures}

Populations: Middle East = mean of the available Middle Eastern sets (Saudi, gnomAD MID; plus GME in exonic regions),
and gnomAD genome groups AFR, AMR, EAS, NFE, SAS. For each SNV:
  In     Rosenberg et al. (2003) informativeness for assignment across the six populations
  ME_gap smallest |AF(Middle East) - AF(other group)|: large only when the Middle East differs from every other group
Quality: autosomal SNVs, AFs within [0, 1], and at least two Middle Eastern sources agreeing within 0.1 (their
disagreement otherwise dominates the ranking). Each list is thinned to one SNV per 500 kb and the top 1,000 kept.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = sys.argv[1]
K = os.environ.get("K", "/home/NASER/annotation/population_study/data/archive") + "/study"
O = os.path.join(D, "ancestry_markers")
T, F = os.path.join(O, "tables"), os.path.join(O, "figures")
os.makedirs(T, exist_ok=True)
os.makedirs(F, exist_ok=True)
OTHER = {"AFR": "gnomad_afr", "AMR": "gnomad_amr", "EAS": "gnomad_eas", "NFE": "gnomad_nfe", "SAS": "gnomad_sas"}


def informativeness(P):
    """Rosenberg In for biallelic markers; P = SNV x population AF matrix."""
    P = np.clip(P, 1e-9, 1 - 1e-9)
    out = 0
    for p in (P, 1 - P):
        m = p.mean(axis=1)
        out = out + (-m * np.log(m) + (p * np.log(p)).mean(axis=1))
    return out


def thin(x, score, bp=500_000):
    x = x.sort_values(score, ascending=False)
    x["bin"] = x.chr.astype(str) + ":" + (x.start // bp).astype(str)
    return x.drop_duplicates("bin").drop(columns="bin")


frames = []
for reg in ("exonic", "non_exonic"):
    x = pd.read_csv(f"{K}/{reg}/pca_snvs.tsv.gz", sep="\t", dtype={"chr": str}, na_values=["NA", "."])
    me_cols = [c for c in ("saudi", "gme", "gnomad_mid") if c in x and x[c].notna().mean() > 0.05]
    x = x.dropna(subset=list(OTHER.values()) + ["gnomad_mid"])
    x = x[x.chr.astype(str).str.fullmatch(r"\d+")]                       # autosomes only (no hemizygous X/Y)
    af_cols = me_cols + list(OTHER.values())
    x = x[((x[af_cols] >= 0) & (x[af_cols] <= 1) | x[af_cols].isna()).all(axis=1)]
    x["ME"] = x[me_cols].mean(axis=1, skipna=True)
    x["ME_sets"] = x[me_cols].notna().sum(axis=1)
    x["ME_spread"] = x[me_cols].max(axis=1) - x[me_cols].min(axis=1)
    # Middle Eastern sources must agree: at least two of them, within 0.1 of each other
    x = x[(x.ME_sets >= 2) & (x.ME_spread <= 0.1)]
    P = x[["ME"] + list(OTHER.values())].to_numpy()
    x["In"] = informativeness(P)
    gaps = np.abs(x[["ME"]].to_numpy() - x[list(OTHER.values())].to_numpy())
    x["ME_gap"] = gaps.min(axis=1)
    x["closest_group"] = np.array(list(OTHER))[gaps.argmin(axis=1)]
    x["region"] = reg
    frames.append(x)
    print(f"{reg}: {len(x):,} SNVs scored; Middle East from {', '.join(me_cols)}")
A = pd.concat(frames, ignore_index=True)
cols = ["region", "chr", "start", "ref", "alt", "ME", "saudi", "gme", "gnomad_mid"] + list(OTHER.values()) + ["In", "ME_gap", "closest_group"]
cols = [c for c in cols if c in A]
panel_in = thin(A, "In").head(1000)
panel_me = thin(A[A.ME_gap > 0], "ME_gap").head(1000)
panel_in[cols].to_csv(os.path.join(T, "panel_continental_top1000_informativeness.tsv"), sep="\t", index=False, float_format="%.4f")
panel_me[cols].to_csv(os.path.join(T, "panel_middle_east_distinct_top1000.tsv"), sep="\t", index=False, float_format="%.4f")
summ = pd.DataFrame({"list": ["continental (In)", "Middle East distinct (ME_gap)"],
                     "SNVs": [len(panel_in), len(panel_me)],
                     "median score": [panel_in.In.median(), panel_me.ME_gap.median()],
                     "exonic": [(panel_in.region == "exonic").sum(), (panel_me.region == "exonic").sum()]})
summ.to_csv(os.path.join(T, "panel_summary.tsv"), sep="\t", index=False, float_format="%.4f")

top = panel_me.head(30)
pops = ["ME"] + list(OTHER.values())
fig, ax = plt.subplots(figsize=(9, 10), constrained_layout=True)
im = ax.imshow(top[pops].to_numpy(), cmap="viridis", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(pops)))
ax.set_xticklabels(["Middle East"] + [f"gnomAD {k}" for k in OTHER], rotation=45, ha="right")
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f"{c}:{s} {r}>{a}" for c, s, r, a in zip(top.chr, top.start, top.ref, top.alt)], fontsize=7)
fig.colorbar(im, ax=ax, fraction=0.04, label="allele frequency")
ax.set_title("Top 30 SNVs distinguishing the Middle East from every other gnomAD group")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(F, f"middle_east_top30_markers.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
ax.hist(A.ME_gap, bins=60, color="#7f7f7f", log=True)
ax.axvline(panel_me.ME_gap.min(), color="#d62728", ls="--", label="top-1,000 cut-off")
ax.set_xlabel("smallest |AF(Middle East) - AF(other gnomAD group)|")
ax.set_ylabel("SNVs (log scale)")
ax.legend()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(F, f"middle_east_gap_distribution.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("ANCESTRY_DONE")
