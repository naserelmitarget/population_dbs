#!/usr/bin/env python3
"""PGR comparison: dataset coding variants absent from gnomAD v4.1 by consequence, next to the Pakistan Genome
Resource (Koch et al., Nature 2026, Table 1: unique relative to all of gnomAD v4.1, MAC >= 1).

usage: python3 scripts/analyses/pgr_comparison_figure.py   (from the repository root; input from pgr_consequence_absent.sh)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = "results/09_pgr_comparison"
w = pd.read_csv(f"{D}/tables/coding_absent_pct_by_consequence_wide.tsv", sep="\t", index_col=0)
pgr = {"synonymous": 491878 / 1976965 * 100, "missense": 1269538 / 4186236 * 100,
       "inframe": (18250 + 9325) / (53716 + 23964) * 100, "pLoF": 180314 / 370128 * 100}
rows = [{"dataset": "PGR (Pakistan, paper Table 1)", **pgr, "coding variants": 6625995}]
for d in ["tommo", "mcps", "allofus", "wbbc", "kova", "sg10k_npm", "topmed", "gme", "saudi", "swegen", "gonl"]:
    rows.append({"dataset": d, **{k: w.loc[d, k] for k in pgr}, "coding variants": w.loc[d, "coding variants"]})
t = pd.DataFrame(rows)
t["pLoF / synonymous"] = t.pLoF / t.synonymous
t.to_csv(f"{D}/tables/pgr_vs_datasets_absent_from_gnomad.tsv", sep="\t", index=False, float_format="%.2f")
fig, ax = plt.subplots(figsize=(14, 6), constrained_layout=True)
x = np.arange(len(t))
for i, (c, col) in enumerate(zip(pgr, ["#9ecae1", "#4292c6", "#2171b5", "#d62728"])):
    ax.bar(x + (i - 1.5) * 0.2, t[c], width=0.2, color=col, label=c)
ax.set_xticks(x)
ax.set_xticklabels(["PGR\n(paper)"] + list(t.dataset[1:]), fontsize=9)
ax.set_ylabel("% of the dataset's coding variants absent from gnomAD v4.1")
ax.set_title("Absent from gnomAD v4.1 by consequence: the Pakistan Genome Resource (Table 1) and the datasets in this study")
ax.legend()
for ext in ("png", "pdf"):
    fig.savefig(f"{D}/figures/pgr_vs_datasets_absent_by_consequence.{ext}", dpi=150, bbox_inches="tight")
