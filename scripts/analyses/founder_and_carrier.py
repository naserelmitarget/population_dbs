#!/usr/bin/env python3
"""Founder / population-enriched pathogenic variants and recessive carrier frequencies (exonic, no raw files).

usage: founder_and_carrier.py <extra_analyses dir>
input : <dir>/shared/clinvar_plp_population_af.tsv.gz   (prep_clinvar_plp.sh)
        ClinGen gene-disease validity (BIAS data release) for the mode of inheritance
output: <dir>/founder_variants/{tables,figures}, <dir>/carrier_frequency/{tables,figures}

Population datasets used: the frequency-study datasets except disease cohorts (SPARK, SFARI WGS, GREGoR, GA4K)
and earlier gnomAD releases (gnomAD v2, ExAC). Dataset AFs come from the study (AF ignored when AN < 100).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = sys.argv[1]
CLINGEN = os.environ.get("CLINGEN", "/home/NASER/annotation/acmg/bias_hg38_current/hg38_clingen_gene_disease_validity.csv")
EXCLUDE = {"spark_wes", "sfari_wgs", "gregor", "ga4k", "gnomad_exome_v2", "gnomad_genome_v2", "exac03"}
LABEL = {"tommo": "Japan (ToMMo)", "kova": "Korea (KOVA)", "wbbc": "China (WBBC)", "sg10k_npm": "Singapore (SG10K)",
         "mcps": "Mexico (MCPS)", "finngen": "Finland (FinnGen)", "swegen": "Sweden (SweGen)", "gonl": "Netherlands (GoNL)",
         "ukbb": "UK (UKBB)", "saudi": "Saudi Arabia", "gme": "Middle East (GME)", "topmed": "TOPMed", "allofus": "All of Us",
         "alfa": "ALFA", "kaviar": "Kaviar", "1000g2015aug": "1000 Genomes", "hgdp1kg": "HGDP+1kG", "hostseq": "HostSeq (Canada)",
         "colorsdb": "CoLoRSdb", "svatalog": "SVatalog"}
MATCH = {"tommo": "eas", "kova": "eas", "wbbc": "eas", "sg10k_npm": "eas", "mcps": "amr", "finngen": "fin",
         "swegen": "nfe", "gonl": "nfe", "ukbb": "nfe", "saudi": "mid", "gme": "mid"}
GROUPS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]
ENRICH, MIN_AF = 10, 0.001


def save(fig, folder, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(folder, f"{name}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def mk(name):
    t, f = os.path.join(D, name, "tables"), os.path.join(D, name, "figures")
    os.makedirs(t, exist_ok=True)
    os.makedirs(f, exist_ok=True)
    return t, f


v = pd.read_csv(os.path.join(D, "shared", "clinvar_plp_population_af.tsv.gz"), sep="\t", dtype={"chr": str})
v["stars"] = pd.to_numeric(v.stars, errors="coerce").fillna(0).astype(int)
v = v[v.exonic_in_study.astype(str) == "True"].copy()

# long table: one row per variant x population dataset with an AF
rows = []
for idx, s in zip(v.index, v.datasets.fillna("")):
    for item in s.split(";"):
        if "=" not in item:
            continue
        d, af = item.split("=", 1)
        if d in EXCLUDE or af in ("", "."):
            continue
        rows.append((idx, d, float(af)))
L = pd.DataFrame(rows, columns=["idx", "dataset", "af"])
L = L.join(v[["gnomad_popmax", "gnomad_overall"] + [f"gnomad_{g}" for g in GROUPS]], on="idx")

# ---------------- founder / population-enriched P/LP ----------------
T, F = mk("founder_variants")
L["fold_vs_popmax"] = L.af / L.gnomad_popmax.replace(0, np.nan)
L["enriched"] = (L.af >= MIN_AF) & (L.gnomad_popmax.isna() | (L.af >= ENRICH * L.gnomad_popmax))
L["matched_group_af"] = np.nan
for d, g in MATCH.items():
    L.loc[L.dataset == d, "matched_group_af"] = L.loc[L.dataset == d, f"gnomad_{g}"]
e = L[L.enriched]
agg = e.groupby("idx").agg(enriched_in=("dataset", lambda x: ";".join(sorted(x))), n_datasets=("dataset", "nunique"),
                           max_af=("af", "max"))
cand = v.join(agg, how="inner")
cand["evidence"] = np.where(cand.n_datasets >= 2, "replicated (>= 2 datasets)", "single dataset")
# technical flags: unreliable-gene list and GIAB segmental duplications / low mappability (unreliable_regions.sh)
ug_fp = os.path.join(D, "unreliable_regions", "tables", "unreliable_genes.tsv")
ug_set = set(pd.read_csv(ug_fp, sep="\t").symbol.dropna()) if os.path.exists(ug_fp) else set()
cand["unreliable_gene"] = cand.gene.isin(ug_set)
cand["giab_difficult"] = False
for b in ("segdups", "lowmappabilityall"):
    fp = os.path.join(D, "unreliable_regions", "tmp", f"{b}.bed")
    if not os.path.exists(fp):
        continue
    bed = pd.read_csv(fp, sep="\t", header=None, names=["chr", "s", "e"], dtype={"chr": str})
    for c, g in cand.groupby("chr"):
        iv = bed[bed.chr == str(c)].sort_values("s")
        if iv.empty:
            continue
        pos = pd.to_numeric(g.start).to_numpy()
        i = np.searchsorted(iv.s.to_numpy(), pos, side="right") - 1
        inside = (i >= 0) & (pos <= iv.e.to_numpy()[np.clip(i, 0, None)])
        cand.loc[g.index[inside], "giab_difficult"] = True
# a P/LP allele above 5% (ACMG BA1 range) in a population is implausible for a penetrant disease allele
cand["implausible_af"] = cand.max_af > 0.05
cand["high_confidence"] = (cand.n_datasets >= 2) & (cand.stars >= 1) & ~cand.unreliable_gene & ~cand.giab_difficult \
                          & ~cand.implausible_af
cols = ["chr", "start", "end", "ref", "alt", "gene", "clnsig", "stars", "consequence", "condition", "max_af", "n_datasets",
        "enriched_in", "evidence", "unreliable_gene", "giab_difficult", "implausible_af", "high_confidence", "gnomad_overall", "gnomad_popmax"] + [f"gnomad_{g}" for g in GROUPS] + ["datasets", "clinvar_id"]
cand = cand.sort_values(["n_datasets", "stars", "max_af"], ascending=[False, False, False])
cand[cols].to_csv(os.path.join(T, "founder_candidates_all.tsv"), sep="\t", index=False)
pairs = e.join(v[["chr", "start", "ref", "alt", "gene", "clnsig", "stars"]], on="idx")
pairs[["chr", "start", "ref", "alt", "gene", "clnsig", "stars", "dataset", "af", "gnomad_popmax", "matched_group_af",
       "fold_vs_popmax"]].sort_values(["dataset", "af"], ascending=[True, False]).to_csv(
    os.path.join(T, "founder_variant_dataset_pairs.tsv"), sep="\t", index=False, float_format="%.6g")
cand[(cand.n_datasets >= 2) & (cand.stars >= 1)][cols].to_csv(os.path.join(T, "founder_candidates_replicated_1star.tsv"), sep="\t", index=False)
cand[cand.high_confidence][cols].to_csv(os.path.join(T, "founder_candidates_high_confidence.tsv"), sep="\t", index=False)
per = (e.merge(v[["stars"]], left_on="idx", right_index=True)
       .assign(replicated=lambda x: x.idx.map(agg.n_datasets) >= 2)
       .groupby("dataset").agg(enriched_variants=("idx", "nunique"),
                               replicated=("replicated", "sum"),
                               with_1plus_star=("stars", lambda s: int((s >= 1).sum())),
                               genes=("idx", lambda i: v.loc[i, "gene"].nunique())))
per = per.sort_values("enriched_variants", ascending=False)
per.to_csv(os.path.join(T, "founder_summary_per_dataset.tsv"), sep="\t")
fig, ax = plt.subplots(figsize=(12, 0.42 * len(per) + 2.5), constrained_layout=True)
p = per.iloc[::-1]
ax.barh([LABEL.get(d, d) for d in p.index], p.enriched_variants, color="#ff9f40", label="enriched (one dataset or more)")
ax.barh([LABEL.get(d, d) for d in p.index], p.replicated, color="#d62728", label="also enriched in another dataset")
ax.set_xlabel(f"ClinVar P/LP variants with AF >= {MIN_AF:.1%} in the dataset and >= {ENRICH}x gnomAD popmax (or absent from gnomAD)")
ax.set_title("Population-enriched pathogenic variants (exonic)")
ax.legend()
save(fig, F, "founder_variants_per_dataset")
print(f"founder: {len(cand)} enriched P/LP variants; replicated {int((cand.n_datasets >= 2).sum())}; "
      f"replicated with >=1 star {int(((cand.n_datasets >= 2) & (cand.stars >= 1)).sum())}; high confidence {int(cand.high_confidence.sum())}")

# ---------------- carrier frequency, autosomal recessive genes ----------------
T, F = mk("carrier_frequency")
cg = pd.read_csv(CLINGEN, skiprows=[0, 1, 2, 3, 5], dtype=str)
cg.columns = [c.strip() for c in cg.columns]
ar = cg[(cg.MOI == "AR") & cg.CLASSIFICATION.isin(["Definitive", "Strong", "Moderate"])]
ar_genes = set(ar["GENE SYMBOL"])
# quality rules for carrier estimates:
#  - only population cohorts: aggregates / tiny / long-read reference sets are left out
#  - genes in the unreliable-gene list (unreliable_regions.py) are left out
#  - a P/LP allele's AF >= 1% counts only if a second source (another dataset or a gnomAD group) also has it >= 1%
CARRIER_SKIP = {"alfa", "kaviar", "svatalog", "colorsdb"}
ug = os.path.join(D, "unreliable_regions", "tables", "unreliable_genes.tsv")
bad_genes = set(pd.read_csv(ug, sep="\t").symbol.dropna()) if os.path.exists(ug) else set()
#  - coding / canonical-splice consequences only (ClinVar MC): intronic repeat alleles such as CFTR poly-T/TG are
#    classed P/LP with variable penetrance and would inflate carrier rates several-fold
CODING = {"nonsense", "frameshift_variant", "missense_variant", "splice_donor_variant", "splice_acceptor_variant",
          "inframe_deletion", "inframe_insertion", "inframe_indel", "initiator_codon_variant", "stop_lost", "start_lost"}
q = v[v.gene.isin(ar_genes - bad_genes) & (v.stars >= 1) & v.consequence.isin(CODING)
      & ~v.clnsig.str.contains("low_penetrance|risk", case=False, na=False)]
Lq = L[L.idx.isin(q.index) & ~L.dataset.isin(CARRIER_SKIP)].copy()
n_hi = Lq[Lq.af >= 0.01].groupby("idx").size().reindex(Lq.idx).fillna(0).to_numpy() \
       + (Lq[[f"gnomad_{g}" for g in GROUPS]].max(axis=1) >= 0.01).to_numpy()
dropped = (Lq.af >= 0.01) & (n_hi < 2)
print(f"carrier: {int(dropped.sum())} unreplicated dataset AFs >= 1% dropped; {len(bad_genes)} unreliable genes excluded")
Lq = Lq[~dropped]
Lq["gene"] = Lq.idx.map(v.gene)
ds_q = Lq.groupby(["gene", "dataset"]).af.sum().unstack(fill_value=0)
gn = pd.DataFrame({f"gnomAD {g.upper()}": q.groupby("gene")[f"gnomad_{g}"].sum(min_count=1) for g in GROUPS}).fillna(0)
Q = ds_q.join(gn, how="outer").fillna(0).clip(upper=1)
CF = 2 * Q * (1 - Q)                                   # carrier frequency (heterozygotes), Hardy-Weinberg
AFF = Q ** 2 * 1e5                                     # expected affected per 100,000 births (random mating)
Q.to_csv(os.path.join(T, "pathogenic_allele_frequency_per_gene.tsv"), sep="\t", float_format="%.6g")
CF.to_csv(os.path.join(T, "carrier_frequency_per_gene.tsv"), sep="\t", float_format="%.6g")
AFF.to_csv(os.path.join(T, "expected_affected_per_100k_births.tsv"), sep="\t", float_format="%.4g")
cum = (1 - (1 - CF).prod(axis=0)).rename("P(carrier of >= 1 AR condition)")
n_genes = (CF > 0.01).sum(axis=0).rename("genes with carrier frequency > 1%")
summ = pd.concat([cum, n_genes, (CF > 0.02).sum(axis=0).rename("genes > 2%")], axis=1).sort_values(cum.name, ascending=False)
summ.to_csv(os.path.join(T, "carrier_summary_per_population.tsv"), sep="\t", float_format="%.4f")
top = []
for c in CF.columns:
    t = CF[c].sort_values(ascending=False).head(15)
    top += [(c, g, f) for g, f in t.items() if f > 0]
pd.DataFrame(top, columns=["population", "gene", "carrier_frequency"]).to_csv(
    os.path.join(T, "top15_genes_per_population.tsv"), sep="\t", index=False, float_format="%.5f")
show = [d for d in ["tommo", "kova", "wbbc", "sg10k_npm", "mcps", "finngen", "swegen", "gonl", "ukbb", "saudi", "gme"] if d in CF] \
       + [c for c in CF.columns if c.startswith("gnomAD")]
genes = CF[show].max(axis=1).sort_values(ascending=False).head(40).index
H = CF.loc[genes, show] * 100
fig, ax = plt.subplots(figsize=(0.6 * len(show) + 5, 0.3 * len(genes) + 3), constrained_layout=True)
im = ax.imshow(np.log10(H.values + 0.01), cmap="YlOrRd", aspect="auto", vmin=-2, vmax=1.3)
ax.set_xticks(range(len(show)))
ax.set_xticklabels([LABEL.get(c, c) for c in show], rotation=90, fontsize=9)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=8)
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        if H.iat[i, j] >= 1:
            ax.text(j, i, f"{H.iat[i, j]:.0f}", ha="center", va="center", fontsize=6)
cb = fig.colorbar(im, ax=ax, fraction=0.03)
cb.set_label("carrier frequency, % (log10 scale; numbers shown when >= 1%)")
ax.set_title("Carrier frequency of autosomal recessive disorders (ClinGen AR genes, ClinVar P/LP >= 1 star, exonic)")
save(fig, F, "carrier_frequency_top_genes")
fig, ax = plt.subplots(figsize=(11, 0.4 * len(summ) + 2), constrained_layout=True)
s = summ.iloc[::-1]
ax.barh([LABEL.get(c, c) for c in s.index], s[cum.name] * 100, color="#1f77b4")
ax.set_xlabel("% of people expected to carry a known pathogenic allele in at least one AR gene")
ax.set_title("Cumulative carrier probability per population")
save(fig, F, "carrier_cumulative_per_population")
print(f"carrier: {len(ar_genes)} AR genes, {q.gene.nunique()} with P/LP >= 1 star among exonic variants")
