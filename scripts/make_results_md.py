#!/usr/bin/env python3
"""Write reports/sections/02_results.md from the result tables in results/.

usage: python3 scripts/make_results_md.py   (run from the repository root)
Every number in the results document comes from a table in results/; edit the narrative here, not in the .md.
"""
import pandas as pd

R = "results"
S = "02_frequency_vs_gnomad"   # step 1-4 results, per region
COUNTRY = {"tommo": "ToMMo (Japan)", "kova": "KOVA (Korea)", "wbbc": "WBBC (China)", "sg10k_npm": "SG10K (Singapore)",
           "mcps": "MCPS (Mexico)", "saudi": "Saudi Arabia", "gme": "GME (Middle East)", "swegen": "SweGen (Sweden)",
           "gonl": "GoNL (Netherlands)", "finngen": "FinnGen (Finland)", "ukbb": "UK Biobank"}
GROUPED = list(COUNTRY) + ["topmed", "allofus"]
out = []


def w(s=""):
    out.append(s)


def table(df, cols, heads, fmt=None):
    fmt = fmt or {}
    w("| " + " | ".join(heads) + " |")
    w("|" + "|".join("---:" if i else "---" for i in range(len(heads))) + "|")
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(fmt[c](v) if c in fmt and pd.notna(v) else ("" if pd.isna(v) else str(v)))
        w("| " + " | ".join(cells) + " |")
    w()


pct = lambda v: f"{v:.2f}%"
num = lambda v: f"{int(v):,}"

w("# Results")
w()
w("All numbers below are generated from the tables in [`results/`](../../results) by `scripts/make_results_md.py`. "
  "Definitions of every class and view are in [Materials and methods](01_materials_and_methods.md).")
w()

# ---------------- union ----------------
hm = pd.read_csv(f"{R}/01_union/summary_how_many.tsv", sep="\t")
tot = hm.variants.sum()
one = hm[hm.how_many == 1].variants.sum()
mid = hm[(hm.how_many >= 2) & (hm.how_many <= 5)].variants.sum()
hi = hm[hm.how_many >= 6].variants.sum()
w("## 1. The union of all datasets")
w()
w(f"The union table holds **{tot:,}** distinct variants (hg38). {one / tot:.1%} are in one dataset only, "
  f"{mid / tot:.1%} in 2–5 datasets and {hi / tot:.1%} in 6 or more. "
  f"The highest number of datasets sharing one variant is {hm.how_many.max()}.")
w()
pdset = pd.read_csv(f"{R}/01_union/summary_per_dataset.tsv", sep="\t").sort_values("variants", ascending=False)
table(pdset, ["dataset", "variants", "only_in_this_dataset"], ["Dataset", "Variants in union", "Only in this dataset"],
      {"variants": num, "only_in_this_dataset": num})

# ---------------- step 1 ----------------
w("## 2. What each dataset adds beyond gnomAD v4.1 (step 1)")
w()
w("Share of each dataset's variants that are absent from gnomAD v4.1, for all variants and for variants with AF ≥ 1% "
  "in the dataset. Exonic: gnomAD exomes + genomes. Non-exonic: gnomAD genomes only.")
w()
for reg, title in [("exonic", "Exonic regions"), ("non_exonic", "Non-exonic regions")]:
    s1 = pd.read_csv(f"{R}/{S}/{reg}/tables/step1_unique_vs_gnomad.tsv", sep="\t")
    for t in ["snv", "indel"]:
        x = s1[s1.type == t].sort_values("pct_not_in_gnomad_af_ge_1pct", ascending=False)
        w(f"### {title}, {t.upper()}s")
        w()
        table(x, ["dataset", "variants", "pct_not_in_gnomad", "variants_af_ge_1pct", "not_in_gnomad_af_ge_1pct",
                  "pct_not_in_gnomad_af_ge_1pct"],
              ["Dataset", "Variants", "% not in gnomAD", "Variants AF ≥ 1%", "AF ≥ 1% not in gnomAD", "% AF ≥ 1% not in gnomAD"],
              {"variants": num, "variants_af_ge_1pct": num, "not_in_gnomad_af_ge_1pct": num,
               "pct_not_in_gnomad": pct, "pct_not_in_gnomad_af_ge_1pct": pct})
    w(f"![{title}: SNVs](../../results/{S}/{reg}/figures/step1_snv_unique_vs_gnomad.png)")
    w()
    w(f"![{title}: indels](../../results/{S}/{reg}/figures/step1_indel_unique_vs_gnomad.png)")
    w()


# ---------------- step 2 ----------------
def step2(reg, group, cls_cols):
    d = pd.read_csv(f"{R}/{S}/{reg}/tables/step2_poly_rare_vs_gnomad.tsv", sep="\t")
    d = d[(d.ds_group == group) & d.dataset.isin(GROUPED) & d.type.isin(["snv", "indel"])]
    rows = []
    for (ds, t), g in d.groupby(["dataset", "type"]):
        r = {"dataset": COUNTRY.get(ds, ds), "type": t.upper()}
        ov = g[g.gnomad_view == "overall"]
        r["n"] = int(ov[["absent", "qcfail", "rare", "low", "poly"]].sum(axis=1).iloc[0])
        for v in ["overall", "popmax", "matched"]:
            x = g[g.gnomad_view == v]
            for c in cls_cols:
                r[f"{v}_{c}"] = (x[f"pct_{c}"].iloc[0] if c != "absqc" else x.pct_absent.iloc[0] + x.pct_qcfail.iloc[0]) if len(x) else None
        rows.append(r)
    return pd.DataFrame(rows).sort_values(["type", "dataset"], ascending=[False, True])


w("## 3. How gnomAD classifies each dataset's polymorphisms (step 2)")
w()
w("Of each dataset's polymorphisms (AF > 10%), the share that gnomAD calls rare (< 5%), and the share absent or QC-failed, "
  "under the three gnomAD views: overall AF, the population maximum (popmax) and the dataset's matching gnomAD group.")
w()
for reg, title in [("exonic", "Exonic"), ("non_exonic", "Non-exonic")]:
    x = step2(reg, "polymorphism (>10%)", ["rare", "absqc"])
    w(f"### {title}")
    w()
    table(x, ["dataset", "type", "n", "overall_rare", "popmax_rare", "matched_rare", "popmax_absqc"],
          ["Dataset", "Type", "Polymorphisms", "Rare in gnomAD (overall)", "Rare (popmax)", "Rare (matched group)", "Absent / QC fail (popmax)"],
          {"n": num, "overall_rare": pct, "popmax_rare": pct, "matched_rare": pct, "popmax_absqc": pct})
w("![Exonic polymorphisms by matched gnomAD group](../../results/02_frequency_vs_gnomad/exonic/figures/step2_poly_matched.png)")
w()
w("![Exonic polymorphisms by gnomAD popmax](../../results/02_frequency_vs_gnomad/exonic/figures/step2_poly_popmax.png)")
w()

w("## 4. The opposite direction: rare in the dataset, polymorphic in gnomAD (step 2)")
w()
w("Of each dataset's rare variants (AF < 5%), the share that gnomAD calls polymorphic (> 10%).")
w()
for reg, title in [("exonic", "Exonic"), ("non_exonic", "Non-exonic")]:
    x = step2(reg, "rare (<5%)", ["poly"])
    w(f"### {title}")
    w()
    table(x, ["dataset", "type", "n", "overall_poly", "matched_poly", "popmax_poly"],
          ["Dataset", "Type", "Rare variants", "Polymorphic in gnomAD (overall)", "Polymorphic (matched group)", "Polymorphic (popmax)"],
          {"n": num, "overall_poly": pct, "matched_poly": pct, "popmax_poly": pct})
w("![Exonic rare variants by matched gnomAD group](../../results/02_frequency_vs_gnomad/exonic/figures/step2_rare_matched.png)")
w()
w("![Exonic rare variants by gnomAD popmax](../../results/02_frequency_vs_gnomad/exonic/figures/step2_rare_popmax.png)")
w()

# ---------------- step 3 ----------------
w("## 5. Population markers and structure (step 3)")
w()
mk = pd.read_csv(f"{R}/{S}/exonic/tables/step3_markers_common_in_dataset_rare_or_absent_in_gnomad.tsv", sep="\t", header=[0, 1, 2], index_col=[0, 1])
w("Exonic SNVs with AF > 10% in the dataset and gnomAD popmax < 5% or absent, per country dataset:")
w()
rows = []
for ds in COUNTRY:
    if (ds, "snv") in mk.index:
        r = mk.loc[(ds, "snv")]
        a = sum(v for k, v in r.items() if k[0] == ">10%" and k[1].startswith("absent"))
        b = sum(v for k, v in r.items() if k[0] == ">10%" and k[1].startswith("rare"))
        rows.append({"dataset": COUNTRY[ds], "absent": a, "rare": b, "total": a + b})
table(pd.DataFrame(rows).sort_values("total", ascending=False), ["dataset", "total", "absent", "rare"],
      ["Dataset", "Marker SNVs", "Absent / QC fail in gnomAD", "Rare (< 5%) in every gnomAD group"],
      {"total": num, "absent": num, "rare": num})
w("![Marker SNVs](../../results/02_frequency_vs_gnomad/exonic/figures/step3_country_marker_snvs.png)")
w()
w("Population PCA on allele frequencies (exonic SNVs with gnomAD genome AF, AF ≥ 1% somewhere, one per 100 kb) "
  "and the mean absolute AF difference between populations:")
w()
w("![PCA](../../results/02_frequency_vs_gnomad/exonic/figures/step3_pca.png)")
w()

# ---------------- step 4 ----------------
w("## 6. ClinVar (step 4)")
w()
w("ClinVar (2026-11-15) variants with AF ≥ 5% in at least one dataset and gnomAD popmax < 5% or absent, "
  "by the number of independent datasets supporting the high frequency (gnomAD v2 and ExAC counted with gnomAD).")
w()
for reg, title in [("exonic", "Exonic"), ("non_exonic", "Non-exonic")]:
    s = pd.read_csv(f"{R}/{S}/{reg}/tables/step4_clinvar_flagged_by_support.tsv", sep="\t").fillna(0)
    cols = [c for c in s.columns if c != "clnsig_group"]
    s["total"] = s[cols].sum(axis=1)
    w(f"### {title}")
    w()
    table(s, ["clnsig_group"] + cols + ["total"], ["ClinVar group"] + [f"{c} dataset(s)" for c in cols] + ["Total"],
          {c: num for c in cols + ["total"]})
rep = pd.read_csv(f"{R}/{S}/exonic/tables/step4_clinvar_flagged_supported_by_2plus_datasets.tsv", sep="\t")
p = rep[rep.clnsig_group == "P/LP"].sort_values("support", ascending=False)
w("Exonic pathogenic / likely pathogenic variants supported by two or more datasets:")
w()
table(p, ["chr", "start", "ref", "alt", "gene", "clnsig", "review_stars", "max_af", "gnomad_popmax", "datasets_af_ge_5pct"],
      ["Chr", "Start", "Ref", "Alt", "Gene", "ClinVar", "Stars", "Max AF", "gnomAD popmax", "Datasets (AF)"])
w("![ClinVar support](../../results/02_frequency_vs_gnomad/exonic/figures/step4_clinvar_flagged_by_support.png)")
w()

# ---------------- gene level ----------------
w("## 7. Gene level (step 5, exonic, protein-coding genes)")
w()
s = pd.read_csv(f"{R}/03_gene_level/tables/A_gene_vs_gnomad_significant.tsv", sep="\t")
t = (s[s["class"].isin(["protein_altering", "synonymous"])]
     .groupby(["dataset", "class", "af", "call"]).size().unstack(["class", "af", "call"], fill_value=0))
t.columns = [f"{a.replace('protein_altering', 'PA').replace('synonymous', 'syn')} {b} {c}" for a, b, c in t.columns]
t = t.reset_index()
keep = ["dataset"] + [c for c in ["PA common high", "PA common low", "syn common high", "PA all high", "PA all low"] if c in t]
w("Genes significantly above (high) or below (low) the gnomAD expectation (FDR < 5%, ≥ 2-fold, ≥ 5 variants); "
  "PA = protein-altering, syn = synonymous (technical control), common = AF ≥ 1%.")
w()
table(t.sort_values("PA common high", ascending=False), keep, ["Dataset"] + keep[1:])
u = pd.read_csv(f"{R}/03_gene_level/tables/C_high_genes_vs_gnomad_by_uniqueness.tsv", sep="\t")
x = u[(u["class"] == "protein_altering") & (u.af == "common")]
w("Uniqueness of genes high for common protein-altering variants:")
w()
vc = x.pattern.value_counts().rename_axis("pattern").reset_index(name="genes")
table(vc, ["pattern", "genes"], ["Pattern", "Genes"])
w("Genes high in the most datasets:")
w()
table(x.sort_values("n_datasets", ascending=False).head(25), ["symbol", "n_datasets", "groups"],
      ["Gene", "Datasets", "Population groups"])
w("![Genes deviating from gnomAD per dataset](../../results/03_gene_level/figures/A_significant_genes_per_dataset_protein_altering_common.png)")
w()
w("![Top genes heatmap](../../results/03_gene_level/figures/A_heatmap_top_genes_protein_altering_common.png)")
w()
w("![Within-dataset density](../../results/03_gene_level/figures/B_heatmap_high_density_genes_protein_altering.png)")
w()

open("reports/sections/02_results.md", "w").write("\n".join(out))
print("wrote reports/sections/02_results.md")
