# Further analyses

Analyses that would add value for clinical and population geneticists, roughly in order of value for effort. Each can
reuse the union, the harmonised per-dataset frequencies and the VEP annotation built here.

**Done** (see [Additional analyses](07_additional_analyses.md)): founder and population-enriched pathogenic variants
(item 5), carrier frequencies (item 4), the "do not trust the frequency" gene and region list (item 7), candidate
ancestry-informative markers (item 10, selection only), and an ACMG reclassification test on three genomes (item 2).
**Need the raw frequency files** (AC, AN or homozygote counts were not kept): items 3 and 8, and the non-coding parts of
items 6 and 9.

## High value

1. **Iranian / Middle Eastern focus.** Add Iranome and in-house Iranian exomes. Repeat steps 1–4 with the Middle East as
   the target population: how many common Iranian variants are rare or absent in gnomAD, GME and Saudi, and how the
   three overlap.
2. **ACMG reclassification impact in real patients.** Annotate a set of patient exomes with gnomAD only, then with gnomAD
   plus the population datasets, and count variants whose BA1, BS1 or PM2 status changes, and whose final class changes.
   Use ClinGen expert-panel gene-specific thresholds where they exist, and the maximum credible population AF (from
   disease prevalence, penetrance and genetic heterogeneity) for BS1.
3. **Filtering allele frequency with confidence intervals.** Most datasets give AC and AN. A 95% Poisson/binomial upper
   bound ("FAF95") per dataset and population would make small cohorts usable for BS1/BA1 without over-trusting them, and
   would replace the fixed AN ≥ 100 rule.
4. **Carrier frequencies for recessive disorders.** Sum the frequencies of ClinVar P/LP (and high-confidence LoF) alleles
   per gene and population to estimate carrier rates and expected affected births, especially for populations with high
   consanguinity. Output: population-specific carrier-screening gene lists.
5. **Founder and population-enriched pathogenic variants.** ClinVar P/LP variants at least 10-fold more frequent in one
   population than in gnomAD's popmax, supported by two datasets, focusing on ACMG secondary-findings genes and newborn
   screening genes.

## Medium value

6. **Pharmacogenomics.** Frequencies of CPIC/PharmGKB actionable alleles (CYP2C19, CYP2D6-tagging SNVs, DPYD, TPMT,
   NUDT15, SLCO1B1, HLA-tagging SNVs) per population, against gnomAD, to show where gnomAD-based dosing assumptions fail.
7. **A curated "do not trust the frequency" resource.** Combine the shared gene outliers from the gene step, the
   single-dataset ClinVar artefacts, and known difficult regions (segmental duplications, GIAB stratifications) into a
   BED and gene list for variant-interpretation pipelines; validate against GIAB truth sets.
8. **Homozygote evidence.** Many datasets give homozygote counts (nhomalt, AC_Hom, HomC). Healthy homozygotes for a
   variant in a recessive gene are evidence against pathogenicity (ACMG BS2); compare homozygote counts across populations.
9. **Non-coding functional variants.** Repeat the "common locally, rare in gnomAD" analysis for splice-altering
   (SpliceAI), UTR and promoter variants, where gnomAD's genome set is smallest.
10. **Ancestry-informative markers for the region.** Select markers from step 3 to build a compact panel distinguishing
    Middle Eastern sub-populations, and test it on HGDP and in-house samples.

## Lower priority or out of scope

- Population-specific gene constraint (observed/expected with a mutation-rate model): needs genotype-level data and
  coverage per dataset to be reliable.
- Structural variants (gnomAD-SV against long-read and national SV catalogues).
- Cross-build (GRCh37/GRCh38) liftover error rates, using the hg19-lifted tables (gnomAD v2, ExAC, Kaviar) as test cases.
