# Data-quality audit

Problems found in public frequency resources and in intermediate processing during this study. Each would have changed
the conclusions if left uncorrected.

## Errors in public resources

| Source | Problem | Evidence | Handling |
|---|---|---|---|
| UCSC ALFA VCF | Indel positions one base too far right | On chr22, REF matches hg38 at POS−1 for 100% of indels and at POS for 5%; 99% of ALFA indels were absent from gnomAD before correction; shifting −1 raised the match to 79% (deletions) and 57% (insertions) | Indels shifted −1 in the union; left out of the frequency analyses (SNVs kept) |
| UCSC GenomeAsia (`varFreqsAll.bb`) | Unreliable frequencies | REF matches hg38, but ~25% of variants at AF 10–90% are absent from both TOPMed and gnomAD | Excluded from the frequency analyses |
| UCSC SPARK, ToMMo (`varFreqsAll.bb`) | Records with tiny allele numbers | e.g. AC = 3 with AF = 1; AC = 1 with AF = 0.5 (AN = 2). These made 960 rare pathogenic variants look common (SPARK) and pushed ToMMo's share of common SNVs absent from gnomAD to 14% (all regions, before the rule; 0.3% exonic and 7.2% non-exonic after it) | AF ignored when AN (or AC/AF) < 100 |
| UCSC `varFreqsAll.bb` | Column order differs from its own schema | Rows with a single source cohort fill another cohort's columns (e.g. SPARK values under `indAC`) | Header rebuilt from the data, verified on 20 M rows |
| UCSC `varFreqsAll.bb` | IndiGen listed in the schema but never a source | Its columns are always empty | Not used |
| gnomAD v4.1 (ANNOVAR tables) | Rows with AF = 0 | Variant seen only in QC-failed genotypes; 60–90% of "common locally, rare in gnomAD" calls were these rows | Treated as no AF; genome AF used in exonic regions |
| gnomAD v4.1 exomes (ANNOVAR table) | Missing group AFs | 12–26% of common sites have no AF for a given group (not assessed in part of the exome data) | Group AF taken from gnomAD genomes |
| ANNOVAR hg38 1000 Genomes | Encoded indel alleles | `0C`, `1AC`, `21` = `<ref length><alt>`; plain SNVs have no prefix | Decoded |
| GME (ANNOVAR) | Indels | 21.6% of GME indel polymorphisms absent from gnomAD, 35.7% rare in every group; several implausible pathogenic indels at ~50% (BRCA2, ASPM) | Reported with a caveat |
| Single-dataset clinical frequencies | Implausible AFs for known pathogenic variants | ALFA MYH7 74%, SVatalog CFTR F508del 77%, GME BRCA2 deletion 47%, UK Biobank ABCC6 deletion 99.99% | ClinVar flags require ≥ 2 independent datasets |

## Errors in the analysis itself, found and fixed

| Problem | Effect before the fix | Fix |
|---|---|---|
| ANNOVAR `.idx` bins are not listed in file order | Chromosome slices read half of gnomAD, All of Us and Kaviar; common SNVs looked 26% "absent from gnomAD" in every dataset | Slice from the smallest start to the largest end offset; verified against full-file counts |
| Insertions written with `-` were shifted one base in an early overlap script | `-` and anchored forms of the same insertion counted twice | Anchor rule fixed before the union was built |
| Allele de-duplication by position only before VEP | ~34 M variants left unannotated | De-duplication on the whole allele; a check reports unannotated variants |
| A file glob picked up the list of SNVs missing from the VEP lookup twice | Those SNVs counted twice per gene | Distinct file names per intermediate list |
