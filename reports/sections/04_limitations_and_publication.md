# Limitations and publication

## Limitations

- **Summary data only.** Allele frequencies, not genotypes: no per-person burden, linkage, relatedness or haplotype
  analysis.
- **Heterogeneous technologies.** Imputed arrays, exomes, short- and long-read genomes, each with its own QC and calling.
  Some differences between datasets reflect technology rather than population (for example long-read cohorts in hard
  regions, imputed cohorts for indels).
- **Third-party conversions.** Several datasets came through ANNOVAR liftovers from hg19 or UCSC conversions; the errors
  found are documented in the [data-quality audit](03_data_quality_audit.md), and others may remain.
- **gnomAD exome and genome not combined by allele number.** The ANNOVAR tables carry AF only, so a joint AF could not be
  computed; the exonic/non-exonic split avoids mixing them.
- **Small cohorts.** Saudi (~300), GoNL (498), GA4K (~550) and SVatalog (101) give noisy frequencies.
- **Region coverage.** No Iranian cohort was included; GME covers the Middle East only partly and has indel quality
  problems.
- **Gene level.** Without a mutation-rate model the comparison is relative to gnomAD, not to an absolute expectation.
- **Additional analyses.** Exonic variants only; ClinVar is biased towards well-studied populations, so founder and
  carrier estimates for other populations are lower bounds; carrier frequencies assume random mating; the ACMG
  comparison uses three genomes; the ancestry panel is untested on genotypes.

## Publication assessment

The study can be presented as a conference abstract now and developed into a journal paper with additions.

**Strengths**

- Scale: 27 datasets, 2.06 billion variants, exonic and non-exonic analysed separately, both directions of disagreement.
- A practical, actionable answer for diagnostic labs: filter on the gnomAD population maximum; add GME for Middle Eastern
  patients; add genome cohorts outside exons.
- A systematic audit of errors in widely used public frequency tracks.
- A clinical angle (ACMG BA1/BS1/PM2) with a replication rule for frequency evidence.

**Needed for a journal paper**

1. A sharper novelty: comparisons of single national databases with gnomAD exist. The strongest angles here are the
   Middle East gap, the error audit, and the popmax rule.
2. An Iranian or Middle Eastern cohort (e.g. Iranome or in-house exomes) to make the Middle East result the centrepiece.
3. Quantified clinical impact on a patient cohort. A test on three genomes is done (23.6–39.7% of VUS resolved,
   [Additional analyses](07_additional_analyses.md#6-acmg-reclassification-in-real-genomes)); a larger cohort is needed.
4. Sensitivity analyses: frequency thresholds, minimum allele number, with and without imputed cohorts.

**Possible venues**

| Venue type | Examples | Status |
|---|---|---|
| Conference | ASHG or ESHG annual meeting (abstract) | Ready with current results |
| Human genetics journal | European Journal of Human Genetics, Human Genetics, Human Genomics, Genetics in Medicine Open | After the additions above |
| Resource or methods journal | NAR Genomics and Bioinformatics, BMC Genomics, Database | If the harmonised union and quality flags are released |

## Data use

No input data are redistributed in this repository. Before releasing any derived per-variant file (for example the
union table), check the terms of each source: gnomAD (open), UK Biobank, All of Us (data and statistics dissemination
policy), FinnGen, ToMMo, TOPMed/BRAVO, SPARK/SFARI, and ANNOVAR's database licence (free for non-commercial use).
