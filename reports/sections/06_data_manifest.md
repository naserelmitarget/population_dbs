# Data manifest

What the study read, what it produced, and what is kept after the raw inputs are deleted.

## Raw inputs (deletable; all can be downloaded again)

Sizes and dates as on the analysis server, 24 November 2026.

### ANNOVAR tables (`/home/NASER/annotation/population_study/data/raw/annovar_popfreq/`, 233 GB)

| File | Size | Date | Re-download |
|---|---:|---|---|
| `hg38_gnomad41_genome.txt` (+ `.idx`) | 102.5 GB | 2024-06-02 | `annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad41_genome humandb/` |
| `hg38_allofus.txt` (+ `.idx`) | 85.5 GB | 2025-03-14 | `… -downdb -webfrom annovar allofus` |
| `hg38_gnomad41_exome.txt` (+ `.idx`) | 18.6 GB | 2024-05-31 | `… gnomad41_exome` |
| `hg38_gnomad_genome.txt` (+ `.idx`) | 17.1 GB | 2019-12-06 | `… gnomad_genome` |
| `hg38_kaviar_20150923.txt` (+ `.idx`) | 12.2 GB | 2019-12-06 | `… kaviar_20150923` |
| `hg38_ALL.sites.2015_08.txt` and AFR, AMR, EAS, EUR, SAS | 9.0 GB | 2015-08-26 | `… 1000g2015aug` |
| `hg38_gnomad_exome.txt` (+ `.idx`) | 1.1 GB | 2019-12-06 | `… gnomad_exome` |
| `hg38_exac03.txt` (+ `.idx`) | 0.6 GB | 2019-12-06 | `… exac03` |
| `hg38_gme.txt` (+ `.idx`) | 0.07 GB | 2019-12-06 | `… gme` |

### UCSC varFreqs (`/home/NASER/annotation/population_study/data/raw/ucsc_varfreqs/`, 122 GB)

Base URL: `http://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/` (each VCF with its `.tbi`). Script: `scripts/get_ucsc_missing.sh`.

| File | Size | Notes |
|---|---:|---|
| `varFreqsAll/varFreqsAll.bb` | 30.3 GB | Last modified at UCSC 2026-01-20; 1,259,969,871 records |
| `varFreqsAll/varFreqsAll.tsv.gz` (+ `.tbi`) | 25.2 GB | Derived: `bb_to_table.sh` + `fix_bb_header.sh` |
| `gregor/gregor.vcf.gz` | 20.6 GB | |
| `npm/SG10K_Health_r5.3.2.sites.vcf.bgz` | 17.6 GB | |
| `hgdp1kFreq/hgdp1k.freq.vcf.gz` | 15.9 GB | |
| `hostseq/hostseq.vcf.gz` | 9.5 GB | |
| `alfa/ALFA.vcf.gz` | 2.8 GB | |
| `wbbc/wbbc.vcf.gz` | 2.2 GB | |
| `sgdpFreq/sgdp.freq.vcf.gz` | 2.0 GB | In the union only |
| `gonl/gonl.vcf.gz` | 1.7 GB | |
| `colorsDb/colorsDbSnv.vcf.gz` | 0.7 GB | |
| `ukbb/ukbb.vcf.gz` | 0.55 GB | |
| `schema/SCHEMA_variant_results_withAF.vcf.gz` | 0.54 GB | In the union only |
| `ga4k/ga4kSnv.vcf.gz` | 0.54 GB | |
| `svatalog/svatalog.vcf.gz` | 0.18 GB | |

### Intermediate files (deletable)

| Location | Size | Content |
|---|---:|---|
| `/home/NASER/annotation/population_study/data/work/freq_union/hg38_freq_union.txt` (+ `.idx`) | 102 GB | Uncompressed union; a compressed copy is kept (below) |
| `/home/NASER/annotation/population_study/data/work/gene_burden/{split,vep,annot,work}` | ~104 GB | Gene-step intermediates; the annotated variant list is kept compressed (below) |

**Before deleting `annovar_popfreq/`:** `/root/naser/annotation/annovar/humandb/` holds symbolic links to these tables.
ANNOVAR runs that use `gnomad41_exome`, `gnomad41_genome`, `allofus`, `kaviar_20150923`, `exac03`, `gme`, `gnomad_exome`,
`gnomad_genome` or `1000g2015aug` from that folder will stop working.

## Kept results

### In this repository

`results/` (tables and figures for both regions and the gene level), `results/02_frequency_vs_gnomad/supplementary/` (per-variant exonic markers
and ClinVar matches), `docs/`, `report/`, `scripts/`.

### On the server: `/home/NASER/annotation/population_study/data/archive/`

| Path | Content |
|---|---|
| `union/hg38_freq_union.txt.gz` (+ `.tbi`) | The union: 2,055,347,367 variants with `is_in` and `how_many` (bgzip, tabix-indexed) |
| `union/summary_*.tsv`, `run.log` | Variants per dataset, per number of datasets, per dataset combination |
| `study/{exonic,non_exonic}/` | `counts.tsv`, per-variant `markers.tsv.gz`, `clinvar_candidates.tsv.gz`, `clinvar_hits.tsv.gz`, `pca_snvs.tsv.gz`, all tables and figures |
| `study/regions/` | gnomAD v4 exome calling intervals used for the split |
| `study/extraction_stats.tsv`, `run.log` | Rows kept and dropped per source |
| `gene_level/exonic_variants_annotated.tsv.gz` | Every exonic variant: gnomAD AF, AF in each dataset, VEP consequence, gene, MANE |
| `gene_level/gene_counts.tsv.gz`, `gene_length.tsv` | Per gene × dataset counts and gene footprints |
| `gene_level/tables/` | Including the two full per-gene tables not stored in this repository |
| `clinvar_20261115.vcf.gz` (+ `.tbi`, `.md5`) | ClinVar release used |

### Additional analyses: `/home/NASER/annotation/population_study/data/work/extra_analyses/`

One folder per analysis (`founder_variants`, `carrier_frequency`, `unreliable_regions`, `ancestry_markers`,
`acmg_impact`, `pgr_comparison`), each with `tables/` and `figures/`; `shared/` holds the ClinVar P/LP table with
population and gnomAD group AFs; `acmg_impact/<sample>/` holds each genome's VEP JSON and three BIAS outputs, and
`acmg_impact/inputs/` the copies of the HUV00002 and HUV00003 genome VCFs (copied read-only from the HUVYA WGS output).
Keep `extra_analyses/acmg_impact/inputs/` if you want to rerun the ACMG stage.

With these files every table and figure can be regenerated without the raw inputs (`freq_study_report.py`,
`gene_burden_report.py`, `make_results_md.py`). Re-running harmonisation, extraction or the union needs the raw inputs again.
