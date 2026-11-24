# Reproducibility protocol

How to regenerate every result in this repository, either from the kept archive (no raw files needed) or from scratch
(re-downloading the raw frequency files). Everything runs through one driver:

```bash
bash scripts/run_all.sh <stage> [<stage> ...]
```

Paths and settings are in [`scripts/config.sh`](../../scripts/config.sh); each stage logs to `$DATA/logs/<stage>.log`.

## 1. Two ways to regenerate

| Route | Needs | Stages | Time (192-core server) | Regenerates |
|---|---|---|---|---|
| **From the kept archive** | `$KEEP` (≈ 21 GB), the reference resources in section 3 | `all-kept` = `extra reports` (+ `acmg`) | ≈ 3 h (+ 20 min) | All tables, figures and documents; all additional analyses |
| **From scratch** | Internet access, ≈ 600 GB free disk, section 3 | `all-raw` = `download prepare union study clinvar gene keep extra reports` | ≈ 12 h (downloads excluded) | Everything, including the union, the harmonised per-dataset frequencies and the kept archive |

The kept archive (`/home/NASER/annotation/population_study/data/archive/`) is what remains after the raw frequency files are
deleted. Its content and checksums are listed in [Data manifest](06_data_manifest.md) and `results/00_provenance/kept_archive.md5`.

## 2. Software

| Software | Version used | Where |
|---|---|---|
| bash, GNU awk, GNU coreutils (sort) | bash 5, gawk 5.1.1, coreutils 9.0 | system |
| bcftools, tabix, bgzip (htslib) | 1.19 | system (`/usr/local/bin`) |
| bedtools | 2.x | system |
| UCSC bigBedToBed, bigBedInfo; aria2 | 482; 1.37.0 | conda env `ucsc` (bioconda, conda-forge) |
| Python with pandas, NumPy, SciPy, matplotlib | 3.10.22; 2.2.2, 1.26.4, 1.14.1, 3.10.9 | conda env `ahm-env` (`$PY`) |
| Ensembl VEP | 116.3 (API 116), Docker image `huvya/vep:116`, `sha256:b8901bb319b3…` | Docker |
| ANNOVAR `annotate_variation.pl` | 2025-03-02 | `$ANNOVAR` (needed only for `download`) |
| BIAS-2015 | git commit `290b99a`, plus `scripts/analyses/bias_pm4_length_chrM.patch` | `$BIAS` (read-only) |

To rebuild the conda environments elsewhere:

```bash
conda create -n ucsc -c conda-forge -c bioconda ucsc-bigbedtobed=482 ucsc-bigbedinfo=482 aria2
conda create -n popfreq -c conda-forge python=3.10 pandas=2.2 numpy=1.26 scipy=1.14 matplotlib=3.10
```

## 3. Reference resources (not deleted with the raw files)

| Resource | Path (config variable) | Used by |
|---|---|---|
| VEP 116 RefSeq cache with gnomAD v4.1 AFs | `/home/NASER_FILES/Naser_VEP/vep_cache` (`VEP_CACHE`) | gene step (indels), ClinVar P/LP group AFs, ACMG |
| VEP annotation of every possible exonic SNV (CADD v1.7 exon set) | `$VEP_ROOT/results/CADD_exons_merged_64cpu/…vep.tsv.gz` (`CADD_VEP_TSV`) | gene step (SNV lookup) |
| hg38 FASTA | `/home/NASER_FILES/annovar/hg38/reference/hg38.fa` (`FA`) | VEP |
| BIAS-2015 data (hg38, 2026-10-15), REVEL, AlphaMissense, ClinVar 2026-10-11 | `/home/NASER/annotation/acmg` (`ACMG`) | ACMG, ClinGen AR genes |
| BIAS-2015 code | `/root/naser/annotation/acmg/BIAS-2015` (`BIAS`) | ACMG |
| Patient genome | `$VEP_ROOT/results/HUV00001_refseq_32cpu/HUV00001.norm.vcf.gz` (`VCF`) | ACMG |
| GIAB v3.5 GRCh38 stratifications | downloaded by the `extra` stage | unreliable regions |

Nothing under `/home/HUVYA` is used. The published gene step used the HUVYA merged VEP cache (read-only); the default
is now the `Naser_VEP` RefSeq cache, which can pick a different transcript for a small share of variants.

## 4. Raw inputs (the `download` stage)

| Source | What | Exact identity |
|---|---|---|
| ANNOVAR | `gnomad41_exome`, `gnomad41_genome`, `gnomad_exome`, `gnomad_genome`, `kaviar_20150923`, `exac03`, `1000g2015aug`, `gme`, `allofus` (hg38) | file sizes and dates in [Data manifest](06_data_manifest.md); MD5 in `results/00_provenance/raw_inputs.md5` |
| UCSC varFreqs | `varFreqsAll.bb` (30,276,835,043 bytes, modified 2026-01-20) and 13 cohort VCFs | `http://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/`; MD5 in `results/00_provenance/raw_inputs.md5` |
| gnomAD | v4 exome calling intervals `ukb.pad50.broad.pad50.union.intervals` | `GNOMAD_INTERVALS_URL` in `config.sh` |
| ClinVar | `clinvar_20261115.vcf.gz` (GRCh38) | NCBI archive URL in `config.sh`; MD5 checked on download |

Public resources change: UCSC and ANNOVAR replace files in place. After `download`, the driver compares every file
with `raw_inputs.md5` and reports whether the inputs are identical to the published run. If they differ, the results
can differ too; keep the kept archive as the reference.

## 5. Stages

| Stage | Script(s) | Output | Time | Check |
|---|---|---|---|---|
| `download` | ANNOVAR `-downdb`, `get_ucsc_missing.sh`, `curl` | `$ANNOVAR_TABLES`, `$UCSC_DIR`, `$CLINVAR`, `$INTERVALS` | hours (network) | log says "raw inputs identical to the published run" |
| `prepare` | `bb_to_table.sh`, `fix_bb_header.sh` | `varFreqsAll.tsv.gz` (+ `.tbi`) | ≈ 1.5 h | 1,259,969,871 rows = `bigBedInfo` itemCount |
| `union` | `freq_union.sh` | `hg38_freq_union.txt` (+ `.idx`, summaries) | ≈ 5 h | 2,055,347,367 variants |
| `study` | `freq_study.sh` | `exonic/`, `non_exonic/` counts and per-variant files | ≈ 1.7 h | 293,015,270 exonic and 3,917,707,046 non-exonic dataset alleles classified |
| `clinvar` | `freq_clinvar.sh`, `freq_study_report.py` | ClinVar matches, tables, figures | ≈ 20 min | — |
| `gene` | `freq_study.sh` (`DUMP=1`), `snv_annot.sh`, `gene_burden.sh`, `gene_burden_report.py` | per-gene counts, gene tables and figures | ≈ 3.5 h | 106,845,209 SNVs by lookup + 89,022,055 by VEP = 195,874,106; no "without VEP annotation" warning |
| `keep` | `make_keep.sh` | `$KEEP` + `provenance/*.md5` | ≈ 1 h | `KEEP_DONE` in the log |
| `extra` | `extra/*` (except ACMG) | `$EXTRA_DIR/<analysis>/{tables,figures}` | ≈ 70 min | numbers in section 6 |
| `acmg` (optional) | `extra/acmg_impact.sh`, `extra/acmg_compare.py` | `$EXTRA_DIR/acmg_impact` | ≈ 20 min for 3 genomes | numbers in section 6 |
| `reports` | `freq_study_report.py`, `gene_burden_report.py` on `$KEEP`; copy; `make_results_md.py`; `build_report_md.sh` | `results/`, `reports/sections/02_results.md`, `reports/population_frequency_report.md` | ≈ 1 h 50 min (the gene-level binomial tests take most of it) | `git diff --stat results/` shows no numeric change |

## 6. Reference values to check a re-run against

| Quantity | Value |
|---|---:|
| Union variants | 2,055,347,367 |
| ToMMo common exonic SNVs absent from gnomAD (AF ≥ 1%) | 0.31% |
| GME common exonic SNVs absent from gnomAD | 7.46% (7,363) |
| TOPMed common non-exonic SNVs absent from gnomAD | 7.65% |
| Exonic ClinVar P/LP flags: total / supported by ≥ 2 datasets | 140 / 12 |
| Founder P/LP variants: enriched / replicated / high confidence | 1,621 / 61 / 32 |
| Carrier frequency CFTR in gnomAD NFE | 4.4% |
| Unreliable genes / with ≥ 2 rules | 57 / 11 |
| ACMG (three genomes): VUS resolved, B1 / B2 | 6,954 / 4,138 of 17,516 (39.7% / 23.6%) |
| PGR comparison: pLoF / synonymous ratio, All of Us | 1.93 |

## 7. Pitfalls already handled (do not undo)

| Pitfall | Where handled |
|---|---|
| ANNOVAR `.idx` bins are not in file order: slice from the smallest to the largest byte offset | `freq_study.sh` (`chunk`) |
| hg38 1000 Genomes indel alleles are encoded `<ref length><alt>` | `freq_union.sh`, `freq_study.sh` |
| UCSC ALFA indel positions are one base too far right | union: shift −1; study: ALFA SNVs only |
| `varFreqsAll.bb` column order differs from its schema | `fix_bb_header.sh` |
| GenomeAsia frequencies unreliable | excluded in `freq_study.sh` |
| Tiny allele numbers (AC = 3, AF = 1) in SPARK, ToMMo | `MINAN=100` |
| gnomAD AF = 0 means QC-failed carriers only | treated as no AF |
| gnomAD exome group AFs missing at 12–26% of common sites | genome fallback per group |
| VEP RefSeq cache is BAM-edited: needs `--fasta` | `config.sh` `FA` |
| VEP gnomAD field names are `gnomADe_<GROUP>_AF` | `prep_clinvar_plp.sh` |
| De-duplicate VEP input on the whole allele, not on position | `gene_burden.sh` |
| Single-dataset frequencies are mostly errors | replication rules in ClinVar, founder and carrier analyses |
| CFTR poly-T/TG intronic alleles inflate carrier rates | coding P/LP only in carrier analysis |

## 8. Verification done (24 November 2026)

`bash scripts/run_all.sh all-kept` was run from the kept archive on a clean copy of the repository. Every result table
under `results/` (study, gene level and additional analyses) was regenerated byte-for-byte identical; only PDF files
differ, through their embedded creation timestamps. The ACMG stage was then run separately on three genomes.

## 9. Before deleting the raw files

1. Check `$KEEP/provenance/raw_inputs.md5` and `kept_archive.md5` exist and are non-empty (they are copied to `results/`).
2. Run `bash scripts/run_all.sh all-kept` once and confirm that `git diff --stat results/` shows no numeric change.
3. Remember that `/root/naser/annotation/annovar/humandb/` holds symbolic links into `$ANNOVAR_TABLES`.
