# Materials and methods

## 1. Questions

1. **Added value beyond gnomAD.** For each dataset, what share of its SNVs and indels is absent from gnomAD v4.1, for all
   variants and for common variants (AF ≥ 1%)?
2. **Classification by gnomAD.** How does gnomAD classify a dataset's polymorphisms (AF > 10%) and rare variants
   (AF < 5%): rare (< 5%), low-frequency (5–10%), polymorphic (> 10%), absent or QC-failed?
3. **Markers and structure.** Which SNVs are common in a country dataset but rare or absent in every gnomAD group, and do
   population frequencies reproduce known population structure?
4. **Clinical relevance.** Which ClinVar variants are common in some population but rare in gnomAD, and would therefore be
   handled wrongly by gnomAD-based ACMG frequency criteria (BA1, BS1, PM2)?
5. **Gene level.** Which genes carry more or fewer variants than expected, compared with gnomAD across populations and
   with other genes within a dataset, and is a deviation specific to one dataset or population?

Every question was answered separately for exonic and non-exonic regions, and for SNVs and indels.

## 2. Data

### 2.1 Population datasets

People per dataset are the largest allele number in the file divided by two where the file gives allele numbers,
otherwise the size stated by the source publication (marked *pub.*). Technology is as described by each source.
All data were downloaded between 19 and 22 November 2026.

| Dataset | Population | People | Technology | Source file |
|---|---|---:|---|---|
| gnomAD v4.1 exomes (reference) | Global | 730,947 pub. | Short-read exome | ANNOVAR `hg38_gnomad41_exome` |
| gnomAD v4.1 genomes (reference) | Global | 76,215 pub. | Short-read genome | ANNOVAR `hg38_gnomad41_genome` |
| gnomAD v2 exomes | Global | 123,136 pub. | Short-read exome, lifted from hg19 | ANNOVAR `hg38_gnomad_exome` |
| gnomAD v2 genomes | Global | 15,496 pub. | Short-read genome, lifted | ANNOVAR `hg38_gnomad_genome` |
| ExAC 0.3 | Global | 60,706 pub. | Short-read exome, lifted | ANNOVAR `hg38_exac03` |
| Kaviar 2015-09-23 | 35 aggregated projects | ~77,300 | Mixed, lifted from hg19 | ANNOVAR `hg38_kaviar_20150923` |
| 1000 Genomes phase 3 | 26 populations | 2,504 pub. | Low-coverage genome + exome | ANNOVAR `hg38_ALL.sites.2015_08` |
| All of Us | USA, diverse | not in file | Short-read genome | ANNOVAR `hg38_allofus` |
| GME Variome | Middle East and North Africa | ~2,500 pub. | Short-read exome | ANNOVAR `hg38_gme` |
| TOPMed Freeze 10 | USA, diverse | ~150,900 | Deep short-read genome | UCSC `varFreqsAll.bb` (`top`) |
| ToMMo 61KJPN | Japan | ~61,500 | Short-read genome | UCSC `varFreqsAll.bb` (`tom`) |
| KOVA v7 | Korea | ~5,300 | Short-read sequencing | UCSC `varFreqsAll.bb` (`kov`) |
| MCPS | Mexico City | ~138,200 | Exome sequencing | UCSC `varFreqsAll.bb` (`mcp`) |
| FinnGen R12 | Finland | ~530,300 | Genotyping array + imputation | UCSC `varFreqsAll.bb` (`fin`) |
| SweGen | Sweden | 1,000 | Short-read genome | UCSC `varFreqsAll.bb` (`swe`) |
| Saudi | Saudi Arabia | ~300 | Short-read sequencing | UCSC `varFreqsAll.bb` (`sau`) |
| SPARK iWES v3 | USA, autism families | ~142,400 | Exome | UCSC `varFreqsAll.bb` (`spk`) |
| SFARI WGS pilot | USA, autism families | 12,519 | Short-read genome | UCSC `varFreqsAll.bb` (`wgs`) |
| ALFA (SNVs only) | dbGaP aggregate | not in file | Mixed | UCSC `alfa/ALFA.vcf.gz` |
| CoLoRSdb | Long-read reference | 1,427 | PacBio HiFi | UCSC `colorsDb/colorsDbSnv.vcf.gz` |
| GA4K | USA, rare disease | ~550 | PacBio HiFi | UCSC `ga4k/ga4kSnv.vcf.gz` |
| GoNL | Netherlands | 498 | Short-read genome | UCSC `gonl/gonl.vcf.gz` |
| GREGoR | USA, rare disease | ~3,300 | Short-read genome/exome | UCSC `gregor/gregor.vcf.gz` |
| HGDP + 1000 Genomes | Worldwide | ~3,400 | High-coverage genome | UCSC `hgdp1kFreq/hgdp1k.freq.vcf.gz` |
| HostSeq | Canada | ~9,200 | Short-read genome | UCSC `hostseq/hostseq.vcf.gz` |
| SG10K Health r5.3.2 | Singapore | ~9,700 | Short-read genome | UCSC `npm/SG10K_Health_r5.3.2.sites.vcf.bgz` |
| SVatalog | Reference panel | 101 | Linked-read genome | UCSC `svatalog/svatalog.vcf.gz` |
| UK Biobank | UK | 361,194 | Genotyping array + imputation | UCSC `ukbb/ukbb.vcf.gz` |
| WBBC | China | 4,480 | Short-read genome | UCSC `wbbc/wbbc.vcf.gz` |

Download locations:

- ANNOVAR tables: `annotate_variation.pl -buildver hg38 -downdb -webfrom annovar <name>` from
  `http://www.openbioinformatics.org/annovar/download/`.
- UCSC: `https://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/` (downloaded over http because the server's CA bundle
  rejected the TLS certificate). `varFreqsAll.bb`: 30,276,835,043 bytes, last modified 2026-01-20, 1,259,969,871 records.

### 2.2 Datasets left out of the frequency analyses

| Dataset | Reason |
|---|---|
| GenomeAsia 100K (UCSC `.bb`) | About 25% of its variants at AF 10–90% are absent from both TOPMed and gnomAD; frequencies judged unreliable |
| ALFA indels | Indel positions in the UCSC file are one base too far right; half remain unmatched after correction |
| IndiGen | No allele frequency in the file |
| SGDP | Allele number 2 per site |
| SCHEMA | Case/control counts per study only |
| UCSC All of Us and ABraOM | All of Us taken from ANNOVAR instead |
| Pan-UK Biobank | GRCh37 only |

GenomeAsia, ALFA indels (position-corrected), IndiGen, SGDP and SCHEMA are present in the union table (presence only).

### 2.3 Reference resources

| Resource | Version | Use |
|---|---|---|
| gnomAD v4 exome calling intervals | `ukb.pad50.broad.pad50.union.intervals` (UK Biobank ∪ Broad targets, ±50 bp; 393,320 intervals, 170.4 Mb) from `gs://gcp-public-data--gnomad/resources/grch38/intervals/` | Exonic / non-exonic split |
| ClinVar | `clinvar_20261115.vcf.gz` (GRCh38, MD5 verified) from `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/` | Step 4 |
| Ensembl VEP | 116.3, merged (Ensembl + RefSeq) cache 116, GRCh38 | Gene and consequence |
| VEP annotation of every possible exonic SNV | Existing run on the CADD v1.7 exon SNV set (exons ±50 bp) | Gene and consequence for 106.8 M SNVs |
| hg38 reference | `Homo_sapiens_assembly38.fasta` / UCSC `hg38.fa` | VEP; left-alignment checks |
| VEP 116 RefSeq cache (`NASER_FILES/Naser_VEP/vep_cache`) | gnomAD v4.1 exome and genome AFs per group | Additional analyses (gnomAD group AFs for ClinVar P/LP; ACMG annotation) |
| ClinGen gene–disease validity | Release in the BIAS-2015 data (2026-10-15) | Autosomal recessive genes for carrier frequencies |
| GIAB genome stratifications v3.5, GRCh38 | `GRCh38_segdups.bed.gz`, `GRCh38_lowmappabilityall.bed.gz` | Validation of unreliable genes |
| BIAS-2015 | Local installation with its hg38 data (2026-10-15), REVEL, AlphaMissense | ACMG impact |

## 3. Processing

### 3.1 Preparing the inputs

- `varFreqsAll.bb` was converted to a bgzipped, tabix-indexed table with UCSC `bigBedToBed` (one job per chromosome).
  The row count matched the file's own count. The bigBed's column order differs from its schema (values are grouped per
  cohort: AC, AF, then that cohort's population AC/AF); the header was rebuilt from the data and verified on 20 million
  rows (`fix_bb_header.sh`).
- Of the 13 cohorts in the bigBed, All of Us and ABraOM were ignored (All of Us comes from ANNOVAR). The other 11 cohorts
  were read from its `sources` column in one pass.
- Only the 13 UCSC cohorts not contained in the bigBed were downloaded as VCFs.

### 3.2 Harmonisation

- Notation: every allele in ANNOVAR form. `chr` prefix removed, `MT` → `M`, shared leading then shared trailing bases
  trimmed, empty allele written `-`, insertions anchored on the preceding base (Start = End = that base).
- Multi-allelic VCF records split; symbolic and `*` alleles dropped; alleles with characters other than A, C, G, T, N
  dropped (counted in `results/01_union/extraction_stats.tsv`).
- Only chromosomes 1–22, X, Y and M kept.
- The hg38 1000 Genomes tables write indel alternate alleles as `<ref length><alt>` (`0C`, `1AC`, `21`); these were decoded.
- ALFA indels: position − 1 (see the [data-quality audit](03_data_quality_audit.md)); in the frequency analyses only ALFA
  SNVs were used.
- Inputs were taken as left-aligned. A check with `bcftools norm -f hg38` on chromosome 22 realigned fewer than 50
  records per source.

### 3.3 Quality rules

- Alleles with AC = 0 dropped; VCF records with FILTER other than `PASS` or `.` dropped.
- A dataset's AF was set to missing when its allele number was below 100 (AN from the file, Kaviar's AN column, or AC/AF).
- gnomAD AF = 0 (carriers only among genotypes failing gnomAD QC) treated as "no AF".
- gnomAD exome group AFs missing for 12–26% of common sites (groups not assessed in part of the exome data) taken from
  gnomAD genomes.

### 3.4 Union

All 33 inputs were merged per chromosome (sorted keys, bitmask of datasets) into `hg38_freq_union.txt`
(`#Chr Start End Ref Alt is_in how_many`), sorted 1–22, X, Y, M, with an ANNOVAR `.idx` index.

### 3.5 Regions and gnomAD reference

| Region | Definition | gnomAD reference AF | Datasets |
|---|---|---|---|
| Exonic | Start inside the gnomAD v4 exome calling intervals | Exome AF; genome AF where the exome table has no AF | All |
| Non-exonic | Everything else | Genome AF only (exome-table calls ignored) | All except exome-only sets (SPARK WES, GME, ExAC, gnomAD v2 exomes) |

Three gnomAD views:

- **overall**: the gnomAD AF.
- **popmax**: the highest of overall AF and the AFs of AFR, AMR, ASJ, EAS, FIN, MID, NFE, SAS. Computed here because
  gnomAD's own `grpmax` excludes MID and FIN.
- **matched**: the AF of the dataset's own gnomAD group: EAS for ToMMo, KOVA, WBBC, SG10K; AMR for MCPS; FIN for
  FinnGen; NFE for SweGen, GoNL, UK Biobank; MID for Saudi and GME.

### 3.6 Classes

| Variable | Classes |
|---|---|
| Dataset AF | < 1%, 1–5%, 5–10% (low), > 10% (polymorphism), no AF |
| gnomAD | absent, QC fail (row without usable AF), rare (< 5%), low (5–10%), polymorphic (> 10%) |
| Variant type | SNV; indel (one allele `-`); MNV (other) |

Rare = < 5% in the dataset; polymorphism = > 10%. In gnomAD, rare = < 5% in exome and genome; polymorphic = > 10% in
exome or genome (one reference AF per region, section 3.5).

## 4. Analyses

| Step | Analysis | Definition |
|---|---|---|
| 1 | Unique to the dataset | Share of the dataset's variants absent from gnomAD, all and AF ≥ 1% |
| 2 | Polymorphisms and rare variants | Cross-tabulation of dataset class × gnomAD class, for each gnomAD view, in both directions |
| 3 | Markers | SNVs > 10% in a country dataset with gnomAD popmax < 5% or absent |
| 3 | Structure | PCA (SVD of the centred population × SNV AF matrix) and mean absolute AF difference; SNVs with gnomAD genome AF, AF ≥ 1% in some country set, thinned to one per 100 kb; country datasets with AF for at least 20% of these SNVs, plus the eight gnomAD genome groups |
| 4 | ClinVar | ClinVar alleles (normalised as in 3.2) with AF ≥ 5% in at least one dataset; grouped as P/LP, VUS, conflicting, B/LB, other; review stars from CLNREVSTAT; support = number of independent datasets with AF ≥ 5% (gnomAD v2 and ExAC counted with gnomAD) |
| 5 | Gene level | See below |
| 6–10 | Additional analyses | Founder variants, carrier frequencies, unreliable genes and regions, ancestry markers, ACMG impact: methods in [Additional analyses](07_additional_analyses.md) |

### 4.1 Gene level (exonic, protein-coding genes)

- **Annotation.** 195,874,106 exonic variants. 106,845,209 SNVs took gene and consequence from the existing VEP
  annotation of every possible exonic SNV; the other 89,022,055 (indels, MNVs and SNVs outside that set) were run
  through VEP 116.3 in its default input format. One consequence per allele (`--pick_allele`, order MANE Select,
  MANE Plus Clinical, canonical, APPRIS, TSL, biotype, CCDS, rank, length).
- **Classes.** Protein-altering = missense, stop gained/lost, start lost, frameshift, in-frame insertion/deletion,
  splice acceptor/donor, protein-altering; synonymous as a technical control. Each for all variants and for AF ≥ 1%.
- **Across populations.** For gene *g* in dataset *d*: expected = N<sub>d</sub> × (n<sub>gnomAD</sub>(g) + 0.5) /
  (N<sub>gnomAD</sub> + 0.5 × genes). Gene length and dataset size cancel. Two-sided exact binomial test, tested when
  observed or expected ≥ 5; Benjamini–Hochberg FDR per dataset, class and AF class. A gene is high or low at FDR < 5%
  and |log<sub>2</sub>(observed/expected)| ≥ 1.
- **Within a dataset.** Variants per kb of the gene's exonic footprint (positions VEP assigns to the gene in exons
  ±50 bp, ≥ 300 bp); robust z-score of log<sub>10</sub> density (median, 1.4826 × MAD); outlier at |z| ≥ 3.
- **Uniqueness.** For each flagged gene: number of datasets and population groups (East Asia, Europe, Middle East,
  Latin America, multi-population) in which it is flagged.

## 5. Checks performed

| Check | Result |
|---|---|
| Union and study row counts against source counts | Match (e.g. bigBed 1,259,969,871 rows) |
| Indel placement per source: REF against hg38 at POS and POS−1; overlap with gnomAD when shifted ±1 | Only ALFA misplaced |
| Normalisation with `bcftools norm -f hg38` (chr22) | < 50 realignments per source |
| Chromosome slices read from ANNOVAR `.idx` against full-file counts | Match after taking the minimum/maximum byte offsets |
| Cross-dataset consistency of rare-allele frequencies | Revealed tiny-AN records (SPARK, ToMMo) |
| ClinVar flags against independent datasets | 127 of 140 P/LP flags rest on one dataset |
| Gene step: every variant annotated once | 106,845,209 + 89,022,055 = 195,874,106 |

## 6. Software and compute

| Software | Version |
|---|---|
| GNU awk | 5.1.1 |
| GNU coreutils (sort) | 9.0 |
| bcftools / tabix (htslib) | 1.19 |
| UCSC bigBedToBed, bigBedInfo (bioconda) | 482 |
| aria2 | 1.37.0 |
| Ensembl VEP (Docker) | 116.3 (Ensembl API 116) |
| Python | 3.10.22 |
| pandas / NumPy / SciPy / matplotlib | 2.2.2 / 1.26.4 / 1.14.1 / 3.10.9 |

Server: 192-core ARM (Kunpeng), 1.5 TB RAM. Extraction runs as 511 jobs (120 in parallel). A full study run takes about
1.7 hours; the gene step about 3 hours, of which VEP takes 2 hours on 120 cores.
