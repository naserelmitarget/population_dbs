# Configuration for scripts/run_all.sh. Every path the pipeline uses is set here; edit, then run a stage.
# Values are exported, so each script picks them up instead of its own default.

# ---- where data live ---------------------------------------------------------------------------------------
# population_study/data: raw/ (deletable inputs), work/ (intermediates), archive/ (kept), logs/
export STUDY_ROOT=${STUDY_ROOT:-/home/NASER/annotation/population_study}
export DATA=${DATA:-$STUDY_ROOT/data}
export ANNOVAR_TABLES=${ANNOVAR_TABLES:-$DATA/raw/annovar_popfreq}   # ANNOVAR frequency tables (raw input)
export UCSC_DIR=${UCSC_DIR:-$DATA/raw/ucsc_varfreqs}            # UCSC varFreqs files (raw input)
export CLINVAR_DIR=${CLINVAR_DIR:-$DATA/raw/clinvar}
export CLINVAR=${CLINVAR:-$CLINVAR_DIR/clinvar_20261115.vcf.gz}
export UNION_DIR=${UNION_DIR:-$DATA/work/freq_union}
export STUDY_DIR=${STUDY_DIR:-$DATA/work/freq_study}
export GENE_DIR=${GENE_DIR:-$DATA/work/gene_burden}
export KEEP=${KEEP:-$DATA/archive}                          # kept archive: enough to regenerate without raw files
export EXTRA_DIR=${EXTRA_DIR:-$DATA/work/extra_analyses}
export INTERVALS=${INTERVALS:-$STUDY_DIR/regions/gnomad_v4_exome_calling.ukb_broad_pad50.intervals}

# ---- tools and references (read-only) -------------------------------------------------------------------
export ANNOVAR=${ANNOVAR:-/root/naser/annotation/annovar/annotate_variation.pl}
export UCSC_BIN=${UCSC_BIN:-/root/miniforge3/envs/ucsc/bin}          # bigBedToBed, bigBedInfo, aria2c (bioconda)
export PY=${PY:-/root/miniforge3/envs/ahm-env/bin/python}            # pandas, numpy, scipy, matplotlib
export VEP_IMAGE=${VEP_IMAGE:-huvya/vep:116}                         # Ensembl VEP 116.3 (sha256:b8901bb3...)
export VEP_ROOT=${VEP_ROOT:-/home/NASER_FILES/Naser_VEP}             # VEP RefSeq 116 cache + CADD-exon VEP table
export VEP_CACHE=${VEP_CACHE:-$VEP_ROOT/vep_cache}
export VEP_CACHE_OPT=${VEP_CACHE_OPT:---refseq}
export VEP_MOUNT=${VEP_MOUNT:-/home/NASER_FILES}
export FA=${FA:-/home/NASER_FILES/annovar/hg38/reference/hg38.fa}
export VEP_FASTA=$FA
export CADD_VEP_TSV=${CADD_VEP_TSV:-$VEP_ROOT/results/CADD_exons_merged_64cpu/hg38_cadd17_exons_pad50.vep.tsv.gz}
export ACMG=${ACMG:-/home/NASER/annotation/acmg}                                      # BIAS-2015 data, REVEL, AlphaMissense
export BIAS=${BIAS:-/root/naser/annotation/acmg/BIAS-2015}           # commit 290b99a + scripts/analyses/bias_pm4_length_chrM.patch
export CLINGEN=${CLINGEN:-$ACMG/bias_hg38_current/hg38_clingen_gene_disease_validity.csv}
# genomes for the ACMG test (stage "acmg"): NAME=GRCh38 VCF, space-separated; copies outside /home/HUVYA
export ACMG_SAMPLES=${ACMG_SAMPLES:-"HUV00001=$VEP_ROOT/results/HUV00001_refseq_32cpu/HUV00001.norm.vcf.gz HUV00002=$EXTRA_DIR/acmg_impact/inputs/HUV00002.hg38.vcf.gz HUV00003=$EXTRA_DIR/acmg_impact/inputs/HUV00003.hg38.vcf.gz"}

# ---- compute ------------------------------------------------------------------------------------------------
export JOBS=${JOBS:-120}        # parallel extraction jobs
export CPUS=${CPUS:-120}        # VEP cores in the gene step
export MINAN=${MINAN:-100}      # minimum allele number for a dataset AF

# ---- download sources (as used in November 2026) ---------------------------------------------------------------
export UCSC_URL=http://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs
export GNOMAD_INTERVALS_URL=https://storage.googleapis.com/gcp-public-data--gnomad/resources/grch38/intervals/ukb.pad50.broad.pad50.union.intervals
export CLINVAR_URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20261115.vcf.gz
export ANNOVAR_DBS="gnomad41_exome gnomad41_genome gnomad_exome gnomad_genome kaviar_20150923 exac03 1000g2015aug gme allofus"

# the per-script names some scripts read
export A=$ANNOVAR_TABLES U=$UCSC_DIR BB=$UCSC_DIR/varFreqsAll/varFreqsAll.tsv.gz K=$KEEP G=$GENE_DIR V=$VEP_CACHE V_ROOT=$VEP_ROOT
