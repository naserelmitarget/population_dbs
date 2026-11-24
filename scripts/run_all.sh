#!/bin/bash
# Regenerate the population-frequency study.
#
#   bash scripts/run_all.sh <stage> [<stage> ...]
#
# Stages, in order. Stages 1-7 need the raw frequency files (download them with stage "download");
# stages "extra" and "reports" need only the kept archive ($KEEP) and run after the raw files are deleted.
#
#   download   ANNOVAR tables, UCSC varFreqs (.bb + 13 VCFs), gnomAD exome intervals, ClinVar 2026-11-15
#   prepare    varFreqsAll.bb -> tabix table with corrected header
#   union      hg38_freq_union.txt with is_in / how_many
#   study      extraction, harmonisation, QC, exonic / non-exonic classification
#   clinvar    ClinVar join (both regions)
#   gene       exonic variant dump, VEP annotation, per-gene counts
#   keep       build the kept archive (+ checksums of raw inputs and archive)
#   ---------- from here on the raw files are not needed ----------
#   extra      additional analyses: ClinVar P/LP table, founder, carrier, unreliable genes, ancestry, PGR comparison
#   acmg       ACMG reclassification test on the genomes in ACMG_SAMPLES (optional; needs their VCFs)
#   reports    tables and figures from the kept archive, copied into results/, then reports/sections/02_results.md and the full reports
#
#   all-raw    download prepare union study clinvar gene keep extra acmg reports
#   all-kept   extra reports            (add "acmg" to include the ACMG test: all-kept acmg, or extra acmg reports)
#
# Paths: scripts/config.sh. Logs: $DATA/logs/<stage>.log.
set -uo pipefail
HERE=$(cd $(dirname $(readlink -f $0)) && pwd)
REPO=$(dirname $HERE)
source $HERE/config.sh
mkdir -p $DATA/logs
log() { echo "[$(date +%F' '%T)] $*"; }

stage_download() {
  source /etc/profile.d/npu-proxy.sh 2>/dev/null || true
  mkdir -p $ANNOVAR_TABLES $CLINVAR_DIR $(dirname $INTERVALS)
  for d in $ANNOVAR_DBS; do
    ( cd $(dirname $ANNOVAR_TABLES) && perl $ANNOVAR -buildver hg38 -downdb -webfrom annovar $d $(basename $ANNOVAR_TABLES) )
  done
  bash $HERE/get_ucsc_missing.sh
  curl -sS -o $INTERVALS $GNOMAD_INTERVALS_URL
  for x in "" .tbi .md5; do curl -sS -o $CLINVAR$x $CLINVAR_URL$x; done
  [ "$(cut -d' ' -f1 $CLINVAR.md5)" = "$(md5sum < $CLINVAR | cut -d' ' -f1)" ] && log "ClinVar MD5 OK" || log "ClinVar MD5 MISMATCH"
  [ -s $KEEP/provenance/raw_inputs.md5 ] && ( cd $(dirname $ANNOVAR_TABLES) && md5sum -c --quiet $KEEP/provenance/raw_inputs.md5 ) \
    && log "raw inputs identical to the published run" || log "raw inputs differ from the published run (sources may have been updated)"
}
stage_prepare() { bash $HERE/bb_to_table.sh && bash $HERE/fix_bb_header.sh; }
stage_union()   { W=$UNION_DIR bash $HERE/freq_union.sh; }
stage_study()   { W=$STUDY_DIR IV=$INTERVALS bash $HERE/freq_study.sh; }
stage_clinvar() {
  for r in exonic non_exonic; do
    W=$STUDY_DIR/$r CV=$CLINVAR bash $HERE/freq_clinvar.sh
    $PY $HERE/freq_study_report.py $STUDY_DIR/$r $r
  done
}
stage_gene() {
  DUMP=1 W=$GENE_DIR/work IV=$INTERVALS bash $HERE/freq_study.sh
  bash $HERE/snv_annot.sh
  bash $HERE/gene_burden.sh
  $PY $HERE/gene_burden_report.py $GENE_DIR
}
stage_keep() { bash $HERE/make_keep.sh; }
stage_extra() {
  X=$HERE/extra
  W=$EXTRA_DIR bash $X/prep_clinvar_plp.sh
  W=$EXTRA_DIR bash $X/unreliable_regions.sh                 # before founder/carrier: they use its gene list
  $PY $X/founder_and_carrier.py $EXTRA_DIR
  $PY $X/ancestry_markers.py $EXTRA_DIR
  mkdir -p $EXTRA_DIR/pgr_comparison/tables
  OUT=$EXTRA_DIR/pgr_comparison/tables bash $X/pgr_consequence_absent.sh
  $PY - <<EOF
import pandas as pd
d = pd.read_csv("$EXTRA_DIR/pgr_comparison/tables/coding_variants_absent_from_gnomad_by_consequence.tsv", sep="\t")
p = d.pivot(index="dataset", columns="class", values="pct_absent")[["synonymous", "missense", "inframe", "pLoF"]]
p["coding variants"] = d.pivot(index="dataset", columns="class", values="variants").sum(axis=1)
p["pLoF/syn ratio"] = p.pLoF / p.synonymous
p.to_csv("$EXTRA_DIR/pgr_comparison/tables/coding_absent_pct_by_consequence_wide.tsv", sep="\t", float_format="%.2f")
EOF
}
stage_acmg() {
  W=$EXTRA_DIR/acmg_impact bash $HERE/extra/acmg_impact.sh
  $PY $HERE/extra/acmg_compare.py $EXTRA_DIR/acmg_impact $(for i in $ACMG_SAMPLES; do printf '%s ' ${i%%=*}; done)
}
stage_reports() {
  # regenerate study and gene-level tables/figures from the kept archive (in place), then copy into the repository
  for r in exonic non_exonic; do $PY $HERE/freq_study_report.py $KEEP/study/$r $r; done
  $PY $HERE/gene_burden_report.py $KEEP/gene_level
  for r in exonic non_exonic; do
    mkdir -p $REPO/results/02_frequency_vs_gnomad/$r/tables $REPO/results/02_frequency_vs_gnomad/$r/figures
    cp $KEEP/study/$r/tables/*.tsv $KEEP/study/$r/counts.tsv $REPO/results/02_frequency_vs_gnomad/$r/tables/
    cp $KEEP/study/$r/figures/* $REPO/results/02_frequency_vs_gnomad/$r/figures/
  done
  mkdir -p $REPO/results/03_gene_level/tables $REPO/results/03_gene_level/figures $REPO/results/01_union $REPO/results/02_frequency_vs_gnomad/supplementary $REPO/results/00_provenance
  cp $KEEP/gene_level/tables/*.tsv $REPO/results/03_gene_level/tables/
  cp $KEEP/gene_level/figures/* $REPO/results/03_gene_level/figures/
  cp $KEEP/union/summary_per_dataset.tsv $KEEP/union/summary_how_many.tsv $REPO/results/01_union/
  cp $KEEP/study/extraction_stats.tsv $REPO/results/01_union/
  cp $KEEP/study/exonic/markers.tsv.gz $REPO/results/02_frequency_vs_gnomad/supplementary/exonic_markers_common_locally_rare_in_gnomad.tsv.gz
  cp $KEEP/study/exonic/clinvar_hits.tsv.gz $REPO/results/02_frequency_vs_gnomad/supplementary/exonic_clinvar_variants_af5pct_somewhere.tsv.gz
  cp $KEEP/study/non_exonic/clinvar_hits.tsv.gz $REPO/results/02_frequency_vs_gnomad/supplementary/non_exonic_clinvar_variants_af5pct_somewhere.tsv.gz
  for a in 04:founder_variants 05:carrier_frequency 06:unreliable_regions 07:ancestry_markers 08:acmg_impact 09:pgr_comparison; do
    n=${a%%:*}; a=${a#*:}
    [ -d $EXTRA_DIR/$a/tables ] || continue
    mkdir -p $REPO/results/${n}_$a/tables $REPO/results/${n}_$a/figures
    cp $EXTRA_DIR/$a/tables/* $REPO/results/${n}_$a/tables/
    cp $EXTRA_DIR/$a/figures/* $REPO/results/${n}_$a/figures/ 2>/dev/null
  done
  cp $KEEP/provenance/*.md5 $REPO/results/00_provenance/ 2>/dev/null
  ( cd $REPO && $PY scripts/analyses/pgr_comparison_figure.py && $PY scripts/make_results_md.py && bash scripts/build_report_md.sh )
}

[ $# -ge 1 ] || { sed -n '2,25p' $0; exit 1; }
stages=()
for s in "$@"; do
  case $s in
    all-raw) stages+=(download prepare union study clinvar gene keep extra acmg reports) ;;
    all-kept) stages+=(extra reports) ;;
    *) stages+=($s) ;;
  esac
done
for s in "${stages[@]}"; do
  declare -F stage_$s >/dev/null || { echo "unknown stage: $s"; exit 1; }
  log "stage $s started"
  stage_$s > $DATA/logs/$s.log 2>&1 && log "stage $s finished" || { log "stage $s FAILED (see $DATA/logs/$s.log)"; exit 1; }
done
