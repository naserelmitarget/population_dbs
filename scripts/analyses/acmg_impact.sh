#!/bin/bash
# ACMG impact of population datasets beyond gnomAD, on real genomes (exonic regions).
# BIAS-2015 (unchanged) classifies each sample's variants three times:
#   A   gnomAD only (BIAS as shipped)
#   B1  frequency raised to the highest AF in any population dataset
#   B2  frequency raised to the second-highest AF among independent datasets (at least two datasets agree)
# Population datasets: the frequency-study datasets except disease cohorts (SPARK, SFARI WGS, GREGoR, GA4K);
# gnomAD v2 and ExAC count as gnomAD (not independent).
#
# Samples: ACMG_SAMPLES="NAME=path.vcf.gz NAME2=path2.vcf.gz ..." (GRCh38 VCFs; normalised here).
# Output: $W/<NAME>/ (VEP JSON, three BIAS tables), $W/population_af_lookup.tsv (all samples, one pass).
set -uo pipefail
ACMG_SAMPLES=${ACMG_SAMPLES:-"HUV00001=/home/NASER_FILES/Naser_VEP/results/HUV00001_refseq_32cpu/HUV00001.norm.vcf.gz"}
K=${K:-/home/NASER/annotation/population_study/data/archive}
W=${W:-/home/NASER/annotation/population_study/data/work/extra_analyses/acmg_impact}
V=${V:-/home/NASER_FILES/Naser_VEP/vep_cache}                       # VEP 116 RefSeq cache (read-only)
FA=${FA:-/home/NASER_FILES/annovar/hg38/reference/hg38.fa}
ACMG=${ACMG:-/home/NASER/annotation/acmg}
BIAS=${BIAS:-/root/naser/annotation/acmg/BIAS-2015}                 # read-only use
PATHS=$ACMG/bias_hg38_current/hg38_vep_required_paths.json
CLINVAR=${BIAS_CLINVAR:-$ACMG/downloads/clinvar_20261011.vcf.gz}   # the ClinVar the BIAS data was built from
VEP_IMAGE=${VEP_IMAGE:-huvya/vep:116}
BIAS_PY=${BIAS_PY:-python3}                                         # BIAS-2015 needs only the standard library
HERE=$(cd $(dirname $0) && pwd)
UNZIP="gzip -dc"; command -v pigz >/dev/null && UNZIP="pigz -dc"
export LC_ALL=C
log() { echo "[$(date +%F' '%T)] $*"; }
mkdir -p $W; cd $W
awk -F'[:-]' -v OFS='\t' '{ print $1, $2 - 1, $3 }' $K/study/regions/gnomad_v4_exome_calling.ukb_broad_pad50.intervals \
  | sort -k1,1 -k2,2n > exome_intervals.bed

# 1. per sample: split multi-allelic records, left-align, keep PASS variants inside the exome calling intervals, VEP
names=""
for item in $ACMG_SAMPLES; do
  s=${item%%=*}; vcf=${item#*=}; names="$names $s"
  mkdir -p $W/$s
  ( cd $W/$s
    bcftools norm -m -any -f $FA $vcf -Oz -o $s.norm.vcf.gz 2> norm.log && tabix -f -p vcf $s.norm.vcf.gz
    bcftools view -f PASS,. -R ../exome_intervals.bed $s.norm.vcf.gz 2>/dev/null | bcftools norm -d exact -Oz -o $s.exonic.vcf.gz \
      && tabix -f -p vcf $s.exonic.vcf.gz && rm -f $s.norm.vcf.gz $s.norm.vcf.gz.tbi
    log "$s exonic PASS variants: $(bcftools view -H $s.exonic.vcf.gz | wc -l)"
    docker run --rm --cpus 48 -u "$(id -u):$(id -g)" -e HOME=/tmp -v /home/NASER_FILES:/home/NASER_FILES:ro -v /home/NASER:/home/NASER \
      -w $W/$s $VEP_IMAGE vep -i $s.exonic.vcf.gz -o $s.vep.json --json --force_overwrite --no_stats --offline --cache \
      --dir_cache $V --refseq --cache_version 116 --assembly GRCh38 --fasta $FA --everything --fork 24 --buffer_size 5000 \
      --custom file=$CLINVAR,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNVC%CLNSIGSCV%CLNSIGCONF \
      --dir_plugins /opt/conda/share/ensembl-vep-116.3-0 \
      --plugin AlphaMissense,file=$ACMG/vep_plugins_hg38/AlphaMissense_hg38.tsv.gz \
      --plugin REVEL,file=$ACMG/vep_plugins_hg38/revel_grch38.tsv.gz > vep.log 2>&1
    log "$s VEP records: $(wc -l < $s.vep.json)" ) &
done
wait

# 2. population AF lookup for all samples' variants in one pass over the kept exonic variant table
for s in $names; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $W/$s/$s.exonic.vcf.gz
done | awk -F'\t' -v OFS='\t' '
  { c = $1; sub(/^chr/, "", c); if (c == "MT") c = "M"; s = $2 + 0; r = $3; a = $4
    if (r !~ /^[ACGTN]+$/ || a !~ /^[ACGTN]+$/) next
    k = 0; while (k < length(r) && k < length(a) && substr(r, k+1, 1) == substr(a, k+1, 1)) k++
    r = substr(r, k+1); a = substr(a, k+1); s = s + k
    while (length(r) > 0 && length(a) > 0 && substr(r, length(r), 1) == substr(a, length(a), 1)) { r = substr(r, 1, length(r)-1); a = substr(a, 1, length(a)-1) }
    if (r == "" && a == "") next
    if (r == "") { s = s - 1; e = s } else e = s + length(r) - 1
    print c, s, e, (r == "" ? "-" : r), (a == "" ? "-" : a) }' | sort -u > sample.keys
$UNZIP $K/gene_level/exonic_variants_annotated.tsv.gz | awk -F'\t' -v OFS='\t' '
  BEGIN { split("spark_wes sfari_wgs gregor ga4k gnomad_exome_v2 gnomad_genome_v2 exac03", x, " "); for (i in x) skip[x[i]] = 1
          print "chr", "start", "end", "ref", "alt", "gnomad_af", "af_b1_max", "af_b2_second", "n_datasets", "datasets" }
  NR == FNR { s[$1 OFS $2 OFS $3 OFS $4 OFS $5] = 1; next }
  FNR == 1 { next }
  ($1 OFS $2 OFS $3 OFS $4 OFS $5) in s {
    n = split($7, d, ";"); m1 = -1; m2 = -1; c = 0; keep = ""
    for (i = 1; i <= n; i++) { split(d[i], kv, "="); if (kv[1] == "" || (kv[1] in skip) || kv[2] == "." || kv[2] == "") continue
      f = kv[2] + 0; c++; keep = keep (keep == "" ? "" : ";") d[i]
      if (f > m1) { m2 = m1; m1 = f } else if (f > m2) m2 = f }
    print $1, $2, $3, $4, $5, $6, (m1 < 0 ? "." : m1), (m2 < 0 ? "." : m2), c, keep }' sample.keys - > population_af_lookup.tsv
log "variants (all samples) with population AF: $(($(wc -l < population_af_lookup.tsv) - 1)) of $(wc -l < sample.keys)"

# 3. BIAS-2015: A, B1, B2 for every sample, in parallel
run_bias() {   # run_bias <sample> <name> <lookup|none> <column>
  local s=$1 name=$2 lk=$3 col=$4
  cd $W/$s
  $BIAS_PY $HERE/bias_with_population_af.py $BIAS $W/$s/$s.vep.json $PATHS $W/$s/bias_$name.raw.tsv $lk $col > bias_$name.log 2>&1
  $BIAS_PY $HERE/bias_add_codes.py bias_$name.raw.tsv bias_$name.tsv > /dev/null 2>&1 && rm -f bias_$name.raw.tsv
  log "$s BIAS $name: $(($(wc -l < bias_$name.tsv) - 1)) rows; $(grep 'population AF' bias_$name.log)"
}
for s in $names; do
  run_bias $s A none "" &
  run_bias $s B1 $W/population_af_lookup.tsv af_b1_max &
  run_bias $s B2 $W/population_af_lookup.tsv af_b2_second &
done
wait
log "ACMG_DONE"
