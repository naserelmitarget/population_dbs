#!/bin/bash
# Genes and variants whose population frequencies should not be trusted, checked against GIAB difficult regions.
# No raw frequency files needed. Output: $W/unreliable_regions/{tables,figures,giab}
set -uo pipefail
K=${K:-/home/NASER/annotation/population_study/data/archive}
W=${W:-/home/NASER/annotation/population_study/data/work/extra_analyses}
O=$W/unreliable_regions
G=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all
export LC_ALL=C K
log() { echo "[$(date +%F' '%T)] $*"; }
mkdir -p $O/giab $O/tables $O/figures $O/tmp
source /etc/profile.d/npu-proxy.sh 2>/dev/null || true

# 1. GIAB v3.5 GRCh38 stratifications: segmental duplications, low mappability (all)
for f in SegmentalDuplications/GRCh38_segdups.bed.gz Mappability/GRCh38_lowmappabilityall.bed.gz; do
  [ -s $O/giab/$(basename $f) ] || curl -sS -o $O/giab/$(basename $f) $G/$f
done
for b in segdups lowmappabilityall; do
  zcat $O/giab/GRCh38_$b.bed.gz | grep -v '^#' | sed 's/^chr//' | sort -k1,1 -k2,2n | bedtools merge > $O/tmp/$b.bed
done
log "GIAB: segdups $(awk '{s+=$3-$2} END {printf "%.1f Mb", s/1e6}' $O/tmp/segdups.bed), low mappability $(awk '{s+=$3-$2} END {printf "%.1f Mb", s/1e6}' $O/tmp/lowmappabilityall.bed)"

# 2. per gene: exonic footprint (distinct positions of observed variants) and the share inside each GIAB set
$(command -v pigz >/dev/null && echo "pigz -dc" || echo "gzip -dc") $K/gene_level/exonic_variants_annotated.tsv.gz | awk -F'\t' -v OFS='\t' 'NR > 1 && $11 != "" && $11 != "-" && $12 == "protein_coding" {
    print $1, $2 - 1, $2, $11 "|" $10 }' | sort -u -S 32G --parallel=16 -k1,1 -k2,2n > $O/tmp/gene_positions.bed
for b in segdups lowmappabilityall; do
  bedtools intersect -sorted -c -a $O/tmp/gene_positions.bed -b $O/tmp/$b.bed \
    | awk -F'\t' -v OFS='\t' '{ n[$4]++; if ($5 > 0) h[$4]++ } END { for (g in n) print g, n[g], h[g] + 0 }' | sort > $O/tmp/gene_$b.tsv
done
awk -F'\t' -v OFS='\t' '{ c[$4] = $1; if (!($4 in s) || $2 < s[$4]) s[$4] = $2; if ($3 > e[$4]) e[$4] = $3 }
  END { for (g in c) print c[g], s[g], e[g], g }' $O/tmp/gene_positions.bed | sort -k1,1 -k2,2n > $O/tmp/gene_span.bed
log "genes with footprint: $(wc -l < $O/tmp/gene_span.bed)"

# 3. evidence, flags, validation, figures
${PY:-/root/miniforge3/envs/ahm-env/bin/python} $(dirname $0)/unreliable_regions.py $W
rm -rf $O/tmp/gene_positions.bed
log "UNRELIABLE_DONE"
