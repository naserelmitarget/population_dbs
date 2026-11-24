#!/bin/bash
# Convert UCSC's combined varFreqsAll.bb (downloaded by get_ucsc_missing.sh) to a bgzipped, tabix-indexed table with a header
# taken from the file's own autoSql. Tools: bioconda ucsc-bigbedtobed/ucsc-bigbedinfo (env "ucsc").
set -uo pipefail
source /etc/profile.d/npu-proxy.sh 2>/dev/null || true
export PATH=/root/miniforge3/envs/ucsc/bin:$PATH
D=${U:-/home/NASER/annotation/population_study/data/raw/ucsc_varfreqs}/varFreqsAll; mkdir -p $D/tmp; cd $D
log() { echo "[$(date +%F' '%T)] $*"; }
[ -s varFreqsAll.bb ] || { log "no varFreqsAll.bb (run get_ucsc_missing.sh)"; exit 1; }
md5sum varFreqsAll.bb > varFreqsAll.bb.md5
bigBedInfo varFreqsAll.bb > varFreqsAll.info.txt
bigBedInfo -as varFreqsAll.bb > varFreqsAll.as
hdr=$(sed -n '/^(/,/^)/p' varFreqsAll.as | awk 'NF > 1 {print $2}' | tr -d ';' | paste -sd'\t')
# one chromosome per job, then concatenated in bigBed order (already sorted)
chroms=$(bigBedInfo -chroms varFreqsAll.bb | awk '/^chromCount/ {f = 1; next} f && NF >= 3 {print $1}')
export -f log
echo "$chroms" | xargs -P 25 -I{} bash -c 'bigBedToBed -chrom={} varFreqsAll.bb tmp/{}.bed && log "converted {}"'
{ printf '#%s\n' "$hdr"; for c in $chroms; do cat tmp/$c.bed; done; } | bgzip -@ 32 > varFreqsAll.tsv.gz \
  && tabix -s1 -b2 -e3 -0 varFreqsAll.tsv.gz && rm -rf tmp
log "rows: $(zcat varFreqsAll.tsv.gz | grep -vc '^#')  (bigBedInfo itemCount: $(awk '/itemCount/ {print $2}' varFreqsAll.info.txt))"
log "BB_DONE"
