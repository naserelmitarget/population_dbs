#!/bin/bash
# UCSC varFreqs cohorts that are NOT inside varFreqsAll.bb (and not ANNOVAR datasets), plus the .bb itself.
# The .bb already holds SPARK, SFARI-WGS, GenomeAsia, KOVA, ToMMo, MCPS, IndiGen, FinnGen, Saudi, SweGen, TOPMed
# (and All of Us / ABraOM, which we take from ANNOVAR). http: the NPU CA bundle rejects hgdownload's TLS cert.
set -uo pipefail
source /etc/profile.d/npu-proxy.sh 2>/dev/null || true
export PATH=/root/miniforge3/envs/ucsc/bin:$PATH
U=${U:-/home/NASER/annotation/population_study/data/raw/ucsc_varfreqs}
B=http://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs
log() { echo "[$(date +%F' '%T)] $*"; }
FILES="varFreqsAll.bb alfa/ALFA.vcf.gz colorsDb/colorsDbSnv.vcf.gz ga4k/ga4kSnv.vcf.gz gonl/gonl.vcf.gz gregor/gregor.vcf.gz
hgdp1kFreq/hgdp1k.freq.vcf.gz hostseq/hostseq.vcf.gz npm/SG10K_Health_r5.3.2.sites.vcf.bgz schema/SCHEMA_variant_results_withAF.vcf.gz
sgdpFreq/sgdp.freq.vcf.gz svatalog/svatalog.vcf.gz ukbb/ukbb.vcf.gz wbbc/wbbc.vcf.gz"
get() { local p=$1 d=$U/$(dirname $1); [ $p = varFreqsAll.bb ] && d=$U/varFreqsAll; mkdir -p $d
  for try in 1 2 3 4 5; do
    aria2c -x8 -s8 -c --file-allocation=none --max-tries=0 --retry-wait=30 --summary-interval=0 --console-log-level=warn \
      -d $d -o $(basename $p) $B/$p >/dev/null 2>&1 && { log "downloaded $p"; return 0; }; sleep 60; done
  log "FAILED $p"; }
export -f get log; export U B
for p in $FILES; do echo $p; [ $p = varFreqsAll.bb ] || echo $p.tbi; done | xargs -P 4 -I{} bash -c 'get {}'
for f in $(find $U -name "*.vcf.gz" -o -name "*.vcf.bgz"); do bgzip -t $f 2>/dev/null && log "OK $f" || log "CORRUPT $f"; done
du -sh $U
log "UCSC_DONE"
