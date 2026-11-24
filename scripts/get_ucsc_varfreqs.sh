#!/bin/bash
# 1. Move the ANNOVAR population-frequency tables from /root (full) to /home, leaving symlinks so
#    ANNOVAR's humandb keeps working.
# 2. Download the per-cohort frequency VCFs of UCSC's hg38 varFreqs collection
#    (https://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs/). Skipped: cohorts we already have from
#    ANNOVAR (allofus, abraom, hrc), the combined varFreqsAll.bb and the ALFA per-population .bb
#    copies of ALFA.vcf.gz, and GenomeAsia's older indel file. http, because the NPU CA bundle
#    rejects the hgdownload certificate.
set -uo pipefail
source /etc/profile.d/npu-proxy.sh 2>/dev/null || true
H=/root/naser/annotation/annovar/humandb
A=/home/NASER/annotation/population_study/data/raw/annovar_popfreq
U=/home/NASER/annotation/population_study/data/raw/ucsc_varfreqs
log() { echo "[$(date +%F' '%T)] $*"; }

log "== moving ANNOVAR frequency tables to $A"
mkdir -p $A
for f in hg38_gnomad41_exome hg38_gnomad41_genome hg38_gnomad_exome hg38_gnomad_genome hg38_kaviar_20150923 \
         hg38_exac03 hg38_gme hg38_ALL.sites.2015_08 hg38_AFR.sites.2015_08 hg38_AMR.sites.2015_08 \
         hg38_EAS.sites.2015_08 hg38_EUR.sites.2015_08 hg38_SAS.sites.2015_08; do
  for x in $f.txt $f.txt.idx; do
    [ -f $H/$x ] && [ ! -L $H/$x ] || continue
    rsync -a $H/$x $A/ && [ "$(md5sum < $H/$x)" = "$(md5sum < $A/$x)" ] && rm $H/$x && ln -s $A/$x $H/$x && log "moved $x" || log "FAILED $x"
  done
done
df -h / | tail -1

log "== downloading UCSC varFreqs VCFs to $U"
B=http://hgdownload.soe.ucsc.edu/gbdb/hg38/varFreqs
FILES="alfa/ALFA.vcf.gz colorsDb/colorsDbSnv.vcf.gz finngen/finnge_R12_annotated_variants_v1.vcf.gz
ga4k/ga4kSnv.vcf.gz ga100k/ga100k.indels.vcf.gz ga100k/ga100k.subst.vcf.gz gonl/gonl.vcf.gz gregor/gregor.vcf.gz
hgdp1kFreq/hgdp1k.freq.vcf.gz hostseq/hostseq.vcf.gz indigenomes/IndiGenomes_Variants.vcf.gz kova/kova.v7.vcf.gz
mcps/mcps.freq.vcf.gz npm/SG10K_Health_r5.3.2.sites.vcf.bgz saudi/saudi.vcf.gz schema/SCHEMA_variant_results_withAF.vcf.gz
sfari/SPARK.iWES_v3.2024_08.deepvariant.norm.vcf.gz sfari/wgs_12519_genome.deepvariant.norm.vcf.gz sgdpFreq/sgdp.freq.vcf.gz
svatalog/svatalog.vcf.gz swefreq/swegen_frequencies_fixploidy_GRCh38_20190204.vcf.gz swefreq/terms_of_use.txt
tommo61kjpn/tommo-61kjpn-20250616-GRCh38-snvindel-af-autosome.vcf.gz topmed/topmed10.vcf.gz ukbb/ukbb.vcf.gz wbbc/wbbc.vcf.gz"
get() { local p=$1; mkdir -p $U/$(dirname $p)
  wget -c -q --tries=0 --waitretry=30 -O $U/$p $B/$p && log "downloaded $p" || log "FAILED $p"; }
export -f get log; export U B
for p in $FILES; do echo $p; case $p in *.txt) ;; *) echo $p.tbi;; esac; done | xargs -P 6 -I{} bash -c 'get {}'

log "== checks"
for f in $(find $U -name "*.vcf.gz" -o -name "*.vcf.bgz"); do bgzip -t $f 2>/dev/null && log "OK $f" || log "CORRUPT $f"; done
du -sh $U
log "UCSC_DONE"
