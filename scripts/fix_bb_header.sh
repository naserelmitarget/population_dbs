#!/bin/bash
# varFreqsAll.bb data columns are grouped per cohort (AC, AF, then that cohort population AC/AF), but the
# autoSql lists all main AC/AF first. Verified on 20M rows: single-source rows fill exactly these columns
# (spk 37-38, wgs 39-40, ga 41-56, kov 57-58, tom 59-60, mcp 61-68, fin 71-72, sau 73-74, swe 75-76, top 77-78).
set -eu
cd ${U:-/home/NASER/annotation/population_study/data/raw/ucsc_varfreqs}/varFreqsAll
h=$(zcat varFreqsAll.tsv.gz | head -1 | cut -f1-22)
for c in "aou:AFR AMR EAS EUR OCE SAS" spk wgs "ga:AFR AMR NEA OCE SAS SEA WER" kov tom "mcp:AFR EUR IMX" ind fin sau swe top abr; do
  k=${c%%:*}; h+=$(printf "\t%sAC\t%sAF" $k $k)
  [ "$k" != "$c" ] && for p in ${c#*:}; do h+=$(printf "\t%sAC_%s\t%sAF_%s" $k $p $k $p); done
done
[ $(echo "$h" | awk -F"\t" "{print NF}") -eq 80 ] || { echo "header has wrong column count"; exit 1; }
echo "$h" | tr "\t" "\n" | awk "{print NR\"\t\"\$0}" > varFreqsAll.columns.tsv
{ echo "$h"; zcat varFreqsAll.tsv.gz | tail -n +2; } | bgzip -@ 32 > varFreqsAll.fixed.tsv.gz
tabix -s1 -b2 -e3 -0 varFreqsAll.fixed.tsv.gz
mv varFreqsAll.fixed.tsv.gz varFreqsAll.tsv.gz && mv varFreqsAll.fixed.tsv.gz.tbi varFreqsAll.tsv.gz.tbi
echo "HEADER_FIXED $(date +%F" "%T)"
