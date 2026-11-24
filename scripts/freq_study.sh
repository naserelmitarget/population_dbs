#!/bin/bash
# Frequency study: every population dataset vs gnomAD v4.1, split into two regions with separate outputs:
#   exonic/      inside gnomAD v4 exome calling intervals (UKB + Broad targets, +-50 bp)
#                gnomAD reference = exome AF; genome AF where the exome table has no AF for the variant
#   non_exonic/  everything else; gnomAD reference = genome AF only (exome-table calls here are ignored)
#                exome-only datasets (wes flag) are left out of this region
# Read-only on all inputs. Variants normalised exactly as in freq_union.sh. (v1, one combined reference:
# freq_study_v1_combined.sh)
#
# Per dataset allele: AF (+ AC / AN / FILTER where the source has them).
#   dropped: AC == 0, FILTER not PASS/"."; AF set to missing (class noaf) when AN < MINAN
#   (AN from INFO, Kaviar_AN, or AC/AF): the UCSC SPARK table has records like AC=3 AF=1.
# Unreliable sources left out: GenomeAsia (UCSC .bb frequencies: ~25% of variants at AF 10-90% absent from both
#   TOPMed and gnomAD), ALFA indels (UCSC ALFA indel positions are off by one and still half unmatched after
#   correction; ALFA SNVs are used), SGDP (AN 2 per site), SCHEMA (no single AF), IndiGen (no AF).
# gnomAD overall AF 0 (carriers only among QC-failed genotypes) is treated as no AF: class qcfail, and in exonic
#   regions the genome AF is used instead when it has one.
# gnomAD views:  overall = AF;  popmax = max(AF, afr amr asj eas fin mid nfe sas);  matched = AF of the dataset's
#   own gnomAD group (country datasets). gnomAD classes: absent | qcfail (row, no AF) | rare <5% | low 5-10% | poly >10%
# Dataset classes: noaf | r1 <1% | r5 1-5% | low 5-10% | poly >10%.  Types: snv | indel | mnv
# Outputs per region: counts.tsv, markers.tsv.gz, clinvar_candidates.tsv.gz, pca_snvs.tsv.gz
set -uo pipefail
A=${A:-/home/NASER/annotation/population_study/data/raw/annovar_popfreq}
U=${U:-/home/NASER/annotation/population_study/data/raw/ucsc_varfreqs}
BB=${BB:-$U/varFreqsAll/varFreqsAll.tsv.gz}
W=${W:-/home/NASER/annotation/population_study/data/work/freq_study}
IV=${IV:-$W/regions/gnomad_v4_exome_calling.ukb_broad_pad50.intervals}
MINAN=${MINAN:-100}
DUMP=${DUMP:-0}            # 1: also write every exonic variant with its AF per dataset (exonic/vars.tsv.gz, gene step)
JOBS=${JOBS:-120}
export LC_ALL=C
log() { echo "[$(date +%F' '%T)] $*"; }
mkdir -p $W/keys $W/tmp $W/regions; cd $W
for r in exonic non_exonic; do mkdir -p $W/$r/agg $W/$r/markers $W/$r/cand $W/$r/pca; done
mkdir -p $W/exonic/vars

# name : file : format : AF field : gnomAD group (- none) : wes (exome-only dataset) or -
#   annovar (AF field = header name) | 1000g | vcf (AF field = INFO key) | vcf_snv (SNVs only) | bb (AF field = cohort code)
DS=(gnomad_exome_v2:$A/hg38_gnomad_exome.txt:annovar:gnomAD_exome_ALL:-:wes
    gnomad_genome_v2:$A/hg38_gnomad_genome.txt:annovar:gnomAD_genome_ALL:-:-
    kaviar:$A/hg38_kaviar_20150923.txt:annovar:Kaviar_AF:-:-
    exac03:$A/hg38_exac03.txt:annovar:ExAC_ALL:-:wes
    1000g2015aug:$A/hg38_ALL.sites.2015_08.txt:1000g:-:-:-
    gme:$A/hg38_gme.txt:annovar:GME_AF:mid:wes
    allofus:$A/hg38_allofus.txt:annovar:gvs_all_af:-:-
    spark_wes:$BB:bb:spk:-:wes
    sfari_wgs:$BB:bb:wgs:-:-
    kova:$BB:bb:kov:eas:-
    tommo:$BB:bb:tom:eas:-
    mcps:$BB:bb:mcp:amr:-
    finngen:$BB:bb:fin:fin:-
    saudi:$BB:bb:sau:mid:-
    swegen:$BB:bb:swe:nfe:-
    topmed:$BB:bb:top:-:-
    alfa:$U/alfa/ALFA.vcf.gz:vcf_snv:AF_GLB:-:-
    colorsdb:$U/colorsDb/colorsDbSnv.vcf.gz:vcf:AF:-:-
    ga4k:$U/ga4k/ga4kSnv.vcf.gz:vcf:AF:-:-
    gonl:$U/gonl/gonl.vcf.gz:vcf:AF:nfe:-
    gregor:$U/gregor/gregor.vcf.gz:vcf:AF:-:-
    hgdp1kg:$U/hgdp1kFreq/hgdp1k.freq.vcf.gz:vcf:AF:-:-
    hostseq:$U/hostseq/hostseq.vcf.gz:vcf:AF:-:-
    sg10k_npm:$U/npm/SG10K_Health_r5.3.2.sites.vcf.bgz:vcf:AF:eas:-
    svatalog:$U/svatalog/svatalog.vcf.gz:vcf:AF:-:-
    ukbb:$U/ukbb/ukbb.vcf.gz:vcf:AF:nfe:-
    wbbc:$U/wbbc/wbbc.vcf.gz:vcf:AF:eas:-)
PCA_DS="tommo kova wbbc sg10k_npm mcps finngen swegen gonl ukbb saudi gme"
GPOPS="afr amr asj eas fin mid nfe sas"
CHROMS=${CHROMS:-"$(seq 1 22 | tr "\n" " ")X Y M"}
NAMES=$(for d in "${DS[@]}"; do printf '%s,' "$(cut -d: -f1 <<< "$d")"; done); NAMES=${NAMES%,}
DGRP=$(for d in "${DS[@]}"; do printf '%s,' "$(cut -d: -f5 <<< "$d")"; done); DGRP=${DGRP%,}
DWES=$(for d in "${DS[@]}"; do printf '%s,' "$(cut -d: -f6 <<< "$d")"; done); DWES=${DWES%,}
BBMAP=$(for i in "${!DS[@]}"; do IFS=: read -r _ _ fmt fld _ <<< "${DS[$i]}"; [ $fmt = bb ] && printf '%s=%s,' $fld $((i + 2)); done); BBMAP=${BBMAP%,}
for d in "${DS[@]}"; do IFS=: read -r _ f _ <<< "$d"; [ -s $f ] || { log "missing $f"; exit 1; }; done
[ -s $IV ] || { log "missing $IV"; exit 1; }

# exome calling intervals: per chromosome, sorted and merged, 1-based inclusive
awk -F'[:-]' -v OFS='\t' '{ c = $1; sub(/^chr/, "", c); if (c == "MT") c = "M"; print c, $2, $3 }' $IV \
  | sort -k1,1 -k2,2n -k3,3n \
  | awk -F'\t' -v d=$W/regions '{ if ($1 != c || $2 > e + 1) { if (c != "") print c, s, e > (d "/iv." c ".txt"); c = $1; s = $2; e = $3 }
                                   else if ($3 > e) e = $3 }
                                 END { print c, s, e > (d "/iv." c ".txt") }'
for c in $CHROMS; do touch regions/iv.$c.txt; done

# shared normaliser: emit(file prefix, chrom, pos, ref, alt, payload) -> "<s> <e> <ref> <alt> <payload>"
NORM='
function emit(pre, c, s, r, a, pay,    ins, k, e) {
  r = toupper(r); a = toupper(a)
  if (r !~ /^([ACGTN]+|-)$/ || a !~ /^([ACGTN]+|-)$/ || r == a) { bad++; return }
  ins = (r == "-"); if (r == "-") r = ""; if (a == "-") a = ""
  k = 0; while (k < length(r) && k < length(a) && substr(r, k+1, 1) == substr(a, k+1, 1)) k++
  r = substr(r, k+1); a = substr(a, k+1); s = s + k
  while (length(r) > 0 && length(a) > 0 && substr(r, length(r), 1) == substr(a, length(a), 1)) {
    r = substr(r, 1, length(r)-1); a = substr(a, 1, length(a)-1) }
  if (r == "" && a == "") { bad++; return }
  if (r == "") { if (!ins) s = s - 1; e = s } else e = s + length(r) - 1
  if (r == "") r = "-"; if (a == "") a = "-"
  print s, e, r, a, pay > (pre "." c ".txt"); n_ok++
}
function chrom(x) { sub(/^chr/, "", x); if (x == "MT") x = "M"; return (x ~ /^([0-9]+|X|Y|M)$/) ? x : "" }
function num(x) { return (x == "" || x == "." || x == "NA" || x == "nan") ? "." : x + 0 }'
export NORM W MINAN GPOPS BB BBMAP

# one chromosome slice of an ANNOVAR table, located through its .idx (1 kb bins, byte offsets)
chunk() {   # chunk <file> <chrom>
  local s e
  # bins are not listed in file order, so take the smallest start and largest end of the chromosome's bins
  read s e < <(awk -F'\t' -v c=$2 '($1 == c || $1 == "chr" c) { if (s == "" || $3 + 0 < s) s = $3 + 0; if ($4 + 0 > e) e = $4 + 0 } END { print s, e }' $1.idx)
  [ -n "$s" ] && tail -c +$((s + 1)) $1 | head -c $((e - s))
}
# header column numbers by name (gnomad41_* prefixes stripped): hcols <file> <name...>
hcols() { local f=$1; shift; head -1 $f | awk -F'\t' -v want="$*" '{ for (j = 1; j <= NF; j++) { h = $j; sub(/^gnomad41_(exome|genome)_/, "", h); col[h] = j }
                                                                     n = split(want, w, " "); for (i = 1; i <= n; i++) printf "%s ", (w[i] in col) ? col[w[i]] : 0 }'; }
export -f chunk hcols log

# 1. extraction jobs, one per (source, chromosome)
job() {   # job <kind> <index> <file> <chrom or contig> [format] [AF field]
  local kind=$1 i=$2 f=$3 c=$4 fmt=${5:-} fld=${6:-} cols
  case $kind in
  g) cols=$(hcols $f AF $(for p in $GPOPS; do printf 'AF_%s ' $p; done))
     chunk $f $c | awk -F'\t' -v OFS='\t' -v pre="$W/keys/g$i" -v tag=$i -v cols="$cols" -v st="$W/keys/g$i.$c.stats" "$NORM"'
       BEGIN { n = split(cols, C, " ") }
       /^#/ { next }
       { cc = chrom($1); if (cc == "") next
         # overall AF 0 = seen only in genotypes failing gnomAD QC: treated as no AF (exonic falls back to genomes)
         pay = tag; for (j = 1; j <= n; j++) { x = num(C[j] ? $(C[j]) : "."); if (j == 1 && x == 0) x = "."; pay = pay OFS x }
         emit(pre, cc, $2 + 0, $4, $5, pay) }
       END { printf "g%d\t%d\t%d\t0\t0\t0\n", tag, n_ok, bad > st }' ;;
  a) local afc=0 anc=0; [ $fmt = annovar ] && read afc anc < <(hcols $f $fld Kaviar_AN)
     chunk $f $c | awk -F'\t' -v OFS='\t' -v pre="$W/keys/d$i" -v tag=$((i + 2)) -v fmt=$fmt -v afc=$afc -v anc=$anc \
                       -v st="$W/keys/d$i.$c.stats" -v MINAN=$MINAN "$NORM"'
       /^#/ { next }
       { cc = chrom($1); if (cc == "") next
         if (fmt == "1000g") {
           if (!match($4, /^[0-9]+/)) { emit(pre, cc, $2 + 0, $3, $4, tag OFS num($5)); next }
           n = substr($4, 1, RLENGTH) + 0; a = substr($4, RLENGTH + 1)
           if (n == 0) r = "-"; else if (length($3) == n) r = $3; else { bad++; next }
           emit(pre, cc, $2 + 0, r, a == "" ? "-" : a, tag OFS num($5)); next }
         af = num($afc); if (af != "." && anc && $anc != "." && $anc + 0 < MINAN) { af = "."; lowan++ }
         emit(pre, cc, $2 + 0, $4, $5, tag OFS af) }
       END { printf "d%d\t%d\t%d\t0\t0\t%d\n", tag - 2, n_ok, bad, lowan > st }' ;;
  b) tabix -h $f chr$c | awk -F'\t' -v OFS='\t' -v W="$W" -v map="$BBMAP" -v c=$c -v MINAN=$MINAN "$NORM"'
       BEGIN { n = split(map, m, ","); for (j = 1; j <= n; j++) { split(m[j], kv, "="); tagof[kv[1]] = kv[2] } }
       /^#/ { sub(/^#/, ""); for (j = 1; j <= NF; j++) col[$j] = j; next }
       { cc = chrom($1); if (cc == "") next
         ns = split($col["sources"], src, ",")
         for (j = 1; j <= ns; j++) { code = src[j]; if (!(code in tagof)) continue
           ac = $col[code "AC"]; if (ac != "" && ac + 0 == 0) { ac0[code]++; continue }
           af = num($col[code "AF"]); if (ac != "" && af != "." && af > 0 && ac / af < MINAN) { af = "."; lowan[code]++ }
           k0 = n_ok; emit(W "/keys/d" (tagof[code] - 2), cc, $col["chromStart"] + 1, $col["ref"], $col["alt"], tagof[code] OFS af)
           if (n_ok > k0) kept[code]++ } }
       END { for (code in tagof) printf "d%d\t%d\t0\t%d\t0\t%d\n", tagof[code] - 2, kept[code], ac0[code], lowan[code] > (W "/keys/bb." code "." c ".stats") }' ;;
  v) tabix $f $c | awk -F'\t' -v OFS='\t' -v pre="$W/keys/d$i" -v tag=$((i + 2)) -v fmt=$fmt -v fld=$fld \
                     -v st="$W/keys/d$i.$c.stats" -v MINAN=$MINAN "$NORM"'
       { cc = chrom($1); if (cc == "") next
         if ($7 != "PASS" && $7 != ".") { filt++; next }
         delete info; m = split($8, kv, ";")
         for (j = 1; j <= m; j++) { p = index(kv[j], "="); if (p) info[substr(kv[j], 1, p-1)] = substr(kv[j], p+1) }
         na = split($5, alts, ","); split(info[fld], afs, ","); split(info["AC"], acs, ",")
         for (j = 1; j <= na; j++) {
           if (fmt == "vcf_snv" && (length($4) != 1 || length(alts[j]) != 1)) { nonsnv++; continue }
           if (("AC" in info) && acs[j] != "" && acs[j] + 0 == 0) { ac0++; continue }
           af = num((fld in info) ? afs[j] : ".")
           an = ("AN" in info) ? info["AN"] + 0 : ((acs[j] != "" && af != "." && af > 0) ? acs[j] / af : -1)
           if (af != "." && an >= 0 && an < MINAN) { af = "."; lowan++ }
           emit(pre, cc, $2 + 0, $4, alts[j], tag OFS af) } }
       END { printf "d%d\t%d\t%d\t%d\t%d\t%d\n", tag - 2, n_ok, bad + nonsnv, ac0, filt, lowan > st }' ;;
  esac
}
export -f job
{
  for c in $CHROMS; do
    echo "g 0 $A/hg38_gnomad41_exome.txt $c"; echo "g 1 $A/hg38_gnomad41_genome.txt $c"; echo "b - $BB $c"
    for i in "${!DS[@]}"; do IFS=: read -r _ f fmt fld _ <<< "${DS[$i]}"
      case $fmt in annovar|1000g) echo "a $i $f $c $fmt $fld";; esac; done
  done
  for i in "${!DS[@]}"; do IFS=: read -r _ f fmt fld _ <<< "${DS[$i]}"
    case $fmt in vcf*) for ct in $(tabix -l $f 2>/dev/null); do x=${ct#chr}; [ "$x" = MT ] && x=M
                         [[ " $CHROMS " == *" $x "* ]] && echo "v $i $f $ct $fmt $fld"; done;; esac
  done
} > tmp/jobs.txt
log "extracting: $(wc -l < tmp/jobs.txt) jobs, $JOBS at a time"
xargs -P $JOBS -L 1 bash -c 'job "$@"' _ < tmp/jobs.txt 2> >(grep -v "index file is older" >&2)
# stats per source: rows_kept rows_dropped AC0 filtered AF_dropped_low_AN
cat keys/*.stats | awk -F'\t' -v OFS='\t' -v names=$NAMES '
  BEGIN { split(names, nm, ","); print "source", "rows_kept", "rows_dropped", "AC0", "filtered", "AF_dropped_low_AN" }
  { for (j = 2; j <= 6; j++) s[$1, j] += $j; seen[$1] = 1 }
  END { for (x in seen) { n = (x ~ /^g/) ? "gnomad41_" (x == "g0" ? "exome" : "genome") : nm[substr(x, 2) + 1]
                          printf "%s", n; for (j = 2; j <= 6; j++) printf "\t%d", s[x, j]; print "" } }' > extraction_stats.tsv
column -t extraction_stats.tsv

# 2. per chromosome: merge sorted, assign region, classify each dataset allele against the region's gnomAD reference
one() {   # one <chrom>
  local c=$1 f
  for f in keys/[gd]*.$c.txt; do [ -f $f ] && sort -S 16G --parallel=4 -T tmp -t $'\t' -k1,1n -k2,2n -k3,3 -k4,4 -k5,5n $f -o $f.s && rm -f $f; done
  ls keys/[gd]*.$c.txt.s >/dev/null 2>&1 || return 0
  sort -m -T tmp -t $'\t' -k1,1n -k2,2n -k3,3 -k4,4 -k5,5n keys/[gd]*.$c.txt.s \
  | awk -F'\t' -v OFS='\t' -v c=$c -v names=$NAMES -v groups=$DGRP -v wes=$DWES -v pops="$GPOPS" -v pcads="$PCA_DS" -v W=$W \
        -v ivf=$W/regions/iv.$c.txt -v dump=$DUMP '
    BEGIN { nd = split(names, NM, ","); split(groups, GR, ","); split(wes, WES, ","); np = split(pops, P, " ")
            for (j = 1; j <= np; j++) PI[P[j]] = j
            npc = split(pcads, PC, " "); for (j = 1; j <= nd; j++) DI[NM[j]] = j
            ni = 0; while ((getline line < ivf) > 0) { split(line, q, " "); IS[++ni] = q[2] + 0; IE[ni] = q[3] + 0 }; ptr = 1 }
    function gcls(x) { return x == "." ? "qcfail" : (x < 0.05 ? "rare" : (x <= 0.10 ? "low" : "poly")) }
    function dcls(x) { return x == "." ? "noaf" : (x < 0.01 ? "r1" : (x < 0.05 ? "r5" : (x <= 0.10 ? "low" : "poly"))) }
    function mx(a, b) { return a == "." ? b : (b == "." ? a : (a > b ? a : b)) }
    function flush(   t, d, R, ex, src, ing, ov, pm, mg, j, g, cls, maxaf, who, row, any, x, pos) {
      if (k == "") return
      split(k, kk, "\t"); pos = kk[1] + 0
      t = (kk[3] == "-" || kk[4] == "-") ? "indel" : (length(kk[3]) == 1 && length(kk[4]) == 1 ? "snv" : "mnv")
      while (ptr <= ni && IE[ptr] < pos) ptr++
      ex = (ptr <= ni && IS[ptr] <= pos); R = ex ? "exonic" : "non_exonic"
      # gnomAD reference for this region
      if (ex) { ing = hasg[0] || hasg[1]; src = (hasg[0] && gaf[0] != ".") ? 0 : ((hasg[1] && gaf[1] != ".") ? 1 : (hasg[0] ? 0 : 1)) }
      else    { ing = hasg[1]; src = 1 }
      # per-group AF: from the reference table, or from genomes where that table did not assess the group
      # (gnomAD exomes leave 12-26% of group AFs empty at common sites captured by only some exome cohorts)
      ov = "."; pm = "."
      for (j = 1; j <= np; j++) pv[j] = (gp[src, j] != ".") ? gp[src, j] : gp[1, j]
      if (ing) { ov = gaf[src]; pm = ov; for (j = 1; j <= np; j++) pm = mx(pm, pv[j]) }
      maxaf = "."; who = ""
      if (dump && ex) {                                     # gene step: variant, gnomAD AF, AF in each dataset
        row = ""; for (d in daf) row = row (row == "" ? "" : ";") NM[d] "=" daf[d]
        print c, k, (ing ? ov : "absent"), row > (W "/exonic/vars/" c ".tsv") }
      for (d in daf) {
        if (!ex && WES[d] == "wes") continue
        cls = dcls(daf[d]); g = GR[d]
        cnt[R OFS NM[d] OFS t OFS cls OFS "overall" OFS (ing ? gcls(ov) : "absent")]++
        cnt[R OFS NM[d] OFS t OFS cls OFS "popmax"  OFS (ing ? gcls(pm) : "absent")]++
        mg = "."
        if (g != "-") { mg = ing ? pv[PI[g]] : "."; cnt[R OFS NM[d] OFS t OFS cls OFS "matched" OFS (ing ? gcls(mg) : "absent")]++ }
        if ((cls == "poly" || cls == "low") && (!ing || pm == "." || pm < 0.05))
          print c, k, NM[d], daf[d], ov, pm, mg > (W "/" R "/markers/" c ".tsv")
        if (daf[d] != "." && daf[d] >= 0.05) { if (maxaf == "." || daf[d] > maxaf) maxaf = daf[d]; who = who (who == "" ? "" : ";") NM[d] ":" daf[d] }
      }
      if (maxaf != ".") print c, k, maxaf, who, ov, pm > (W "/" R "/cand/" c ".tsv")
      if (t == "snv" && hasg[1] && gaf[1] != ".") {        # PCA matrix: SNVs with gnomAD genome AF, AF >= 1% in a country set
        row = ""; any = 0
        for (j = 1; j <= npc; j++) { d = DI[PC[j]]; x = ((d in daf) && (ex || WES[d] != "wes")) ? daf[d] : "NA"
                                     row = row OFS x; if (x != "NA" && x != "." && x >= 0.01) any = 1 }
        if (any) { for (j = 1; j <= np; j++) row = row OFS gp[1, j]; print c, k row > (W "/" R "/pca/" c ".tsv") }
      }
    }
    { key = $1 OFS $2 OFS $3 OFS $4
      if (key != k) { flush(); k = key; delete daf; hasg[0] = hasg[1] = 0; gaf[0] = gaf[1] = "."
                      for (j = 1; j <= np; j++) gp[0, j] = gp[1, j] = "." }
      if ($5 < 2) { hasg[$5] = 1; gaf[$5] = mx(gaf[$5], $6); for (j = 1; j <= np; j++) gp[$5, j] = mx(gp[$5, j], $(6 + j)) }
      else { d = $5 - 1; daf[d] = (d in daf) ? mx(daf[d], $6) : $6 } }
    END { flush(); for (x in cnt) { split(x, y, OFS); print substr(x, length(y[1]) + 2), cnt[x] > (W "/" y[1] "/agg/" c ".tsv") } }'
  rm -f keys/[gd]*.$c.txt.s
}
export -f one; export NAMES DGRP DWES PCA_DS DUMP
log "classifying"
echo $CHROMS | tr ' ' '\n' | xargs -P 25 -I{} bash -c 'one {}'

# 3. collect, one folder per region
hdr="chr\tstart\tend\tref\talt"
for r in exonic non_exonic; do
  cd $W/$r
  cat agg/*.tsv 2>/dev/null | awk -F'\t' -v OFS='\t' '{ n[$1 OFS $2 OFS $3 OFS $4 OFS $5] += $6 } END { for (x in n) print x, n[x] }' \
    | sort > counts.body && { printf 'dataset\ttype\tds_class\tgnomad_view\tgnomad_class\tvariants\n'; cat counts.body; } > counts.tsv && rm counts.body
  { printf "$hdr\tdataset\tds_af\tgnomad_af\tgnomad_popmax\tgnomad_matched\n"; for c in $CHROMS; do [ -f markers/$c.tsv ] && cat markers/$c.tsv; done; } | bgzip -@ 16 > markers.tsv.gz
  { printf "$hdr\tmax_af\tdatasets_af_ge_5pct\tgnomad_af\tgnomad_popmax\n"; for c in $CHROMS; do [ -f cand/$c.tsv ] && cat cand/$c.tsv; done; } | bgzip -@ 16 > clinvar_candidates.tsv.gz
  { printf "$hdr"; for d in $PCA_DS; do printf "\t%s" $d; done; for p in $GPOPS; do printf "\tgnomad_%s" $p; done; printf '\n'
    for c in $CHROMS; do [ -f pca/$c.tsv ] && cat pca/$c.tsv; done; } | bgzip -@ 16 > pca_snvs.tsv.gz
  if [ $r = exonic ] && [ $DUMP = 1 ]; then
    { printf "$hdr\tgnomad_af\tdatasets\n"; for c in $CHROMS; do [ -f vars/$c.tsv ] && cat vars/$c.tsv; done; } | bgzip -@ 16 > vars.tsv.gz
  fi
  rm -rf agg markers cand pca vars
  log "$r: $(awk -F'\t' 'NR > 1 && $4 == "overall" {n += $6} END {print n}' counts.tsv) dataset alleles classified"
done
cd $W; rm -rf keys tmp
log "STUDY_DONE"
