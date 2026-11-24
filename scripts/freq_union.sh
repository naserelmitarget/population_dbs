#!/bin/bash
# Union of ANNOVAR frequency databases (read-only on the databases).
# Output: one ANNOVAR filter table, sorted by chromosome (1..22, X, Y, M) and position:
#   #Chr  Start  End  Ref  Alt  is_in  how_many
#   is_in    = datasets holding the variant, ';'-separated, in the order of DBS below
#   how_many = number of those datasets
# Variants are assumed left-aligned. Every row is only trimmed to ANNOVAR notation, so the same
# variant written as "-" (ANNOVAR) or anchored (VCF) collapses to one row:
#   chromosome: "chr" removed, MT -> M
#   shared leading bases trimmed, then shared trailing bases; empty allele -> "-"
#   insertion Start = End = the base before it (ANNOVAR convention)
#   rows whose alleles are not A/C/G/T/N or "-" are dropped (counted in keys/*.stats)
# 1000 Genomes hg38 files write the alt as <ref length><alt bases> ("0C", "1AC", "21"); decoded here.
# VCF inputs: multi-allelic ALTs are split; symbolic/"*" alleles dropped; every listed allele counts
# (no AC/FILTER filtering, same as the ANNOVAR tables). Only chromosomes 1-22, X, Y, M are kept.
# UCSC varFreqsAll (bb_to_table.sh output): one pass; each cohort code in the "sources" column becomes its
# own dataset (codes not listed in DBS, i.e. aou/abr which we take from ANNOVAR, are ignored).
# ALFA (vcf_alfa): UCSC's ALFA VCF has indel POS one base too far right (REF matches hg38 at POS-1
# for 100% of indels on chr22, at POS for 5%); indel POS is shifted -1. SNVs are correct.
set -uo pipefail
A=${A:-/home/NASER/annotation/population_study/data/raw/annovar_popfreq}           # ANNOVAR tables (moved from /root)
H=${H:-$A}                                                # allofus table (same folder as the others by default)
U=${U:-/home/NASER/annotation/population_study/data/raw/ucsc_varfreqs}              # UCSC varFreqs VCFs
W=${W:-/home/NASER/annotation/population_study/data/work/freq_union}
OUT=${OUT:-hg38_freq_union.txt}
BB=${BB:-$U/varFreqsAll/varFreqsAll.tsv.gz}
# name:file[,file...]:format[:bb source code]   (format: annovar | 1000g | vcf | bb)
DBS=(gnomad41_exome:$A/hg38_gnomad41_exome.txt:annovar
     gnomad41_genome:$A/hg38_gnomad41_genome.txt:annovar
     gnomad_exome:$A/hg38_gnomad_exome.txt:annovar
     gnomad_genome:$A/hg38_gnomad_genome.txt:annovar
     kaviar_20150923:$A/hg38_kaviar_20150923.txt:annovar
     exac03:$A/hg38_exac03.txt:annovar
     1000g2015aug:$A/hg38_ALL.sites.2015_08.txt:1000g
     gme:$A/hg38_gme.txt:annovar
     allofus:$A/hg38_allofus.txt:annovar
     spark_wes:$BB:bb:spk
     sfari_wgs:$BB:bb:wgs
     genomeasia:$BB:bb:ga
     kova:$BB:bb:kov
     tommo:$BB:bb:tom
     mcps:$BB:bb:mcp
     indigen:$BB:bb:ind
     finngen:$BB:bb:fin
     saudi:$BB:bb:sau
     swegen:$BB:bb:swe
     topmed:$BB:bb:top
     alfa:$U/alfa/ALFA.vcf.gz:vcf_alfa
     colorsdb:$U/colorsDb/colorsDbSnv.vcf.gz:vcf
     ga4k:$U/ga4k/ga4kSnv.vcf.gz:vcf
     gonl:$U/gonl/gonl.vcf.gz:vcf
     gregor:$U/gregor/gregor.vcf.gz:vcf
     hgdp1kg:$U/hgdp1kFreq/hgdp1k.freq.vcf.gz:vcf
     hostseq:$U/hostseq/hostseq.vcf.gz:vcf
     sg10k_npm:$U/npm/SG10K_Health_r5.3.2.sites.vcf.bgz:vcf
     schema:$U/schema/SCHEMA_variant_results_withAF.vcf.gz:vcf
     sgdp:$U/sgdpFreq/sgdp.freq.vcf.gz:vcf
     svatalog:$U/svatalog/svatalog.vcf.gz:vcf
     ukbb:$U/ukbb/ukbb.vcf.gz:vcf
     wbbc:$U/wbbc/wbbc.vcf.gz:vcf)
mkdir -p $W/keys $W/chr $W/tmp; cd $W
log() { echo "[$(date +%F' '%T)] $*"; }
export LC_ALL=C
NAMES=$(for d in "${DBS[@]}"; do printf '%s,' "${d%%:*}"; done); NAMES=${NAMES%,}

for d in "${DBS[@]}"; do IFS=: read -r _ fs _ <<< "$d"; for f in ${fs//,/ }; do [ -s $f ] || { log "missing $f"; exit 1; }; done; done

# 1. normalised keys (start, end, ref, alt), one file per database x chromosome
norm() {   # norm <db index, or code=index,... for bb> <file[,file...]> <format>
  local i=$1 fs=$2 fmt=$3 f
  for f in ${fs//,/ }; do case $f in *.gz|*.bgz) gzip -dc $f;; *) cat $f;; esac; done \
  | awk -F'\t' -v OFS='\t' -v dir="$W/keys" -v db=$i -v name="$fs" -v fmt=$fmt '
    BEGIN { if (fmt == "bb") { n = split(db, m, ","); for (j = 1; j <= n; j++) { split(m[j], kv, "="); code[kv[1]] = kv[2] }; db = "bb" } }
    function emit(c, s, r, a, ins, k, e) {      # prints the key to every file in out[]
      r = toupper(r); a = toupper(a)
      if (r !~ /^([ACGTN]+|-)$/ || a !~ /^([ACGTN]+|-)$/ || r == a) { bad++; return }
      ins = (r == "-")                 # ANNOVAR insertion: Start is already the base before it
      if (r == "-") r = ""; if (a == "-") a = ""
      k = 0; while (k < length(r) && k < length(a) && substr(r, k+1, 1) == substr(a, k+1, 1)) k++
      r = substr(r, k+1); a = substr(a, k+1); s = s + k
      while (length(r) > 0 && length(a) > 0 && substr(r, length(r), 1) == substr(a, length(a), 1)) {
        r = substr(r, 1, length(r)-1); a = substr(a, 1, length(a)-1) }
      if (r == "" && a == "") { bad++; return }
      if (r == "") { if (!ins) s = s - 1; e = s } else e = s + length(r) - 1
      if (r == "") r = "-"; if (a == "") a = "-"
      if (fmt == "bb") { for (k in hit) print s, e, r, a > (dir "/" hit[k] "." c ".txt") }
      else print s, e, r, a > (dir "/" db "." c ".txt")
      n_ok++
    }
    fmt == "bb" && /^#/ { sub(/^#/, ""); for (j = 1; j <= NF; j++) col[$j] = j; next }
    /^#/ { next }
    {
      c = $1; sub(/^chr/, "", c); if (c == "MT") c = "M"
      if (c !~ /^([0-9]+|X|Y|M)$/) { other++; next }
      if (fmt == "bb") {               # 0-based chromStart, VCF-style ref/alt, sources = cohort codes
        delete hit; n = split($col["sources"], src, ","); h = 0
        for (j = 1; j <= n; j++) if (src[j] in code) hit[++h] = code[src[j]]
        if (h) emit(c, $col["chromStart"] + 1, $col["ref"], $col["alt"]); else skipped++
      } else if (fmt ~ /^vcf/) {      # CHROM POS ID REF ALT ... (vcf_alfa: indel POS shifted -1)
        n = split($5, alts, ",")
        for (j = 1; j <= n; j++) emit(c, $2 + 0 - (fmt == "vcf_alfa" && length($4) != length(alts[j])), $4, alts[j])
      } else if (fmt == "1000g") {     # chr pos ref <n><alt> af rsid
        if (!match($4, /^[0-9]+/)) { emit(c, $2 + 0, $3, $4); next }   # SNVs have no prefix
        n = substr($4, 1, RLENGTH) + 0; a = substr($4, RLENGTH + 1)
        if (n == 0) r = "-"
        else if (length($3) == n) r = $3
        else { bad++; next }
        emit(c, $2 + 0, r, a == "" ? "-" : a)
      } else emit(c, $2 + 0, $4, $5)
    }
    END { printf "%s\trows_kept\t%d\trows_dropped\t%d\tother_contigs\t%d\tonly_ignored_sources\t%d\n", name, n_ok, bad, other, skipped > (dir "/" db ".stats") }'
}
log "normalising ${#DBS[@]} databases: $NAMES"
bbmap=; bbfile=
for i in "${!DBS[@]}"; do
  IFS=: read -r _ f fmt cd <<< "${DBS[$i]}"
  if [ $fmt = bb ]; then bbmap+="${bbmap:+,}$cd=$i"; bbfile=$f; else norm $i $f $fmt & fi
done
[ -n "$bbmap" ] && norm $bbmap $bbfile bb &
wait
cat keys/*.stats

# 2. per chromosome: sorted unique keys per database, merged; one row per variant with its dataset set
one() {   # one <chrom>
  local c=$1 i
  local -a parts=()
  for i in $(seq 0 $((ND - 1))); do
    [ -f keys/$i.$c.txt ] || continue
    sort -u -S 8G -T tmp -t $'\t' -k1,1n -k2,2n -k3,3 -k4,4 keys/$i.$c.txt \
      | awk -v b=$i '{print $0 "\t" b}' > keys/$i.$c.sorted && rm -f keys/$i.$c.txt
    parts+=(keys/$i.$c.sorted)
  done
  sort -m -T tmp -t $'\t' -k1,1n -k2,2n -k3,3 -k4,4 -k5,5n "${parts[@]}" \
    | awk -F'\t' -v OFS='\t' -v c=$c -v names=$NAMES '
        BEGIN { split(names, nm, ",") }
        function flush() { if (k != "") print c, k, set, cnt }
        { key = $1 "\t" $2 "\t" $3 "\t" $4
          if (key != k) { flush(); k = key; set = ""; cnt = 0 }
          set = set (cnt ? ";" : "") nm[$5 + 1]; cnt++ }
        END { flush() }' > chr/$c.txt
  rm -f "${parts[@]}"
}
export -f one; export W NAMES ND=${#DBS[@]}
chroms=$(ls keys | sed -nE 's/^[0-9]+\.([^.]+)\.txt$/\1/p' | sort -u)
log "union for $(echo $chroms | wc -w) chromosomes"
echo $chroms | tr ' ' '\n' | xargs -P 25 -I{} bash -c 'one {}'

# 3. one file, chromosomes in 1..22, X, Y, M order (anything else after, alphabetically)
order=$(echo $chroms | tr ' ' '\n' | awk '{ k = ($0 ~ /^[0-9]+$/) ? sprintf("%03d", $0) : ($0 == "X" ? "100" : $0 == "Y" ? "101" : $0 == "M" ? "102" : "200" $0); print k "\t" $0 }' | sort | cut -f2)
{ printf '#Chr\tStart\tEnd\tRef\tAlt\tis_in\thow_many\n'; for c in $order; do cat chr/$c.txt; done; } > $OUT
log "wrote $OUT: $(($(wc -l < $OUT) - 1)) variants"
rm -f chr/*.txt

# 4. summary: variants per dataset, per how_many, and per exact dataset combination
awk -F'\t' 'NR > 1 { h[$7]++; combo[$6]++; n = split($6, x, ";"); for (j = 1; j <= n; j++) { per[x[j]]++; if (n == 1) only[x[j]]++ } }
  END {
    print "dataset\tvariants\tonly_in_this_dataset" > "summary_per_dataset.tsv"
    for (d in per) print d "\t" per[d] "\t" only[d] + 0 >> "summary_per_dataset.tsv"
    print "how_many\tvariants" > "summary_how_many.tsv"
    for (v in h) print v "\t" h[v] >> "summary_how_many.tsv"
    print "is_in\tvariants" > "summary_combinations.tsv"
    for (v in combo) print v "\t" combo[v] >> "summary_combinations.tsv" }' $OUT
sort -n summary_how_many.tsv -o summary_how_many.tsv
(head -1 summary_combinations.tsv; tail -n +2 summary_combinations.tsv | sort -t $'\t' -k2,2nr) > tmp/c && mv tmp/c summary_combinations.tsv
column -t -s $'\t' summary_per_dataset.tsv; echo; column -t -s $'\t' summary_how_many.tsv

# 5. ANNOVAR index (1 kb bins)
IDX=${IDX:-$(dirname $(readlink -f $0))/make_annovar_idx.py}
[ -f $IDX ] && python3 $IDX $OUT && log "index written"
mv keys/*.stats . 2>/dev/null; rmdir keys tmp chr 2>/dev/null
log "UNION_DONE"
