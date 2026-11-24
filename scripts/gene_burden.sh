#!/bin/bash
# Gene-level variant burden, exonic regions (gnomAD v4 exome calling intervals).
# Inputs (read-only):
#   $G/work/exonic/vars.tsv.gz     every exonic variant: chr start end ref alt gnomad_af datasets(name=AF;...)
#                                  (freq_study.sh with DUMP=1)
#   $G/annot/snv.<chr>.tsv         VEP of every possible exonic SNV (Naser_VEP CADD exon run; snv_annot.sh)
#   VEP 116 docker image + a VEP cache (mounted read-only) for the other variants. The published run used the
#   HUVYA merged cache; the default below is the Naser_VEP RefSeq cache (results can differ slightly in gene picks).
# Output: $G/gene_counts.tsv  gene  symbol  biotype  dataset  class  af  variants
#           class: all | protein_altering | synonymous ;  af: all | common (AF >= 1%)
#           dataset "gnomad41" = gnomAD v4.1 exome (genome where exome has no AF) = the reference
#         $G/gene_length.tsv   gene  symbol  exonic_positions (positions VEP assigns to the gene, exons +-50 bp)
set -uo pipefail
G=${G:-/home/NASER/annotation/population_study/data/work/gene_burden}
VEP_CACHE=${VEP_CACHE:-/home/NASER_FILES/Naser_VEP/vep_cache}     # directory holding homo_sapiens_<type>/116_GRCh38
VEP_CACHE_OPT=${VEP_CACHE_OPT:---refseq}                          # --refseq | --merged | "" (Ensembl)
VEP_FASTA=${VEP_FASTA:-/home/NASER_FILES/annovar/hg38/reference/hg38.fa}
VEP_MOUNT=${VEP_MOUNT:-/home/NASER_FILES}                          # mounted read-only into the container
CPUS=${CPUS:-120}
export LC_ALL=C
log() { echo "[$(date +%F' '%T)] $*"; }
cd $G; mkdir -p split vep
CHROMS="$(seq 1 22 | tr '\n' ' ')X Y"

# 1. split the exonic variant list per chromosome: SNVs vs everything else
log "splitting variants"
zcat work/exonic/vars.tsv.gz | awk -F'\t' -v OFS='\t' -v d=$G/split 'NR > 1 {
    snv = (length($4) == 1 && length($5) == 1 && $4 != "-" && $5 != "-")
    print > (d "/" (snv ? "snv" : "nonsnv") "." $1 ".tsv") }'

# 2. SNVs: look up the VEP annotation; SNVs missing from the table go to VEP with the rest
one() {   # one <chrom>
  local c=$1
  [ -f split/snv.$c.tsv ] || return 0
  awk -F'\t' -v OFS='\t' -v miss=split/snvmiss.$c.tsv '
    NR == FNR { v[$2 ":" $4 ":" $5] = $0; next }
    { k = $1 ":" $2 ":" $3; if (k in v) { print v[k], $4, $5, $6, $7, $8, $9; delete v[k] } }
    END { for (k in v) print v[k] > miss }' split/snv.$c.tsv annot/snv.$c.tsv > split/snvannot.$c.tsv
}
export -f one
log "SNV lookup"
echo $CHROMS | tr ' ' '\n' | xargs -P 24 -I{} bash -c 'one {}'
# indels/MNVs + SNVs missing from the lookup table (distinct file names: a glob must not pick up other files)
cat split/nonsnv.*.tsv split/snvmiss.*.tsv 2>/dev/null > vep/input.tsv
log "variants: $(cat split/snv.*.tsv split/nonsnv.*.tsv | grep -c .) | SNVs annotated by lookup: $(cat split/snvannot.*.tsv | wc -l) | sent to VEP: $(wc -l < vep/input.tsv)"

# 3. VEP for indels / MNVs / unmatched SNVs, in VEP's default (ensembl) format: chr start end REF/ALT strand id
#    ANNOVAR insertion (start = end = base before it) -> VEP start = end + 1
awk -F'\t' -v OFS='\t' '{ s = $2; e = $3; if ($4 == "-") { s = $2 + 1; e = $2 }
                          print "chr" $1, s, e, $4 "/" $5, "+", $1 ":" $2 ":" $3 ":" $4 ":" $5 }' vep/input.tsv \
  | sort -u | sort -k1,1V -k2,2n > vep/input.ens     # unique whole lines (several alleles can share a position)
rm -f vep/chunk_*; split -n l/$((CPUS / 8)) -d -a 3 vep/input.ens vep/chunk_
log "VEP on $(wc -l < vep/input.ens) variants, $((CPUS / 8)) x 8 forks"
docker run --rm --cpus $CPUS -u "$(id -u):$(id -g)" -e HOME=/tmp -v $VEP_MOUNT:$VEP_MOUNT:ro -v $G:$G -w $G/vep huvya/vep:116 bash -c "
  ls chunk_??? | xargs -P $((CPUS / 8)) -I{} vep --input_file {} --output_file {}.out --format ensembl --tab --force_overwrite \
    --fields Uploaded_variation,Consequence,IMPACT,SYMBOL,Gene,BIOTYPE,MANE_SELECT \
    --cache --offline --dir_cache $VEP_CACHE $VEP_CACHE_OPT --species homo_sapiens --assembly GRCh38 --cache_version 116 \
    --fasta $VEP_FASTA --symbol --mane --canonical --appris --tsl --ccds --biotype \
    --pick_allele --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq \
    --fork 8 --no_stats --buffer_size 5000" > vep/vep.log 2>&1
n_out=$(cat vep/chunk_???.out 2>/dev/null | grep -vc '^#')
log "VEP rows: $n_out"
# attach VEP annotation to the variant rows: key -> consequence impact symbol gene biotype mane
cat vep/chunk_???.out | awk -F'\t' -v OFS='\t' '!/^#/ { print $1, $2, $3, $4, $5, $6, $7 }' > vep/annot.tsv
awk -F'\t' -v OFS='\t' 'NR == FNR { a[$1] = $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7; next }
  { k = $1 ":" $2 ":" $3 ":" $4 ":" $5; if (k in a) print $0, a[k]; else nomatch++ }
  END { if (nomatch) print "variants without VEP annotation: " nomatch > "/dev/stderr" }' vep/annot.tsv vep/input.tsv > split/vep.annot.tsv

# 4. per gene x dataset counts (protein-coding and other genes; one gene per variant = VEP pick)
#    variant row: chr start end ref alt gnomad_af datasets consequence impact symbol gene biotype mane
log "counting"
cat split/snvannot.*.tsv split/vep.annot.tsv | awk -F'\t' -v OFS='\t' '
  BEGIN { split("missense_variant stop_gained stop_lost start_lost frameshift_variant inframe_insertion inframe_deletion splice_acceptor_variant splice_donor_variant protein_altering_variant", P, " ")
          for (i in P) PA[P[i]] = 1 }
  function add(ds, af,    k) {
    for (k = 1; k <= nc; k++) { n[gene OFS ds OFS cl[k] OFS "all"]++; if (af != "." && af != "absent" && af + 0 >= 0.01) n[gene OFS ds OFS cl[k] OFS "common"]++ } }
  { if ($11 == "" || $11 == "-" || $11 == ".") next
    gene = $11 OFS ($10 == "" || $10 == "-" ? "." : $10) OFS $12
    nt = split($8, cq, "&|,"); pa = 0; syn = 0
    for (j = 1; j <= nt; j++) { if (cq[j] in PA) pa = 1; if (cq[j] == "synonymous_variant") syn = 1 }
    nc = 1; cl[1] = "all"; if (pa) cl[++nc] = "protein_altering"; else if (syn) cl[++nc] = "synonymous"
    if ($6 != "absent") add("gnomad41", $6)
    m = split($7, ds, ";"); for (j = 1; j <= m; j++) { split(ds[j], kv, "="); if (kv[1] != "") add(kv[1], kv[2]) } }
  END { for (x in n) print x, n[x] }' | sort > gene_counts.body
{ printf 'gene\tsymbol\tbiotype\tdataset\tclass\taf\tvariants\n'; cat gene_counts.body; } > gene_counts.tsv && rm gene_counts.body

# 5. gene length: exonic positions (exons +-50 bp) VEP assigns to the gene, from the all-SNV table
log "gene lengths"
for c in $CHROMS; do [ -f annot/snv.$c.tsv ] && awk -F'\t' -v OFS='\t' -v c=$c '$7 != "" && $7 != "-" { k = $7 OFS $6 OFS c ":" $1; if (!(k in s)) { s[k] = 1; print $7, $6 } }' annot/snv.$c.tsv; done \
  | awk -F'\t' -v OFS='\t' '{ n[$1 OFS $2]++ } END { for (x in n) print x, n[x] }' | sort > gene_length.body
{ printf 'gene\tsymbol\texonic_positions\n'; cat gene_length.body; } > gene_length.tsv && rm gene_length.body
log "genes counted: $(cut -f1 gene_counts.tsv | sort -u | wc -l), with length: $(($(wc -l < gene_length.tsv) - 1))"
log "GENE_DONE"
