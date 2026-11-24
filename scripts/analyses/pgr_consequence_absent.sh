#!/bin/bash
# Share of each dataset exonic coding variants absent from gnomAD v4.1, by consequence (as Table 1 of the PGR paper)
$(command -v pigz >/dev/null && echo "pigz -dc" || echo "gzip -dc") ${K:-/home/NASER/annotation/population_study/data/archive}/gene_level/exonic_variants_annotated.tsv.gz | awk -F"\t" -v OFS="\t" "
  NR == 1 { next }
  \$12 != \"protein_coding\" { next }
  { c = \$8; k = \"\"
    if (c ~ /stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant/) k = \"pLoF\"
    else if (c ~ /missense_variant/) k = \"missense\"
    else if (c ~ /inframe_insertion|inframe_deletion/) k = \"inframe\"
    else if (c ~ /synonymous_variant/) k = \"synonymous\"
    if (k == \"\") next
    ab = (\$6 == \"absent\")
    n = split(\$7, d, \";\"); for (i = 1; i <= n; i++) { split(d[i], kv, \"=\"); if (kv[1] == \"\") continue
      t[kv[1] OFS k]++; if (ab) a[kv[1] OFS k]++ } }
  END { print \"dataset\", \"class\", \"variants\", \"absent_from_gnomad\", \"pct_absent\"
        for (x in t) printf \"%s\t%d\t%d\t%.2f\n\", x, t[x], a[x], 100 * a[x] / t[x] }" > ${OUT:-tables}/coding_variants_absent_from_gnomad_by_consequence.tsv
echo CONS_DONE
