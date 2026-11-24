#!/bin/bash
# compact SNV annotation from Naser_VEP CADD exon run (read-only): per chromosome, sorted by position
# columns: pos ref alt consequence impact symbol gene biotype mane_select
set -uo pipefail
T=${CADD_VEP_TSV:-/home/NASER_FILES/Naser_VEP/results/CADD_exons_merged_64cpu/hg38_cadd17_exons_pad50.vep.tsv.gz}
O=${G:-/home/NASER/annotation/population_study/data/work/gene_burden}/annot
cols=$(zcat $T | head -1 | awk -F"\t" "{for(j=1;j<=NF;j++) c[\$j]=j; print c[\"CHROM\"], c[\"POS\"], c[\"REF\"], c[\"ALT\"], c[\"Consequence\"], c[\"IMPACT\"], c[\"SYMBOL\"], c[\"Gene\"], c[\"BIOTYPE\"], c[\"MANE_SELECT\"]}")
echo "columns: $cols"
zcat $T | awk -F"\t" -v OFS="\t" -v cols="$cols" -v O=$O "BEGIN{split(cols,C,\" \")} NR>1 { c=\$C[1]; sub(/^chr/,\"\",c); if (c==\"MT\") c=\"M\"; print \$C[2], \$C[3], \$C[4], \$C[5], \$C[6], \$C[7], \$C[8], \$C[9], \$C[10] > (O \"/snv.\" c \".tsv\") }"
ls $O | wc -l; echo SNV_ANNOT_DONE
