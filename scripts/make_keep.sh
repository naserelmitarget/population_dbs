#!/bin/bash
# Compact copy of everything derived from the raw inputs (the "kept archive"). Copies and compresses; deletes nothing.
# With this archive every table, figure and additional analysis can be regenerated without the raw frequency files.
# Paths from scripts/config.sh (or the environment).
set -uo pipefail
HERE=$(cd $(dirname $0) && pwd)
[ -n "${KEEP:-}" ] || source $HERE/config.sh
log() { echo "[$(date +%F' '%T)] $*"; }
mkdir -p $KEEP/union $KEEP/study/exonic $KEEP/study/non_exonic $KEEP/gene_level $KEEP/provenance

# 1. union: bgzip + tabix (ANNOVAR-style chr, start, end; header line starts with #)
bgzip -@ 48 -c $UNION_DIR/hg38_freq_union.txt > $KEEP/union/hg38_freq_union.txt.gz \
  && tabix -f -s1 -b2 -e3 -c "#" $KEEP/union/hg38_freq_union.txt.gz && log "union compressed"
cp $UNION_DIR/summary_*.tsv $UNION_DIR/run.log $KEEP/union/

# 2. study: per-variant outputs, tables and figures, both regions
for r in exonic non_exonic; do
  cp $STUDY_DIR/$r/*.tsv.gz $STUDY_DIR/$r/counts.tsv $KEEP/study/$r/
  cp -r $STUDY_DIR/$r/tables $STUDY_DIR/$r/figures $KEEP/study/$r/
done
cp $STUDY_DIR/extraction_stats.tsv $STUDY_DIR/run.log $KEEP/study/
cp -r $STUDY_DIR/regions $KEEP/study/

# 3. gene level: counts, lengths, tables, figures and the annotated exonic variant list
cp $GENE_DIR/gene_length.tsv $GENE_DIR/gene_burden.log $KEEP/gene_level/
bgzip -@ 16 -c $GENE_DIR/gene_counts.tsv > $KEEP/gene_level/gene_counts.tsv.gz
cp -r $GENE_DIR/tables $GENE_DIR/figures $KEEP/gene_level/
{ printf "chr\tstart\tend\tref\talt\tgnomad_af\tdatasets\tconsequence\timpact\tsymbol\tgene\tbiotype\tmane_select\n"
  cat $GENE_DIR/split/snvannot.*.tsv $GENE_DIR/split/vep.annot.tsv; } | bgzip -@ 48 > $KEEP/gene_level/exonic_variants_annotated.tsv.gz \
  && log "annotated exonic variants compressed"

# 4. ClinVar release used
cp $CLINVAR* $KEEP/

# 5. provenance: checksums of the raw inputs (while they exist) and of the archive
( cd $(dirname $ANNOVAR_TABLES) && find $(basename $ANNOVAR_TABLES) $(basename $UCSC_DIR) -type f \( -name "*.txt" -o -name "*.idx" -o -name "*.vcf.gz" \
    -o -name "*.vcf.bgz" -o -name "*.tbi" -o -name "*.bb" -o -name "*.tsv.gz" \) 2>/dev/null | sort | xargs -r -P 16 -I{} md5sum {} | sort -k2 ) \
  > $KEEP/provenance/raw_inputs.md5
( cd $KEEP && find . -type f ! -path "./provenance/*" | sort | xargs -P 16 -I{} md5sum {} | sort -k2 ) > $KEEP/provenance/kept_archive.md5
du -sh $KEEP/*; log KEEP_DONE
