#!/bin/bash
# Build reports/population_frequency_report.md (the whole study in one file) from the README key findings and
# reports/sections/*.md, then the HTML versions of the Markdown reports with pandoc.
# usage: bash scripts/build_report_md.sh   (from the repository root; run make_results_md.py first)
set -euo pipefail
S=reports/sections
OUT=reports/population_frequency_report.md
{
  echo "# Beyond gnomAD: population allele-frequency resources compared with gnomAD v4.1"
  echo
  echo "Full study report. Sections are also available separately in [sections/](sections/)."
  echo
  sed -n '/^## Key findings/,/^## Documents/p' README.md | sed '$d' \
    | sed -e 's#](reports/sections/#](sections/#g' -e 's#](reports/#](#g' -e 's#](results/#](../results/#g' -e 's#](scripts/#](../scripts/#g'
  for f in $S/01_materials_and_methods.md $S/02_results.md $S/03_data_quality_audit.md \
           $S/07_additional_analyses.md $S/08_comparison_pakistan_genome_resource.md \
           $S/04_limitations_and_publication.md $S/05_further_analyses.md $S/06_data_manifest.md $S/09_reproducibility.md; do
    echo
    # one level down so the report keeps a single top-level title; links are rewritten for reports/
    sed -e 's/^#/##/' -e 's#](\.\./\.\./#](../#g' -e 's#](\([0-9][0-9]_[a-z_]*\.md\)#](sections/\1#g' "$f"
  done
} > $OUT
echo "wrote $OUT ($(wc -l < $OUT) lines)"
# Markdown -> HTML (pandoc), both reports; links to results/ and sections/ stay live inside the repository
if command -v pandoc >/dev/null; then
  pandoc reports/publication_report.md -s --toc --toc-depth=2 -c publication.css \
    -o reports/publication_report.html && echo "wrote reports/publication_report.html"
  pandoc reports/population_frequency_report.md -s --toc --toc-depth=2 -c publication.css \
    --metadata pagetitle="Population frequency report" \
    -o reports/population_frequency_report.html && echo "wrote reports/population_frequency_report.html"
fi
# single-file copies (figures embedded) for sharing outside the repository: set SHARE_DIR to build them
if [ -n "${SHARE_DIR:-}" ]; then
  mkdir -p "$SHARE_DIR"
  python3 scripts/build_artifact_html.py reports/publication_report.md "Population Frequency Study" "$SHARE_DIR/publication_report.html"
  python3 scripts/build_artifact_html.py reports/population_frequency_report.md "Population Frequency Report" "$SHARE_DIR/population_frequency_report.html"
fi
