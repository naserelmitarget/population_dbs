#!/bin/bash
# ClinVar angle: ClinVar variants that are >= 5% in at least one population dataset (freq_study.sh
# clinvar_candidates.tsv.gz), with their classification, review stars and gnomAD v4.1 frequency.
# Same normalisation as freq_union.sh / freq_study.sh. Read-only on inputs.
set -uo pipefail
CV=${CV:-/home/NASER/annotation/population_study/data/raw/clinvar/clinvar_20261115.vcf.gz}
W=${W:-/home/NASER/annotation/population_study/data/work/freq_study}
export LC_ALL=C
cd $W
{ zcat $CV; echo "#CANDIDATES"; zcat clinvar_candidates.tsv.gz; } \
| awk -F'\t' -v OFS='\t' '
  function key(c, s, r, a,    ins, k, e) {
    r = toupper(r); a = toupper(a)
    if (r !~ /^([ACGTN]+|-)$/ || a !~ /^([ACGTN]+|-)$/ || r == a) return ""
    ins = (r == "-"); if (r == "-") r = ""; if (a == "-") a = ""
    k = 0; while (k < length(r) && k < length(a) && substr(r, k+1, 1) == substr(a, k+1, 1)) k++
    r = substr(r, k+1); a = substr(a, k+1); s = s + k
    while (length(r) > 0 && length(a) > 0 && substr(r, length(r), 1) == substr(a, length(a), 1)) {
      r = substr(r, 1, length(r)-1); a = substr(a, 1, length(a)-1) }
    if (r == "" && a == "") return ""
    if (r == "") { if (!ins) s = s - 1; e = s } else e = s + length(r) - 1
    if (r == "") r = "-"; if (a == "") a = "-"
    return c OFS s OFS e OFS r OFS a
  }
  function grp(x) {
    if (x ~ /^(Pathogenic|Likely_pathogenic|Pathogenic\/Likely_pathogenic)/) return "P/LP"
    if (x ~ /^(Uncertain_significance|VUS)/) return "VUS"
    if (x ~ /^Conflicting/) return "Conflicting"
    if (x ~ /^(Benign|Likely_benign|Benign\/Likely_benign)/) return "B/LB"
    return "Other"
  }
  function stars(x) {
    if (x ~ /practice_guideline/) return 4
    if (x ~ /reviewed_by_expert_panel/) return 3
    if (x ~ /multiple_submitters,_no_conflicts/) return 2
    if (x ~ /criteria_provided/) return 1
    return 0
  }
  $0 == "#CANDIDATES" { cand = 1; next }
  !cand && /^#/ { next }
  !cand {
    c = $1; sub(/^chr/, "", c); if (c == "MT") c = "M"
    delete info; m = split($8, kv, ";")
    for (j = 1; j <= m; j++) { p = index(kv[j], "="); if (p) info[substr(kv[j], 1, p-1)] = substr(kv[j], p+1) }
    sig = ("CLNSIG" in info) ? info["CLNSIG"] : "no_classification"
    g = info["GENEINFO"]; sub(/:.*/, "", g); gsub(/\|.*/, "", g)
    mc = info["MC"]; sub(/^[^|]*\|/, "", mc); sub(/,.*/, "", mc)
    n = split($5, alts, ",")
    for (j = 1; j <= n; j++) { k = key(c, $2 + 0, $4, alts[j]); if (k == "") continue
      cv[k] = $3 OFS grp(sig) OFS sig OFS stars(info["CLNREVSTAT"]) OFS (g == "" ? "." : g) OFS (mc == "" ? "." : mc) OFS info["CLNDN"]; ncv++ }
    next
  }
  FNR == 1 || $1 == "chr" { print "chr", "start", "end", "ref", "alt", "clinvar_id", "clnsig_group", "clnsig", "review_stars",
                            "gene", "consequence", "condition", "max_af", "datasets_af_ge_5pct", "gnomad_af", "gnomad_popmax"; next }
  { k = $1 OFS $2 OFS $3 OFS $4 OFS $5; if (k in cv) { print k, cv[k], $6, $7, $8, $9; hit++ } }
  END { printf "clinvar alleles indexed: %d; matched candidates: %d\n", ncv, hit > "/dev/stderr" }' \
| bgzip > clinvar_hits.tsv.gz
echo CLINVAR_DONE
