#!/bin/bash
# Shared input for the founder-variant and carrier-frequency analyses (no raw frequency files needed).
# ClinVar P/LP alleles with
#   - AF in every population dataset (exonic variants; data/archive/gene_level/exonic_variants_annotated.tsv.gz)
#   - gnomAD v4.1 exome and genome AF per genetic-ancestry group, from the Naser_VEP VEP 116 RefSeq cache
#     (--af_gnomade --af_gnomadg). Nothing under /home/HUVYA is used.
# Output: $W/shared/clinvar_plp_population_af.tsv.gz
set -uo pipefail
K=${K:-/home/NASER/annotation/population_study/data/archive}
W=${W:-/home/NASER/annotation/population_study/data/work/extra_analyses}
CV=${CV:-$K/clinvar_20261115.vcf.gz}
V=${V_ROOT:-/home/NASER_FILES/Naser_VEP}      # VEP 116 RefSeq cache with gnomAD v4.1 (read-only)
FA=${FA:-/home/NASER_FILES/annovar/hg38/reference/hg38.fa}   # the RefSeq cache is BAM-edited: VEP needs a FASTA
export LC_ALL=C
log() { echo "[$(date +%F' '%T)] $*"; }
mkdir -p $W/shared/vep; cd $W/shared

# 1. ClinVar P/LP records (CLNSIG starting with Pathogenic or Likely_pathogenic, incl. low-penetrance / risk suffixes)
bcftools view -i 'INFO/CLNSIG ~ "^Pathogenic" || INFO/CLNSIG ~ "^Likely_pathogenic"' $CV -Oz -o clinvar_plp.vcf.gz && tabix -f -p vcf clinvar_plp.vcf.gz
log "ClinVar P/LP records: $(bcftools view -H clinvar_plp.vcf.gz | wc -l)"

# 2. normalised keys (same normaliser as freq_study.sh) + ClinVar fields
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\t%INFO/MC\t%INFO/CLNDN\n' clinvar_plp.vcf.gz \
| awk -F'\t' -v OFS='\t' '
  function key(c, s, r, a,    ins, k, e) {
    r = toupper(r); a = toupper(a)
    if (r !~ /^[ACGTN]+$/ || a !~ /^[ACGTN]+$/ || r == a) return ""
    k = 0; while (k < length(r) && k < length(a) && substr(r, k+1, 1) == substr(a, k+1, 1)) k++
    r = substr(r, k+1); a = substr(a, k+1); s = s + k
    while (length(r) > 0 && length(a) > 0 && substr(r, length(r), 1) == substr(a, length(a), 1)) {
      r = substr(r, 1, length(r)-1); a = substr(a, 1, length(a)-1) }
    if (r == "" && a == "") return ""
    if (r == "") { s = s - 1; e = s } else e = s + length(r) - 1
    if (r == "") r = "-"; if (a == "") a = "-"
    return c OFS s OFS e OFS r OFS a
  }
  function stars(x) { return x ~ /practice_guideline/ ? 4 : x ~ /reviewed_by_expert_panel/ ? 3 :
                             x ~ /multiple_submitters,_no_conflicts/ ? 2 : x ~ /criteria_provided/ ? 1 : 0 }
  { c = $1; sub(/^chr/, "", c); if (c == "MT") c = "M"; k = key(c, $2 + 0, $4, $5); if (k == "") next
    g = $8; sub(/:.*/, "", g); mc = $9; sub(/^[^|]*\|/, "", mc); sub(/,.*/, "", mc)
    print k, $3, $6, stars($7), (g == "." ? "" : g), (mc == "." ? "" : mc), $10 }' > clinvar_plp.keys.tsv
log "P/LP alleles with a key: $(wc -l < clinvar_plp.keys.tsv)"

# 3. AF in every population dataset: join with the kept exonic variant table
$(command -v pigz >/dev/null && echo "pigz -dc" || echo "gzip -dc") $K/gene_level/exonic_variants_annotated.tsv.gz | awk -F'\t' -v OFS='\t' '
  NR == FNR { cv[$1 OFS $2 OFS $3 OFS $4 OFS $5] = 1; next }
  FNR == 1 { next }
  ($1 OFS $2 OFS $3 OFS $4 OFS $5) in cv { print $1, $2, $3, $4, $5, $6, $7, $8, $10 }' clinvar_plp.keys.tsv - > clinvar_plp.dataset_af.tsv
log "P/LP alleles found among exonic population variants: $(wc -l < clinvar_plp.dataset_af.tsv)"

# 4. gnomAD v4.1 per-group AF from the VEP cache (exome: AFR AMR ASJ EAS FIN MID NFE SAS; genome: + AMI)
# VEP output names: gnomADe_AF and gnomADe_<GROUP>_AF (likewise gnomADg)
F="Uploaded_variation,gnomADe_AF$(for g in AFR AMR ASJ EAS FIN MID NFE SAS; do printf ',gnomADe_%s_AF' $g; done),gnomADg_AF$(for g in AFR AMR ASJ EAS FIN MID NFE SAS; do printf ',gnomADg_%s_AF' $g; done)"
docker run --rm --cpus 64 -u "$(id -u):$(id -g)" -e HOME=/tmp -v /home/NASER_FILES:/home/NASER_FILES:ro -v $W:$W -w $W/shared/vep huvya/vep:116 \
  vep -i $W/shared/clinvar_plp.vcf.gz -o $W/shared/vep/clinvar_plp.gnomad.tsv --tab --force_overwrite --no_stats --offline --cache \
      --dir_cache $V/vep_cache --refseq --cache_version 116 --assembly GRCh38 --fasta $FA --check_existing --af_gnomade --af_gnomadg \
      --pick --fork 32 --buffer_size 5000 --fields "$F" > $W/shared/vep/vep.log 2>&1
[ -s vep/clinvar_plp.gnomad.tsv ] || { log "VEP failed, see vep/vep.log"; exit 1; }
log "VEP rows: $(grep -vc '^#' vep/clinvar_plp.gnomad.tsv)"

# 5. one table: ClinVar fields, VEP gnomAD groups, dataset AFs (exonic only; empty when not exonic / not observed)
${PY:-/root/miniforge3/envs/ahm-env/bin/python} - <<'EOF'
import pandas as pd
k = pd.read_csv("clinvar_plp.keys.tsv", sep="\t", header=None, dtype=str,
                names=["chr", "start", "end", "ref", "alt", "clinvar_id", "clnsig", "stars", "gene", "consequence", "condition"])
v = pd.read_csv("vep/clinvar_plp.gnomad.tsv", sep="\t", comment="#", dtype=str, header=None)
v.columns = ["clinvar_id", "gnomade_af", "e_afr", "e_amr", "e_asj", "e_eas", "e_fin", "e_mid", "e_nfe", "e_sas",
             "gnomadg_af", "g_afr", "g_amr", "g_asj", "g_eas", "g_fin", "g_mid", "g_nfe", "g_sas"]
v = v.drop_duplicates("clinvar_id").replace("-", pd.NA)
d = pd.read_csv("clinvar_plp.dataset_af.tsv", sep="\t", header=None, dtype=str,
                names=["chr", "start", "end", "ref", "alt", "gnomad_ref_af", "datasets", "vep_consequence", "vep_gene"])
m = k.merge(v, on="clinvar_id", how="left").merge(d, on=["chr", "start", "end", "ref", "alt"], how="left")
for g in ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]:   # group AF: exome, genome where exome is missing
    e, gg = pd.to_numeric(m["e_" + g], errors="coerce"), pd.to_numeric(m["g_" + g], errors="coerce")
    m["gnomad_" + g] = e.where(e.notna(), gg)
grp = [c for c in m.columns if c.startswith("gnomad_") and c[7:] in ("afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas")]
m["gnomad_popmax"] = m[grp].max(axis=1)
m["gnomad_overall"] = pd.to_numeric(m.gnomade_af, errors="coerce").fillna(pd.to_numeric(m.gnomadg_af, errors="coerce"))
m["exonic_in_study"] = m.datasets.notna() | m.gnomad_ref_af.notna()
keep = ["chr", "start", "end", "ref", "alt", "clinvar_id", "clnsig", "stars", "gene", "consequence", "condition",
        "gnomad_overall", "gnomad_popmax"] + grp + ["exonic_in_study", "gnomad_ref_af", "datasets"]
m[keep].to_csv("clinvar_plp_population_af.tsv.gz", sep="\t", index=False, compression="gzip")
print(f"rows {len(m)}; exonic in study {int(m.exonic_in_study.sum())}; with gnomAD AF {int(m.gnomad_overall.notna().sum())}")
EOF
log "PREP_DONE"
