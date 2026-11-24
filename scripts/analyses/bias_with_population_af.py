#!/usr/bin/env python3
"""Run BIAS-2015 unchanged, but let population datasets beyond gnomAD raise a variant's frequency.

usage: bias_with_population_af.py <BIAS dir> <vep.json> <required_paths.json> <out.tsv> <lookup.tsv|none> <column>

BIAS's frequency criteria (BA1, BS1, PM2 and expert-panel rules) read the gnomAD record that
extract_from_vep_json.process_variant builds (allAf, per-group AFs, popmax). This wrapper wraps that function:
after BIAS builds the variant, allAf and popmax are raised to the population-dataset AF from <lookup.tsv> when it is
higher. If gnomAD has no record, one is created, so "absent from controls" (PM2) no longer applies.
<lookup.tsv>: chr start end ref alt <columns...> in ANNOVAR notation (as in the frequency study); <column> picks the AF.
With "none" BIAS runs exactly as shipped (scenario A). The BIAS source files are not modified.
"""
import csv
import os
import sys

bias_dir, vep_json, paths, out, lookup_fp = sys.argv[1:6]
column = sys.argv[6] if len(sys.argv) > 6 else None

lookup = {}
if lookup_fp != "none":
    with open(lookup_fp) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            v = row.get(column, "")
            if v not in ("", ".", "NA"):
                lookup[(row["chr"], row["start"], row["end"], row["ref"], row["alt"])] = float(v)

os.chdir(bias_dir)
sys.path.insert(0, bias_dir)
from src.bias_2015 import extract_from_vep_json as X  # noqa: E402
import bias_2015  # noqa: E402

stats = {"seen": 0, "matched": 0, "raised": 0, "created": 0}
_orig = X.process_variant


def annovar_key(entry):
    """VEP JSON entry -> ANNOVAR-style key used by the frequency study."""
    chrom = str(entry.get("seq_region_name", "")).replace("chr", "")
    chrom = "M" if chrom == "MT" else chrom
    a = entry.get("allele_string", "").split("/")
    if len(a) != 2:
        return None
    ref, alt = a
    start, end = int(entry["start"]), int(entry["end"])
    if ref == "-":                      # VEP insertion: start = end + 1; ANNOVAR anchors on the base before
        return (chrom, str(end), str(end), "-", alt)
    return (chrom, str(start), str(end), ref, alt)


def process_variant(vep_entry, hgnc_to_gene_data):
    v = _orig(vep_entry, hgnc_to_gene_data)
    if v is None or not lookup:
        return v
    stats["seen"] += 1
    k = annovar_key(vep_entry)
    af = lookup.get(k) if k else None
    if af is None:
        return v
    stats["matched"] += 1
    g = v.gnomad if v.gnomad else {}
    if not g:
        stats["created"] += 1
    if af > float(g.get("allAf", 0) or 0) or af > float(g.get("popmax", 0) or 0):
        stats["raised"] += 1
        g["allAf"] = max(af, float(g.get("allAf", 0) or 0))
        g["popmax"] = max(af, float(g.get("popmax", 0) or 0))
    v.gnomad = g
    return v


X.process_variant = process_variant
sys.argv = ["bias_2015.py", vep_json, paths, out, "--annotator", "vep"]
bias_2015.main()
if lookup:
    print(f"population AF ({column}): variants {stats['seen']}, in lookup {stats['matched']}, "
          f"AF raised {stats['raised']}, gnomAD record created {stats['created']}", file=sys.stderr)
