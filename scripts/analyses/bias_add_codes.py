#!/usr/bin/env python3
"""Add a readable acmgCodes column (e.g. "PVS1,PM2_Supporting,PP5_VeryStrong") to a BIAS-2015 TSV.

BIAS writes its evidence as JSON in the 'rationale' column: {"pvs": {"pvs1": [score, text]}, ...}.
A code is applied when score > 0. The strength suffix is added when it differs from the code's
default strength (BIAS scale: 1 supporting, 2 moderate, 3 strong, 4 very strong, 5 stand-alone).
usage: bias_add_codes.py in.tsv out.tsv
"""
import csv
import json
import sys

STRENGTH = {1: "Supporting", 2: "Moderate", 3: "Strong", 4: "VeryStrong", 5: "StandAlone"}
DEFAULT = {"pvs": 4, "ps": 3, "pm": 2, "pp": 1, "ba": 5, "bs": 3, "bp": 1}

csv.field_size_limit(sys.maxsize)
with open(sys.argv[1], newline="") as fi, open(sys.argv[2], "w", newline="") as fo:
    r = csv.DictReader(fi, delimiter="\t")
    cols = r.fieldnames[:]
    cols.insert(cols.index("acmgClassification") + 1, "acmgCodes")
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for row in r:
        codes = []
        try:
            rat = json.loads(row.get("rationale") or "{}")
        except ValueError:
            rat = {}
        for group, items in rat.items():
            for code, val in (items or {}).items():
                score = val[0] if isinstance(val, list) and val else 0
                if isinstance(score, (int, float)) and score > 0:
                    name = code.upper()
                    if DEFAULT.get(group) != score:
                        name += "_" + STRENGTH.get(int(score), str(score))
                    codes.append(name)
        row["acmgCodes"] = ",".join(codes)
        w.writerow(row)
