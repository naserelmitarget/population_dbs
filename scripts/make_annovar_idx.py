#!/usr/bin/env python3
"""Build an ANNOVAR filter-database .idx (1 kb bins) for a sorted generic table.

Format (same as ANNOVAR-distributed .idx files):
    #BIN<TAB>1000<TAB><file size>
    <chr><TAB><bin start><TAB><byte start><TAB><byte end>
The table must be sorted by chr then start so every bin is one contiguous byte range.

usage:
  make_annovar_idx.py table.txt                 index one table -> table.txt.idx
  make_annovar_idx.py --part part.txt           relative index of one headerless chunk -> part.txt.pidx
  make_annovar_idx.py --merge out.txt p1 p2 ... out.txt is the concatenation of header + p1 + p2...;
                                                writes out.txt.idx from the .pidx files, no rescan
"""
import os
import sys

BIN = 1000


def scan(path, base=0):
    """Yield (chr, bin, start, end) byte ranges of `path`, offsets shifted by `base`."""
    cur, cur_start, seen, offset = None, 0, set(), 0
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b"#"):
                offset += len(line)
                continue
            f = line.split(b"\t", 2)
            key = (f[0].decode(), int(f[1]) // BIN * BIN)
            if key != cur:
                if cur is not None:
                    yield cur[0], cur[1], base + cur_start, base + offset
                if key in seen:
                    sys.exit(f"{path} not sorted: bin {key} reappears at byte {offset}")
                seen.add(key)
                cur, cur_start = key, offset
            offset += len(line)
    if cur is not None:
        yield cur[0], cur[1], base + cur_start, base + offset


def main(argv):
    if argv[0] == "--part":
        with open(argv[1] + ".pidx", "w") as out:
            for r in scan(argv[1]):
                out.write("\t".join(map(str, r)) + "\n")
    elif argv[0] == "--merge":
        table, parts = argv[1], argv[2:]
        with open(table, "rb") as fh:
            first = fh.readline()
            base = len(first) if first.startswith(b"#") else 0
        size = os.path.getsize(table)
        with open(table + ".idx", "w") as out:
            out.write(f"#BIN\t{BIN}\t{size}\n")
            for p in parts:
                for line in open(p + ".pidx"):
                    c, b, s, e = line.rstrip("\n").split("\t")
                    out.write(f"{c}\t{b}\t{int(s) + base}\t{int(e) + base}\n")
                base += os.path.getsize(p)
        if base != size:
            sys.exit(f"size mismatch: parts+header={base} table={size}")
    else:
        size = os.path.getsize(argv[0])
        with open(argv[0] + ".idx", "w") as out:
            out.write(f"#BIN\t{BIN}\t{size}\n")
            for r in scan(argv[0]):
                out.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
