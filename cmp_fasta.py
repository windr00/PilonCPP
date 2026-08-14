#!/usr/bin/env python3
"""Compare two FASTA files ignoring line-wrapping and order."""

import sys

def parse_fasta(path):
    """Return {name: seq} mapping with sequence concatenated and uppercased."""
    records = {}
    with open(path) as f:
        name = None
        seqs = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(seqs).upper()
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
        if name is not None:
            records[name] = ''.join(seqs).upper()
    return records

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} file_a.fasta file_b.fasta", file=sys.stderr)
        sys.exit(2)

    a = parse_fasta(sys.argv[1])
    b = parse_fasta(sys.argv[2])

    ok = True

    # Contigs only in one file
    only_a = set(a) - set(b)
    only_b = set(b) - set(a)
    for name in sorted(only_a):
        print(f"MISSING  {name}: only in {sys.argv[1]}", file=sys.stderr)
        ok = False
    for name in sorted(only_b):
        print(f"MISSING  {name}: only in {sys.argv[2]}", file=sys.stderr)
        ok = False

    # Sequence content comparison
    for name in sorted(set(a) & set(b)):
        sa = a[name]
        sb = b[name]
        if sa == sb:
            continue
        ok = False
        minlen = min(len(sa), len(sb))
        diff_pos = next((i for i in range(minlen) if sa[i] != sb[i]), None)
        if diff_pos is not None:
            # Show surrounding context
            lo = max(diff_pos - 10, 0)
            hi = min(diff_pos + 11, minlen)
            print(f"DIFF     {name}: length_a={len(sa)} length_b={len(sb)}", file=sys.stderr)
            print(f"         a[{lo}:{hi}] = {repr(sa[lo:hi])}", file=sys.stderr)
            print(f"         b[{lo}:{hi}] = {repr(sb[lo:hi])}", file=sys.stderr)
        else:
            print(f"LENGTH   {name}: a={len(sa)} b={len(sb)}", file=sys.stderr)

    if ok:
        print("OK: files are identical (ignoring line-wrapping and order)", file=sys.stderr)

    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
