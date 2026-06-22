#!/usr/bin/env python3
"""
Lab 6 helper — summarize_dtl.py

Reads a GeneRax `<family>_speciesEventCounts.txt` file and prints a short,
readable summary: total D/T/L events, and which species-tree branches
(if any) show something other than plain speciation.

Usage:
    python3 scripts/summarize_dtl.py path/to/<family>_speciesEventCounts.txt
"""
import sys

def parse_event_counts(path):
    """
    GeneRax writes one line per species-tree node:
        species  speciations  duplications  losses  transfers
    The exact column separator/header can vary slightly between
    GeneRax versions, so this parser is deliberately tolerant: it skips
    any line that doesn't end in 4 numeric columns, and treats the first
    column as the species/node name.
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            name = parts[0]
            try:
                spec, dup, loss, transfer = (float(x) for x in parts[1:5])
            except ValueError:
                continue  # header row
            rows.append((name, spec, dup, loss, transfer))
    return rows


def main():
    if len(sys.argv) != 2:
        print(f"Usage: python3 {sys.argv[0]} <speciesEventCounts.txt>")
        sys.exit(1)

    path = sys.argv[1]
    rows = parse_event_counts(path)

    if not rows:
        print(f"Could not parse any event-count rows from {path}")
        print("(Check that the file path is correct and GeneRax finished successfully.)")
        sys.exit(1)

    total_spec = sum(r[1] for r in rows)
    total_dup = sum(r[2] for r in rows)
    total_loss = sum(r[3] for r in rows)
    total_transfer = sum(r[4] for r in rows)

    print(f"=== {path} ===")
    print(f"Total speciations : {total_spec:.0f}")
    print(f"Total duplications: {total_dup:.0f}")
    print(f"Total losses      : {total_loss:.0f}")
    print(f"Total transfers   : {total_transfer:.0f}")
    print()

    interesting = [r for r in rows if r[2] > 0 or r[3] > 0 or r[4] > 0]

    if not interesting:
        print("No branch shows duplication, loss, or transfer events.")
        print("This gene family is fully explained by speciation alone —")
        print("i.e. a 'clean', vertically-inherited gene family.")
    else:
        print("Branches with non-speciation events:")
        print(f"{'species/node':32s} {'dup':>5s} {'loss':>5s} {'transfer':>9s}")
        for name, spec, dup, loss, transfer in interesting:
            print(f"{name:32s} {dup:5.0f} {loss:5.0f} {transfer:9.0f}")


if __name__ == "__main__":
    main()
