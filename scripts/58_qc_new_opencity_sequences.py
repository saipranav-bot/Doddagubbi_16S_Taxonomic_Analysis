#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
from collections import Counter

BASE = Path.home() / "Doddagubbi_Master"

INPUT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "ncbi_check"
    / "OPENCITY_MISSING_KJU0224_KJU0225.fasta"
)

OUT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "ncbi_check"
    / "NEW_SEQUENCE_QC_REPORT.txt"
)

records = list(SeqIO.parse(INPUT, "fasta"))

lines = []

lines.append("=" * 80)
lines.append("QC OF NEW OPENCITY SEQUENCES")
lines.append("=" * 80)

for record in records:

    seq = str(record.seq).upper()

    length = len(seq)
    counts = Counter(seq)

    valid = set("ACGTN")
    invalid = sorted(set(seq) - valid)

    n_count = counts.get("N", 0)
    gc_count = counts.get("G", 0) + counts.get("C", 0)

    gc_percent = (
        (gc_count / length) * 100
        if length > 0
        else 0
    )

    lines.append("")
    lines.append("-" * 80)
    lines.append(f"ID: {record.id}")
    lines.append("-" * 80)

    lines.append(f"Length: {length}")
    lines.append(f"A: {counts.get('A', 0)}")
    lines.append(f"C: {counts.get('C', 0)}")
    lines.append(f"G: {counts.get('G', 0)}")
    lines.append(f"T: {counts.get('T', 0)}")
    lines.append(f"N: {n_count}")
    lines.append(f"GC%: {gc_percent:.2f}")
    lines.append(f"Invalid characters: {invalid}")

    if invalid:
        status = "FAIL - invalid nucleotide characters"
    elif length < 500:
        status = "REVIEW - short sequence"
    elif n_count > 0:
        status = "REVIEW - ambiguous bases present"
    else:
        status = "PASS - basic sequence QC"

    lines.append(f"QC STATUS: {status}")

lines.append("")
lines.append("=" * 80)
lines.append("END OF QC")
lines.append("=" * 80)

OUT.write_text("\n".join(lines))

print("\n".join(lines))
print("\nReport saved to:")
print(OUT)
