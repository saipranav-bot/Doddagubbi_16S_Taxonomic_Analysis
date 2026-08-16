#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
import subprocess
import shutil

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

TAXDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

NEW_FASTA = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "ncbi_check"
    / "OPENCITY_MISSING_KJU0224_KJU0225.fasta"
)

COMBINED = TAXDIR / "NEW_KJUS_PLUS_NCBI_REFERENCES.fasta"

ALIGNMENT = TAXDIR / "NEW_KJUS_PLUS_NCBI_REFERENCES_MAFFT.fasta"

REPORT = TAXDIR / "ALIGNMENT_REPORT.txt"

REFERENCE = TAXDIR / "ALL_NCBI_TAXONOMIC_REFERENCES.fasta"

# ============================================================
# CHECK MAFFT
# ============================================================

print("=" * 80)
print("STEP 62 — ALIGN NEW KJU SEQUENCES WITH NCBI TAXONOMIC REFERENCES")
print("=" * 80)

if shutil.which("mafft") is None:
    raise SystemExit(
        "\nERROR: MAFFT is not installed or not in PATH.\n"
        "Install it with:\n"
        "conda install -c bioconda mafft\n"
    )

# ============================================================
# LOAD SEQUENCES
# ============================================================

new_records = list(SeqIO.parse(NEW_FASTA, "fasta"))
ref_records = list(SeqIO.parse(REFERENCE, "fasta"))

print("\nNew KJU sequences:", len(new_records))
print("NCBI references  :", len(ref_records))

# ============================================================
# CREATE COMBINED FASTA
# ============================================================

# Make sure the KJU identifiers are clean and unique.
records = []

for record in new_records:

    record.id = record.id.split()[0]
    record.name = record.id
    record.description = ""

    records.append(record)

for record in ref_records:

    record.id = record.id.split()[0]
    record.name = record.id

    records.append(record)

SeqIO.write(records, COMBINED, "fasta")

print("\nCombined FASTA:")
print(COMBINED)

print("\nTotal sequences:", len(records))

# ============================================================
# RUN MAFFT
# ============================================================

print("\n" + "=" * 80)
print("RUNNING MAFFT")
print("=" * 80)

cmd = [
    "mafft",
    "--auto",
    str(COMBINED)
]

print("\nCommand:")
print(" ".join(cmd))

with open(ALIGNMENT, "w") as outfile:

    result = subprocess.run(
        cmd,
        stdout=outfile,
        stderr=subprocess.PIPE,
        text=True
    )

if result.returncode != 0:

    print("\nMAFFT ERROR:")
    print(result.stderr)

    raise SystemExit("MAFFT failed.")

# ============================================================
# VERIFY ALIGNMENT
# ============================================================

aligned = list(SeqIO.parse(ALIGNMENT, "fasta"))

if not aligned:
    raise SystemExit("ERROR: Alignment file is empty.")

alignment_lengths = sorted(
    set(len(record.seq) for record in aligned)
)

print("\nAlignment completed successfully.")

print("Aligned sequences:", len(aligned))
print("Alignment length :", len(aligned[0].seq))

if len(alignment_lengths) == 1:
    print("All sequences have the same aligned length.")
else:
    print("WARNING: Different alignment lengths detected:")
    print(alignment_lengths)

# ============================================================
# KJU CHECK
# ============================================================

kju_ids = {"KJU0224", "KJU0225"}

found_kjus = []

for record in aligned:

    if record.id in kju_ids:
        found_kjus.append(record.id)

print("\nNew KJU sequences found in alignment:")

for kju in sorted(found_kjus):
    print(" -", kju)

missing = kju_ids - set(found_kjus)

if missing:
    print("\nWARNING — Missing KJU sequences:")
    for kju in sorted(missing):
        print(" -", kju)

# ============================================================
# REPORT
# ============================================================

report = []

report.append("=" * 80)
report.append("STEP 62 — ALIGNMENT REPORT")
report.append("=" * 80)
report.append("")
report.append(f"New KJU sequences: {len(new_records)}")
report.append(f"NCBI references: {len(ref_records)}")
report.append(f"Combined sequences: {len(records)}")
report.append(f"Aligned sequences: {len(aligned)}")
report.append(f"Alignment length: {len(aligned[0].seq)}")
report.append("")
report.append("KJU sequences included:")
report.append("KJU0224")
report.append("KJU0225")
report.append("")
report.append("MAFFT status: SUCCESS")
report.append("")
report.append(f"Combined FASTA: {COMBINED}")
report.append(f"Alignment FASTA: {ALIGNMENT}")

REPORT.write_text("\n".join(report))

# ============================================================
# DONE
# ============================================================

print("\n" + "=" * 80)
print("STEP 62 COMPLETE")
print("=" * 80)

print("\nGenerated:")
print(COMBINED)
print(ALIGNMENT)
print(REPORT)

