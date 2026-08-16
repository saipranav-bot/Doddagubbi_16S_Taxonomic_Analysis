#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from pathlib import Path

# ------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------

Entrez.email = "YOUR_REAL_EMAIL@example.com"

BASE = Path.home() / "Doddagubbi_Master"

OUT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "ncbi_check"
)

OUT.mkdir(parents=True, exist_ok=True)

ACCESSIONS = {
    "PZ740679.1": "KJU0224",
    "PZ740680.1": "KJU0225",
}

# ------------------------------------------------------------
# DOWNLOAD
# ------------------------------------------------------------

print("=" * 80)
print("DOWNLOAD MISSING OPENCITY KJU SEQUENCES")
print("=" * 80)

records = []

for accession, kju in ACCESSIONS.items():

    print("\n" + "-" * 80)
    print(f"KJU: {kju}")
    print(f"Accession: {accession}")
    print("-" * 80)

    # Download GenBank
    handle = Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    gb_text = handle.read()
    handle.close()

    gb_file = OUT / f"{kju}_{accession}.gb"
    gb_file.write_text(gb_text)

    # Parse
    record = SeqIO.read(str(gb_file), "genbank")

    print("Accession :", record.id)
    print("Length    :", len(record.seq))
    print("Description:")
    print(record.description)

    # Extract organism
    organism = "UNKNOWN"

    for feature in record.features:
        if feature.type == "source":
            organism = feature.qualifiers.get(
                "organism",
                ["UNKNOWN"]
            )[0]
            break

    print("Organism  :", organism)

    # Save FASTA with KJU ID
    fasta_file = OUT / f"{kju}_{accession}.fasta"

    record.id = kju
    record.name = kju
    record.description = (
        f"{kju} | {accession} | {organism}"
    )

    SeqIO.write(record, str(fasta_file), "fasta")

    print("GenBank   :", gb_file)
    print("FASTA     :", fasta_file)

    records.append(record)

# ------------------------------------------------------------
# COMBINED FASTA
# ------------------------------------------------------------

combined = OUT / "OPENCITY_MISSING_KJU0224_KJU0225.fasta"

SeqIO.write(records, str(combined), "fasta")

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)

print("Combined FASTA:")
print(combined)

print("\nSequences downloaded:", len(records))
