#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from pathlib import Path
import time

# ============================================================
# SETTINGS
# ============================================================

Entrez.email = "saipranavkar@gmail.com"

BASE = Path.home() / "Doddagubbi_Master"

OUTDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# TARGETS
# ============================================================

targets = {
    "KJU0224": {
        "accession": "PZ740679.1",
        "claimed_species": "Pseudomonas putida",
        "genera": ["Pseudomonas"],
    },

    "KJU0225": {
        "accession": "PZ740680.1",
        "claimed_species": "Chryseobacterium indologenes",
        "genera": ["Chryseobacterium"],
    },
}

# ============================================================
# SEARCH FUNCTION
# ============================================================

def search_ncbi(term, retmax=20):

    print(f"\nNCBI SEARCH:")
    print(term)

    handle = Entrez.esearch(
        db="nuccore",
        term=term,
        retmax=retmax,
        sort="relevance"
    )

    result = Entrez.read(handle)
    handle.close()

    return result["IdList"]


# ============================================================
# DOWNLOAD FUNCTION
# ============================================================

def download_records(ids):

    if not ids:
        return []

    handle = Entrez.efetch(
        db="nuccore",
        id=",".join(ids),
        rettype="fasta",
        retmode="text"
    )

    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    return records


# ============================================================
# MAIN
# ============================================================

print("=" * 80)
print("STEP 61 — NCBI TAXONOMIC REFERENCE COLLECTION")
print("=" * 80)

all_reference_records = []

for kju, info in targets.items():

    print("\n" + "=" * 80)
    print(kju)
    print("=" * 80)

    species = info["claimed_species"]
    genus = info["genera"][0]

    # --------------------------------------------------------
    # Search claimed species
    # --------------------------------------------------------

    species_term = (
        f'"{species}"[Organism] AND '
        f'16S[Title] AND '
        f'("type strain"[Title] OR "reference"[Title] OR '
        f'"ribosomal RNA"[Title])'
    )

    species_ids = search_ncbi(species_term, retmax=15)

    species_records = download_records(species_ids)

    print("Species reference records:", len(species_records))

    # --------------------------------------------------------
    # Search genus
    # --------------------------------------------------------

    genus_term = (
        f'"{genus}"[Organism] AND '
        f'16S[Title] AND '
        f'("ribosomal RNA"[Title] OR "16S ribosomal RNA"[Title])'
    )

    genus_ids = search_ncbi(genus_term, retmax=20)

    genus_records = download_records(genus_ids)

    print("Genus reference records:", len(genus_records))

    # --------------------------------------------------------
    # Remove target sequence itself
    # --------------------------------------------------------

    target_acc = info["accession"].split(".")[0]

    filtered = []

    for record in species_records + genus_records:

        record_acc = record.id.split(".")[0]

        if record_acc == target_acc:
            continue

        filtered.append(record)

    # --------------------------------------------------------
    # Deduplicate
    # --------------------------------------------------------

    unique = {}

    for record in filtered:
        unique[record.id] = record

    filtered = list(unique.values())

    print("Unique references:", len(filtered))

    # --------------------------------------------------------
    # Save per-target FASTA
    # --------------------------------------------------------

    outfile = OUTDIR / f"{kju}_NCBI_REFERENCES.fasta"

    SeqIO.write(filtered, outfile, "fasta")

    print("Saved:", outfile)

    # --------------------------------------------------------
    # Add target sequence label later
    # --------------------------------------------------------

    for record in filtered:
        record.description = f"{kju}_REFERENCE | {record.description}"
        all_reference_records.append(record)

    time.sleep(1)


# ============================================================
# COMBINED REFERENCE FASTA
# ============================================================

combined = OUTDIR / "ALL_NCBI_TAXONOMIC_REFERENCES.fasta"

SeqIO.write(all_reference_records, combined, "fasta")

print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

print("Total reference sequences:", len(all_reference_records))
print("Combined FASTA:", combined)

print("\nDONE.")
