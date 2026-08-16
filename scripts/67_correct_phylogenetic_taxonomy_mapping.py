#!/usr/bin/env python3

from pathlib import Path
from Bio import Phylo, SeqIO
import pandas as pd
import re

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

TREEDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
    / "phylogenetic_trees"
)

TAXDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

REFERENCE_FASTA = TAXDIR / "ALL_NCBI_TAXONOMIC_REFERENCES.fasta"

OUTFILE = TAXDIR / "CORRECTED_PHYLOGENETIC_PLACEMENT.csv"
SUMMARY = TAXDIR / "CORRECTED_PHYLOGENETIC_PLACEMENT_SUMMARY.txt"

targets = {
    "KJU0224": {
        "claimed_species": "Pseudomonas putida",
        "tree": TREEDIR / "KJU0224_FastTree.nwk",
    },
    "KJU0225": {
        "claimed_species": "Chryseobacterium indologenes",
        "tree": TREEDIR / "KJU0225_FastTree.nwk",
    },
}

# ============================================================
# TAXONOMY PARSER
# ============================================================

def clean_taxonomy(description):

    text = description.strip()

    # Remove common sequence-description suffixes
    text = re.sub(
        r"\s+16S ribosomal RNA.*$",
        "",
        text,
        flags=re.IGNORECASE
    )

    text = re.sub(
        r"\s+16S rRNA.*$",
        "",
        text,
        flags=re.IGNORECASE
    )

    # Remove "strain ..."
    text = re.sub(
        r"\s+strain\s+.*$",
        "",
        text,
        flags=re.IGNORECASE
    )

    text = text.strip()

    # Extract species-like name
    parts = text.split()

    if len(parts) >= 2:
        genus = parts[0]
        species = parts[0] + " " + parts[1]
    elif len(parts) == 1:
        genus = parts[0]
        species = parts[0]
    else:
        genus = ""
        species = ""

    return genus, species


# ============================================================
# LOAD NCBI TAXONOMIC INFORMATION
# ============================================================

print("=" * 90)
print("STEP 67 — CORRECTED PHYLOGENETIC TAXONOMY MAPPING")
print("=" * 90)

reference_taxonomy = {}

for record in SeqIO.parse(REFERENCE_FASTA, "fasta"):

    accession = record.id.split(".")[0]

    description = record.description

    genus, species = clean_taxonomy(description)

    reference_taxonomy[record.id] = {
        "accession": record.id,
        "genus": genus,
        "species": species,
        "description": description,
    }

    reference_taxonomy[accession] = reference_taxonomy[record.id]


print("\nReference taxonomy records loaded:",
      len(reference_taxonomy))

# ============================================================
# ANALYSIS
# ============================================================

results = []
summary_lines = []

summary_lines.append("=" * 90)
summary_lines.append("CORRECTED PHYLOGENETIC PLACEMENT SUMMARY")
summary_lines.append("=" * 90)

for kju, info in targets.items():

    claimed_species = info["claimed_species"]
    claimed_genus = claimed_species.split()[0]

    print("\n" + "=" * 90)
    print("TARGET:", kju)
    print("CLAIMED SPECIES:", claimed_species)
    print("=" * 90)

    tree = Phylo.read(str(info["tree"]), "newick")

    terminals = tree.get_terminals()

    target = None

    for terminal in terminals:
        if terminal.name == kju:
            target = terminal
            break

    if target is None:
        print("ERROR: target not found in tree")
        continue

    # --------------------------------------------------------
    # Calculate distances
    # --------------------------------------------------------

    neighbors = []

    for terminal in terminals:

        if terminal.name == kju:
            continue

        distance = tree.distance(target, terminal)

        accession = terminal.name
        tax = reference_taxonomy.get(accession)

        if tax is None:
            tax = reference_taxonomy.get(
                accession.split(".")[0],
                {
                    "accession": accession,
                    "genus": "UNKNOWN",
                    "species": "UNKNOWN",
                    "description": "",
                }
            )

        neighbors.append({
            "accession": accession,
            "distance": distance,
            "genus": tax["genus"],
            "species": tax["species"],
            "description": tax["description"],
        })

    neighbors = sorted(
        neighbors,
        key=lambda x: x["distance"]
    )

    top10 = neighbors[:10]

    print("\nTOP 10 PHYLOGENETIC NEIGHBORS")
    print("-" * 90)

    genus_matches = 0
    species_matches = 0

    for i, n in enumerate(top10, 1):

        genus_match = (
            n["genus"].lower() == claimed_genus.lower()
        )

        species_match = (
            n["species"].lower() == claimed_species.lower()
        )

        if genus_match:
            genus_matches += 1

        if species_match:
            species_matches += 1

        print(
            f"{i:2d}. "
            f"{n['accession']} | "
            f"{n['species']} | "
            f"distance={n['distance']:.6f} | "
            f"genus_match={genus_match} | "
            f"species_match={species_match}"
        )

        results.append({
            "query": kju,
            "claimed_species": claimed_species,
            "reference_accession": n["accession"],
            "reference_genus": n["genus"],
            "reference_species": n["species"],
            "distance": n["distance"],
            "genus_match": genus_match,
            "species_match": species_match,
            "reference_description": n["description"],
        })

    # --------------------------------------------------------
    # Best neighbor
    # --------------------------------------------------------

    best = top10[0]

    best_genus_match = (
        best["genus"].lower() == claimed_genus.lower()
    )

    best_species_match = (
        best["species"].lower() == claimed_species.lower()
    )

    print("\nBEST PHYLOGENETIC NEIGHBOR")
    print("-" * 90)
    print("Reference :", best["accession"])
    print("Organism  :", best["species"])
    print("Distance  :", f"{best['distance']:.6f}")
    print("Genus     :", best_genus_match)
    print("Species   :", best_species_match)

    # --------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------

    if species_matches >= 3:
        interpretation = (
            "STRONG PHYLOGENETIC SUPPORT FOR CLAIMED SPECIES"
        )
    elif genus_matches >= 3:
        interpretation = (
            "GENUS SUPPORTED; SPECIES-LEVEL SUPPORT LIMITED"
        )
    elif best_genus_match:
        interpretation = (
            "BEST NEIGHBOR SUPPORTS CLAIMED GENUS"
        )
    else:
        interpretation = (
            "PHYLOGENETIC PLACEMENT DOES NOT SUPPORT CLAIMED GENUS"
        )

    print("\nTOP-10 TAXONOMIC AGREEMENT")
    print("-" * 90)
    print("Genus matches   :", genus_matches, "/ 10")
    print("Species matches :", species_matches, "/ 10")

    print("\nINTERPRETATION")
    print("-" * 90)
    print(interpretation)

    summary_lines.append("")
    summary_lines.append("=" * 90)
    summary_lines.append(kju)
    summary_lines.append("=" * 90)
    summary_lines.append(
        f"Claimed species: {claimed_species}"
    )
    summary_lines.append(
        f"Best neighbor: {best['accession']}"
    )
    summary_lines.append(
        f"Best neighbor organism: {best['species']}"
    )
    summary_lines.append(
        f"Best distance: {best['distance']:.6f}"
    )
    summary_lines.append(
        f"Genus matches among top 10: {genus_matches}/10"
    )
    summary_lines.append(
        f"Species matches among top 10: {species_matches}/10"
    )
    summary_lines.append(
        f"Interpretation: {interpretation}"
    )

# ============================================================
# SAVE
# ============================================================

df = pd.DataFrame(results)

df.to_csv(OUTFILE, index=False)

SUMMARY.write_text(
    "\n".join(summary_lines) + "\n"
)

print("\n" + "=" * 90)
print("STEP 67 COMPLETE")
print("=" * 90)

print("\nResults:")
print(OUTFILE)

print("\nSummary:")
print(SUMMARY)

