#!/usr/bin/env python3

from pathlib import Path
from Bio import Entrez, SeqIO, Phylo
import pandas as pd
import time
import re

# ============================================================
# SETTINGS
# ============================================================

Entrez.email = "saipranavkar@gmail.com"

BASE = Path.home() / "Doddagubbi_Master"

TREEDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
    / "phylogenetic_trees"
)

OUTDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

OUTFILE = OUTDIR / "NCBI_BACKED_PHYLOGENETIC_PLACEMENT.csv"

SUMMARY = OUTDIR / "NCBI_BACKED_PHYLOGENETIC_PLACEMENT_SUMMARY.txt"

# ============================================================
# TARGETS
# ============================================================

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
# HELPER
# ============================================================

def accession_from_name(name):

    """
    Convert tree labels such as:

        PZ794310.1
        PZ794310.1 KJU0224_REFERENCE

    into:

        PZ794310.1
    """

    if not name:
        return None

    name = name.strip()

    return name.split()[0]


# ============================================================
# NCBI TAXONOMY FETCH
# ============================================================

def fetch_taxonomy(accession):

    print(f"  NCBI lookup: {accession}")

    try:

        handle = Entrez.efetch(
            db="nuccore",
            id=accession,
            rettype="gb",
            retmode="text"
        )

        record = SeqIO.read(handle, "genbank")

        handle.close()

        organism = record.annotations.get(
            "organism",
            ""
        )

        description = record.description

        # Extract taxonomy lineage if available
        taxonomy = record.annotations.get(
            "taxonomy",
            []
        )

        genus = ""
        species = ""

        if organism:

            parts = organism.split()

            if len(parts) >= 1:
                genus = parts[0]

            if len(parts) >= 2:
                species = parts[0] + " " + parts[1]

        return {
            "accession": record.id,
            "organism": organism,
            "genus": genus,
            "species": species,
            "description": description,
            "taxonomy": "; ".join(taxonomy),
        }

    except Exception as e:

        print(
            f"  ERROR retrieving {accession}: {e}"
        )

        return {
            "accession": accession,
            "organism": "UNKNOWN",
            "genus": "UNKNOWN",
            "species": "UNKNOWN",
            "description": "",
            "taxonomy": "",
        }


# ============================================================
# COLLECT TREE ACCESSIONS
# ============================================================

print("=" * 90)
print("STEP 68 — NCBI-BACKED PHYLOGENETIC TAXONOMY VALIDATION")
print("=" * 90)

all_taxonomy = {}

for kju, info in targets.items():

    tree = Phylo.read(
        str(info["tree"]),
        "newick"
    )

    for terminal in tree.get_terminals():

        if terminal.name == kju:
            continue

        accession = accession_from_name(
            terminal.name
        )

        if accession:
            all_taxonomy[accession] = None


print(
    "\nUnique reference accessions requiring NCBI validation:",
    len(all_taxonomy)
)

# ============================================================
# FETCH NCBI TAXONOMY
# ============================================================

print("\n" + "=" * 90)
print("FETCHING NCBI TAXONOMY")
print("=" * 90)

for accession in sorted(all_taxonomy):

    all_taxonomy[accession] = fetch_taxonomy(
        accession
    )

    time.sleep(0.35)

# ============================================================
# ANALYZE EACH TREE
# ============================================================

results = []

summary_lines = []

summary_lines.append(
    "=" * 90
)

summary_lines.append(
    "NCBI-BACKED PHYLOGENETIC PLACEMENT SUMMARY"
)

summary_lines.append(
    "=" * 90
)

for kju, info in targets.items():

    claimed_species = info["claimed_species"]

    claimed_genus = claimed_species.split()[0]

    print("\n" + "=" * 90)
    print("TARGET:", kju)
    print("CLAIMED SPECIES:", claimed_species)
    print("=" * 90)

    tree = Phylo.read(
        str(info["tree"]),
        "newick"
    )

    target = None

    for terminal in tree.get_terminals():

        if terminal.name == kju:

            target = terminal
            break

    if target is None:

        print("ERROR: target not found")

        continue

    neighbors = []

    for terminal in tree.get_terminals():

        if terminal.name == kju:
            continue

        accession = accession_from_name(
            terminal.name
        )

        tax = all_taxonomy.get(
            accession,
            {
                "organism": "UNKNOWN",
                "genus": "UNKNOWN",
                "species": "UNKNOWN",
                "description": "",
                "taxonomy": "",
            }
        )

        distance = tree.distance(
            target,
            terminal
        )

        genus_match = (
            tax["genus"].lower()
            == claimed_genus.lower()
        )

        species_match = (
            tax["species"].lower()
            == claimed_species.lower()
        )

        neighbors.append({

            "query": kju,

            "claimed_species":
                claimed_species,

            "reference_accession":
                accession,

            "reference_organism":
                tax["organism"],

            "reference_genus":
                tax["genus"],

            "reference_species":
                tax["species"],

            "distance":
                distance,

            "genus_match":
                genus_match,

            "species_match":
                species_match,

            "description":
                tax["description"],

            "taxonomy":
                tax["taxonomy"],
        })

    neighbors.sort(
        key=lambda x: x["distance"]
    )

    top10 = neighbors[:10]

    print("\nTOP 10 PHYLOGENETIC NEIGHBORS")
    print("-" * 90)

    genus_matches = 0
    species_matches = 0

    for i, n in enumerate(top10, 1):

        if n["genus_match"]:
            genus_matches += 1

        if n["species_match"]:
            species_matches += 1

        print(
            f"{i:2d}. "
            f"{n['reference_accession']} | "
            f"{n['reference_organism']} | "
            f"distance={n['distance']:.6f} | "
            f"genus={n['genus_match']} | "
            f"species={n['species_match']}"
        )

        results.append(n)

    best = top10[0]

    # --------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------

    if species_matches >= 3:

        interpretation = (
            "STRONG SPECIES-LEVEL PHYLOGENETIC SUPPORT"
        )

    elif species_matches >= 1:

        interpretation = (
            "SPECIES-LEVEL PHYLOGENETIC SUPPORT PRESENT"
        )

    elif genus_matches >= 3:

        interpretation = (
            "GENUS-LEVEL PHYLOGENETIC SUPPORT"
        )

    elif best["genus_match"]:

        interpretation = (
            "BEST PHYLOGENETIC NEIGHBOR SUPPORTS CLAIMED GENUS"
        )

    else:

        interpretation = (
            "NO SUPPORT FROM TOP PHYLOGENETIC NEIGHBORS"
        )

    print("\n" + "-" * 90)
    print("BEST PHYLOGENETIC NEIGHBOR")
    print("-" * 90)

    print(
        "Accession :",
        best["reference_accession"]
    )

    print(
        "Organism  :",
        best["reference_organism"]
    )

    print(
        "Distance  :",
        f"{best['distance']:.6f}"
    )

    print(
        "Genus     :",
        best["genus_match"]
    )

    print(
        "Species   :",
        best["species_match"]
    )

    print("\nTOP-10 TAXONOMIC AGREEMENT")
    print("-" * 90)

    print(
        "Genus matches   :",
        genus_matches,
        "/ 10"
    )

    print(
        "Species matches :",
        species_matches,
        "/ 10"
    )

    print("\nINTERPRETATION")
    print("-" * 90)

    print(interpretation)

    summary_lines.append("")
    summary_lines.append(
        "=" * 90
    )

    summary_lines.append(kju)

    summary_lines.append(
        "=" * 90
    )

    summary_lines.append(
        f"Claimed species: {claimed_species}"
    )

    summary_lines.append(
        f"Best neighbor accession: "
        f"{best['reference_accession']}"
    )

    summary_lines.append(
        f"Best neighbor organism: "
        f"{best['reference_organism']}"
    )

    summary_lines.append(
        f"Best phylogenetic distance: "
        f"{best['distance']:.6f}"
    )

    summary_lines.append(
        f"Genus matches: "
        f"{genus_matches}/10"
    )

    summary_lines.append(
        f"Species matches: "
        f"{species_matches}/10"
    )

    summary_lines.append(
        f"Interpretation: "
        f"{interpretation}"
    )

# ============================================================
# SAVE
# ============================================================

df = pd.DataFrame(results)

df.to_csv(
    OUTFILE,
    index=False
)

SUMMARY.write_text(
    "\n".join(summary_lines)
    + "\n"
)

print("\n" + "=" * 90)
print("STEP 68 COMPLETE")
print("=" * 90)

print("\nResults:")
print(OUTFILE)

print("\nSummary:")
print(SUMMARY)

