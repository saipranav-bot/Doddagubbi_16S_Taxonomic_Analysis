#!/usr/bin/env python3

from pathlib import Path
from Bio import Phylo
from Bio import SeqIO
import pandas as pd

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

OUTDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

RESULTS = []

targets = {
    "KJU0224": {
        "claimed_species": "Pseudomonas putida",
        "tree": TREEDIR / "KJU0224_FastTree.nwk",
        "alignment": TREEDIR / "KJU0224_MAFFT.fasta",
    },
    "KJU0225": {
        "claimed_species": "Chryseobacterium indologenes",
        "tree": TREEDIR / "KJU0225_FastTree.nwk",
        "alignment": TREEDIR / "KJU0225_MAFFT.fasta",
    },
}

# ============================================================
# FUNCTION TO EXTRACT TAXONOMY
# ============================================================

def clean_name(name):

    name = str(name)

    if "|" in name:
        parts = name.split("|")

        if len(parts) >= 3:
            return parts[2].replace("_", " ")

    return name.replace("_", " ")


# ============================================================
# ANALYZE EACH TREE
# ============================================================

print("=" * 90)
print("STEP 66 — PHYLOGENETIC PLACEMENT ANALYSIS")
print("=" * 90)

for target, info in targets.items():

    print("\n" + "=" * 90)
    print("TARGET:", target)
    print("CLAIMED SPECIES:", info["claimed_species"])
    print("=" * 90)

    tree_file = info["tree"]

    tree = Phylo.read(tree_file, "newick")

    terminals = tree.get_terminals()

    print("Tree terminals:", len(terminals))

    target_node = None

    for node in terminals:

        if node.name == target:
            target_node = node
            break

    if target_node is None:

        print("ERROR: Target not found in tree")
        continue

    # --------------------------------------------------------
    # FIND CLOSEST TERMINAL NEIGHBORS
    # --------------------------------------------------------

    distances = []

    for node in terminals:

        if node == target_node:
            continue

        try:
            distance = tree.distance(target_node, node)
        except Exception:
            continue

        distances.append(
            {
                "target": target,
                "reference": node.name,
                "taxonomy": clean_name(node.name),
                "distance": distance,
            }
        )

    distances = sorted(
        distances,
        key=lambda x: x["distance"]
    )

    # --------------------------------------------------------
    # TOP 10
    # --------------------------------------------------------

    print("\nTOP 10 PHYLOGENETIC NEIGHBORS")
    print("-" * 90)

    for i, item in enumerate(distances[:10], 1):

        print(
            f"{i:2d}. "
            f"{item['reference']} | "
            f"{item['taxonomy']} | "
            f"distance={item['distance']:.6f}"
        )

    # --------------------------------------------------------
    # TAXONOMIC AGREEMENT
    # --------------------------------------------------------

    claimed_genus = info["claimed_species"].split()[0]
    claimed_species = info["claimed_species"]

    genus_matches = 0
    species_matches = 0

    for item in distances[:10]:

        tax = item["taxonomy"]

        if tax.startswith(claimed_genus):
            genus_matches += 1

        if tax == claimed_species:
            species_matches += 1

    # --------------------------------------------------------
    # BEST NEIGHBOR
    # --------------------------------------------------------

    best = distances[0]

    best_taxonomy = best["taxonomy"]

    if best_taxonomy == claimed_species:

        placement = "STRONG SPECIES SUPPORT"

    elif best_taxonomy.startswith(claimed_genus):

        placement = "GENUS SUPPORT — SPECIES NOT FULLY RESOLVED"

    else:

        placement = "CONFLICT WITH CLAIMED TAXONOMY"

    print("\nBEST PHYLOGENETIC NEIGHBOR")
    print("-" * 90)

    print("Reference :", best["reference"])
    print("Taxonomy  :", best_taxonomy)
    print("Distance  :", f"{best['distance']:.6f}")

    print("\nTOP-10 TAXONOMIC AGREEMENT")
    print("-" * 90)

    print("Genus matches   :", genus_matches, "/ 10")
    print("Species matches :", species_matches, "/ 10")

    print("\nPLACEMENT INTERPRETATION")
    print("-" * 90)
    print(placement)

    # --------------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------------

    for rank, item in enumerate(distances[:10], 1):

        RESULTS.append(
            {
                "query": target,
                "claimed_species": info["claimed_species"],
                "rank": rank,
                "reference": item["reference"],
                "taxonomy": item["taxonomy"],
                "phylogenetic_distance": item["distance"],
                "genus_match": item["taxonomy"].startswith(
                    claimed_genus
                ),
                "species_match": item["taxonomy"] == claimed_species,
                "placement_interpretation": placement,
            }
        )

# ============================================================
# SAVE CSV
# ============================================================

df = pd.DataFrame(RESULTS)

outfile = OUTDIR / "PHYLOGENETIC_PLACEMENT_RESULTS.csv"

df.to_csv(outfile, index=False)

# ============================================================
# SUMMARY
# ============================================================

summary_file = OUTDIR / "PHYLOGENETIC_PLACEMENT_SUMMARY.txt"

with open(summary_file, "w") as f:

    f.write("=" * 90 + "\n")
    f.write("PHYLOGENETIC PLACEMENT VALIDATION SUMMARY\n")
    f.write("=" * 90 + "\n\n")

    for target in targets:

        sub = df[df["query"] == target]

        if len(sub) == 0:
            continue

        first = sub.iloc[0]

        f.write(f"{target}\n")
        f.write("-" * 90 + "\n")
        f.write(
            f"Claimed species: "
            f"{first['claimed_species']}\n"
        )

        f.write(
            f"Closest reference: "
            f"{first['reference']}\n"
        )

        f.write(
            f"Closest taxonomy: "
            f"{first['taxonomy']}\n"
        )

        f.write(
            f"Phylogenetic distance: "
            f"{first['phylogenetic_distance']:.6f}\n"
        )

        f.write(
            f"Placement interpretation: "
            f"{first['placement_interpretation']}\n\n"
        )

print("\n" + "=" * 90)
print("STEP 66 COMPLETE")
print("=" * 90)

print("\nResults:")
print(outfile)

print("\nSummary:")
print(summary_file)

