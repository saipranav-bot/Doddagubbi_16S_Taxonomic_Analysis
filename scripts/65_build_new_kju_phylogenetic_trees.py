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

COMBINED = TAXDIR / "NEW_KJUS_PLUS_NCBI_REFERENCES.fasta"

TREEDIR = TAXDIR / "phylogenetic_trees"

TREEDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# CHECK TOOLS
# ============================================================

print("=" * 80)
print("STEP 65 — PHYLOGENETIC VALIDATION")
print("=" * 80)

for tool in ["mafft", "FastTree"]:

    path = shutil.which(tool)

    if path:
        print(f"{tool}: FOUND -> {path}")
    else:
        print(f"{tool}: NOT FOUND")
        raise SystemExit(
            f"\nERROR: {tool} is required."
        )

# ============================================================
# LOAD COMBINED DATASET
# ============================================================

records = list(
    SeqIO.parse(
        COMBINED,
        "fasta"
    )
)

print("\nTotal sequences:", len(records))

# ============================================================
# TARGETS
# ============================================================

targets = [
    "KJU0224",
    "KJU0225"
]

# ============================================================
# BUILD ONE TREE PER KJU
# ============================================================

for target in targets:

    print("\n" + "=" * 80)
    print("TARGET:", target)
    print("=" * 80)

    target_record = None

    for record in records:

        if record.id == target:
            target_record = record
            break

    if target_record is None:

        print(
            f"ERROR: {target} not found."
        )

        continue

    # --------------------------------------------------------
    # Determine relevant genus
    # --------------------------------------------------------

    if target == "KJU0224":

        genus = "Pseudomonas"

    elif target == "KJU0225":

        genus = "Chryseobacterium"

    # --------------------------------------------------------
    # Select references
    # --------------------------------------------------------

    selected = [
        target_record
    ]

    for record in records:

        if record.id == target:
            continue

        description = (
            record.description
            .lower()
        )

        if genus.lower() in description:

            selected.append(record)

    print(
        "Target sequence:",
        target_record.id
    )

    print(
        "Selected genus references:",
        len(selected) - 1
    )

    print(
        "Total sequences for tree:",
        len(selected)
    )

    # --------------------------------------------------------
    # Save FASTA
    # --------------------------------------------------------

    fasta_file = (
        TREEDIR
        / f"{target}_PHYLOGENETIC_DATASET.fasta"
    )

    with open(
        fasta_file,
        "w"
    ) as handle:

        for record in selected:

            SeqIO.write(
                record,
                handle,
                "fasta"
            )

    print(
        "\nDataset saved:",
        fasta_file
    )

    # --------------------------------------------------------
    # MAFFT
    # --------------------------------------------------------

    alignment_file = (
        TREEDIR
        / f"{target}_MAFFT.fasta"
    )

    print("\nRunning MAFFT...")

    with open(
        alignment_file,
        "w"
    ) as outfile:

        subprocess.run(
            [
                "mafft",
                "--auto",
                str(fasta_file)
            ],
            stdout=outfile,
            check=True
        )

    print(
        "Alignment saved:",
        alignment_file
    )

    # --------------------------------------------------------
    # FastTree
    # --------------------------------------------------------

    tree_file = (
        TREEDIR
        / f"{target}_FastTree.nwk"
    )

    print("\nRunning FastTree...")

    with open(
        tree_file,
        "w"
    ) as outfile:

        subprocess.run(
            [
                "FastTree",
                "-nt",
                str(alignment_file)
            ],
            stdout=outfile,
            check=True
        )

    print(
        "Tree saved:",
        tree_file
    )

# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("STEP 65 COMPLETE")
print("=" * 80)

print("\nGenerated files:")

for file in sorted(
    TREEDIR.glob("*")
):

    print(
        " -",
        file.name
    )

print("\nTrees are ready for inspection.")

