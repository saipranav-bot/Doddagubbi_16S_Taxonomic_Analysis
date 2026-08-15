from pathlib import Path
from Bio import Phylo
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

TREE = META / "FINAL_109_PHYLOGENY_V1" / "KJU109_16S.treefile"
QC = META / "FINAL_109_SEQUENCE_QC_V2" / "FINAL_109_SEQUENCE_QC_V2.csv"

OUT = META / "FINAL_109_PHYLOGENY_V1" / \
      "SPECIES_RESOLUTION_REVIEW_V1.csv"

SUMMARY = META / "FINAL_109_PHYLOGENY_V1" / \
          "SPECIES_RESOLUTION_REVIEW_SUMMARY_V1.txt"

tree = Phylo.read(TREE, "newick")
qc = pd.read_csv(QC, dtype=str).fillna("")

def species_name(x):
    p = str(x).strip().split()
    return " ".join(p[:2]) if len(p) >= 2 else str(x).strip()

qc["species"] = qc["organism"].map(species_name)
qc["genus"] = (
    qc["organism"]
    .astype(str)
    .str.strip()
    .str.split()
    .str[0]
)

# ------------------------------------------------------------
# Tree map
# ------------------------------------------------------------

tips = {}

for tip in tree.get_terminals():

    parts = str(tip.name).split("|")

    if len(parts) < 3:
        continue

    kju = parts[0]
    accession = parts[1]
    organism = parts[2].replace("_", " ")

    tips[kju] = {
        "tip": tip,
        "accession": accession,
        "organism": organism,
        "species": species_name(organism),
    }

rows = []

# ------------------------------------------------------------
# Analyze every species represented >=2 times
# ------------------------------------------------------------

species_counts = qc["species"].value_counts()

for species, count in species_counts.items():

    if count < 2:
        continue

    members = qc[
        qc["species"] == species
    ]["kju_id"].tolist()

    tree_members = [
        tips[k]["tip"]
        for k in members
        if k in tips
    ]

    if len(tree_members) < 2:
        continue

    # --------------------------------------------------------
    # MRCA of all same-species members
    # --------------------------------------------------------

    try:
        mrca = tree.common_ancestor(tree_members)
        terminals = mrca.get_terminals()
    except Exception:
        continue

    # Species composition inside MRCA
    terminal_species = []

    for terminal in terminals:

        p = str(terminal.name).split("|")

        if len(p) < 3:
            continue

        org = p[2].replace("_", " ")
        terminal_species.append(
            species_name(org)
        )

    species_purity = (
        sum(
            x == species
            for x in terminal_species
        )
        / len(terminal_species)
        if terminal_species
        else 0
    )

    exclusive = all(
        x == species
        for x in terminal_species
    )

    # --------------------------------------------------------
    # Pairwise within-species distances
    # --------------------------------------------------------

    within_distances = []

    for i in range(len(tree_members)):

        for j in range(i + 1, len(tree_members)):

            d = tree.distance(
                tree_members[i],
                tree_members[j]
            )

            within_distances.append(d)

    min_within = min(within_distances)
    mean_within = sum(within_distances) / len(within_distances)
    max_within = max(within_distances)

    # --------------------------------------------------------
    # Closest different-species neighbour for the group
    # --------------------------------------------------------

    different_distances = []

    for kju in members:

        if kju not in tips:
            continue

        tip = tips[kju]["tip"]

        for other_kju, other in tips.items():

            if other_kju == kju:
                continue

            if other["species"] == species:
                continue

            d = tree.distance(
                tip,
                other["tip"]
            )

            different_distances.append(
                (
                    d,
                    other_kju,
                    other["organism"]
                )
            )

    different_distances.sort(
        key=lambda x: x[0]
    )

    nearest_different = (
        different_distances[0]
        if different_distances
        else None
    )

    # --------------------------------------------------------
    # Classification
    # --------------------------------------------------------

    if exclusive:

        classification = "SPECIES_MONOPHYLY_OBSERVED"

    elif species_purity >= 0.90:

        classification = "MOSTLY_SPECIES_COHERENT"

    else:

        classification = "SPECIES_MIXING_OBSERVED"

    rows.append({
        "species": species,
        "kju_count": len(members),
        "kju_ids": ";".join(members),
        "mrca_terminal_count": len(terminals),
        "species_purity_percent": round(
            species_purity * 100,
            2
        ),
        "min_within_species_distance": round(
            min_within,
            6
        ),
        "mean_within_species_distance": round(
            mean_within,
            6
        ),
        "max_within_species_distance": round(
            max_within,
            6
        ),
        "nearest_different_species_kju":
            nearest_different[1]
            if nearest_different else "",
        "nearest_different_species":
            nearest_different[2]
            if nearest_different else "",
        "nearest_different_species_distance":
            round(nearest_different[0], 6)
            if nearest_different else "",
        "classification":
            classification
    })

result = pd.DataFrame(rows)

if not result.empty:
    result = result.sort_values(
        ["classification", "species"]
    )

result.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = [
    "DODDAGGUBBI - SPECIES RESOLUTION REVIEW",
    "=" * 110,
    "",
    f"Species groups represented by >=2 KJUs: {len(result)}",
    "",
    "CLASSIFICATION",
    "-" * 110,
]

if not result.empty:

    for key, value in result[
        "classification"
    ].value_counts().items():

        summary.append(
            f"{key}: {value}"
        )

summary += [
    "",
    "INTERPRETATION",
    "-" * 110,
    "SPECIES_MONOPHYLY_OBSERVED:",
    "all observed members of the species occur in an exclusive clade.",
    "",
    "MOSTLY_SPECIES_COHERENT:",
    "the species occupies a highly pure clade, but the MRCA also",
    "contains one or more other taxa.",
    "",
    "SPECIES_MIXING_OBSERVED:",
    "the species is substantially intermixed with other species.",
    "",
    "IMPORTANT:",
    "16S species-level non-monophyly does not automatically mean",
    "the metadata annotation is wrong. Closely related bacteria",
    "can share highly similar 16S sequences.",
]

summary.append("")
summary.append("MULTI-KJU SPECIES")
summary.append("-" * 110)

if result.empty:
    summary.append("NONE")
else:
    for _, r in result.iterrows():

        summary.append(
            f"{r['species']} | "
            f"KJUs={r['kju_count']} | "
            f"purity={r['species_purity_percent']}% | "
            f"within_mean={r['mean_within_species_distance']} | "
            f"nearest_other="
            f"{r['nearest_different_species']} "
            f"({r['nearest_different_species_kju']}) | "
            f"other_distance="
            f"{r['nearest_different_species_distance']} | "
            f"{r['classification']}"
        )

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("=" * 110)
print(" SPECIES RESOLUTION REVIEW COMPLETE")
print("=" * 110)

print("")
print(result.to_string(index=False))

print("")
print("Output:")
print(OUT)

print("")
print("Summary:")
print(SUMMARY)

print("")
print("NO SOURCE DATA WAS MODIFIED.")
