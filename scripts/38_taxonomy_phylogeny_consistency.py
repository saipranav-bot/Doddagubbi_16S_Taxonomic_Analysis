from pathlib import Path
from Bio import Phylo
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

TREE = META / "FINAL_109_PHYLOGENY_V1" / "KJU109_16S.treefile"
QC = META / "FINAL_109_SEQUENCE_QC_V2" / "FINAL_109_SEQUENCE_QC_V2.csv"

OUT = META / "FINAL_109_PHYLOGENY_V1" / \
      "TAXONOMY_PHYLOGENY_CONSISTENCY_V1.csv"

SUMMARY = META / "FINAL_109_PHYLOGENY_V1" / \
          "TAXONOMY_PHYLOGENY_CONSISTENCY_SUMMARY_V1.txt"

print("=" * 100)
print(" DODDAGGUBBI - TAXONOMY ↔ PHYLOGENY CONSISTENCY")
print("=" * 100)

tree = Phylo.read(TREE, "newick")
qc = pd.read_csv(QC, dtype=str).fillna("")

# ------------------------------------------------------------
# Parse tree tips
# ------------------------------------------------------------

tip_info = []

for tip in tree.get_terminals():

    name = tip.name

    parts = name.split("|")

    kju = parts[0] if len(parts) > 0 else ""
    accession = parts[1] if len(parts) > 1 else ""

    organism = (
        parts[2].replace("_", " ")
        if len(parts) > 2
        else ""
    )

    tip_info.append({
        "kju_id": kju,
        "accession": accession,
        "organism": organism,
        "tip_name": name
    })

tips = pd.DataFrame(tip_info)

# ------------------------------------------------------------
# Basic validation
# ------------------------------------------------------------

qc_sub = qc[
    [
        "kju_id",
        "accession",
        "organism"
    ]
].copy()

qc_sub["organism"] = (
    qc_sub["organism"]
    .str.strip()
)

merged = tips.merge(
    qc_sub,
    on=["kju_id", "accession"],
    how="left",
    suffixes=("_tree", "_metadata")
)

# ------------------------------------------------------------
# Genus / species parsing
# ------------------------------------------------------------

def genus(x):
    parts = str(x).split()
    return parts[0] if parts else ""

def species(x):
    parts = str(x).split()
    return parts[1] if len(parts) > 1 else ""

merged["tree_genus"] = merged["organism_tree"].map(genus)
merged["metadata_genus"] = merged["organism_metadata"].map(genus)

merged["tree_species"] = merged["organism_tree"].map(species)
merged["metadata_species"] = merged["organism_metadata"].map(species)

merged["annotation_self_match"] = (
    merged["organism_tree"].str.lower()
    ==
    merged["organism_metadata"].str.lower()
)

# ------------------------------------------------------------
# Tree-nearest-neighbor analysis
# ------------------------------------------------------------

rows = []

for idx, row in merged.iterrows():

    kju = row["kju_id"]
    tip_name = row["tip_name"]

    tip = next(
        t for t in tree.get_terminals()
        if t.name == tip_name
    )

    distances = []

    for other in tree.get_terminals():

        if other.name == tip_name:
            continue

        try:
            d = tree.distance(
                tip,
                other
            )
        except Exception:
            continue

        other_row = merged[
            merged["tip_name"] == other.name
        ]

        if other_row.empty:
            continue

        orow = other_row.iloc[0]

        distances.append({
            "other_kju": orow["kju_id"],
            "other_organism": orow["organism_metadata"],
            "other_genus": orow["metadata_genus"],
            "distance": d
        })

    distances.sort(
        key=lambda x: x["distance"]
    )

    nearest = distances[:5]

    same_genus_nearest = None
    first_different_genus = None

    for d in distances:

        if (
            d["other_genus"]
            ==
            row["metadata_genus"]
        ):
            same_genus_nearest = d
            break

    for d in distances:

        if (
            d["other_genus"]
            !=
            row["metadata_genus"]
        ):
            first_different_genus = d
            break

    nearest_kjus = ";".join(
        x["other_kju"]
        for x in nearest
    )

    nearest_organisms = ";".join(
        x["other_organism"]
        for x in nearest
    )

    nearest_same_genus = (
        same_genus_nearest["other_kju"]
        if same_genus_nearest
        else ""
    )

    nearest_same_genus_distance = (
        same_genus_nearest["distance"]
        if same_genus_nearest
        else ""
    )

    nearest_different_genus = (
        first_different_genus["other_kju"]
        if first_different_genus
        else ""
    )

    nearest_different_genus_distance = (
        first_different_genus["distance"]
        if first_different_genus
        else ""
    )

    # Flag possible concern if nearest neighbor belongs
    # to a different genus.
    nearest_neighbor_genus = (
        nearest[0]["other_genus"]
        if nearest
        else ""
    )

    possible_phylogenetic_conflict = (
        nearest_neighbor_genus != ""
        and
        nearest_neighbor_genus
        !=
        row["metadata_genus"]
    )

    rows.append({
        "kju_id": kju,
        "accession": row["accession"],
        "annotated_organism":
            row["organism_metadata"],
        "annotated_genus":
            row["metadata_genus"],
        "nearest_neighbor_kju":
            nearest[0]["other_kju"]
            if nearest else "",
        "nearest_neighbor_organism":
            nearest[0]["other_organism"]
            if nearest else "",
        "nearest_neighbor_genus":
            nearest_neighbor_genus,
        "nearest_neighbor_distance":
            nearest[0]["distance"]
            if nearest else "",
        "nearest_5_kjus":
            nearest_kjus,
        "nearest_5_organisms":
            nearest_organisms,
        "nearest_same_genus_kju":
            nearest_same_genus,
        "nearest_same_genus_distance":
            nearest_same_genus_distance,
        "nearest_different_genus_kju":
            nearest_different_genus,
        "nearest_different_genus_distance":
            nearest_different_genus_distance,
        "possible_phylogenetic_conflict":
            possible_phylogenetic_conflict
    })

result = pd.DataFrame(rows)

result.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

conflicts = result[
    result["possible_phylogenetic_conflict"]
]

same_genus = (
    result["nearest_neighbor_genus"]
    ==
    result["annotated_genus"]
).sum()

summary = []

summary.append(
    "DODDAGGUBBI - TAXONOMY ↔ PHYLOGENY CONSISTENCY"
)

summary.append("=" * 100)
summary.append("")

summary.append(
    f"Tree tips analyzed: {len(result)}"
)

summary.append(
    f"Nearest-neighbor same-genus: {same_genus}"
)

summary.append(
    f"Nearest-neighbor different-genus: "
    f"{len(result) - same_genus}"
)

summary.append(
    f"Potential phylogenetic conflicts: "
    f"{len(conflicts)}"
)

summary.append("")
summary.append(
    "IMPORTANT:"
)
summary.append(
    "A different-genus nearest neighbor is a screening signal only."
)
summary.append(
    "It does NOT by itself establish misidentification."
)
summary.append(
    "Taxonomic conclusions require clade-level inspection and,"
)
summary.append(
    "where necessary, additional markers or genome-scale evidence."
)

summary.append("")
summary.append("POTENTIAL SCREENING CANDIDATES")
summary.append("-" * 100)

if conflicts.empty:
    summary.append("NONE")
else:

    for _, r in conflicts.sort_values(
        "nearest_neighbor_distance"
    ).iterrows():

        summary.append(
            f"{r['kju_id']} | "
            f"{r['annotated_organism']} | "
            f"nearest={r['nearest_neighbor_kju']} "
            f"({r['nearest_neighbor_organism']}) | "
            f"distance={r['nearest_neighbor_distance']}"
        )

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("")
print("=" * 100)
print("CONSISTENCY ANALYSIS COMPLETE")
print("=" * 100)

print("Tree tips:", len(result))
print("Nearest neighbor same genus:", same_genus)
print(
    "Potential screening conflicts:",
    len(conflicts)
)

print("")
print("Output:", OUT)
print("Summary:", SUMMARY)
