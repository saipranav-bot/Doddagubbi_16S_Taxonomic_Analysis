from pathlib import Path
from Bio import Phylo
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

TREE = META / "FINAL_109_PHYLOGENY_V1" / "KJU109_16S.treefile"
QC = META / "FINAL_109_SEQUENCE_QC_V2" / "FINAL_109_SEQUENCE_QC_V2.csv"

OUT = META / "FINAL_109_PHYLOGENY_V1" / \
      "GENUS_CLADE_CONSISTENCY_V1.csv"

SUMMARY = META / "FINAL_109_PHYLOGENY_V1" / \
          "GENUS_CLADE_CONSISTENCY_SUMMARY_V1.txt"

tree = Phylo.read(TREE, "newick")
qc = pd.read_csv(QC, dtype=str).fillna("")

qc["genus"] = qc["organism"].str.split().str[0]

tip_map = {}

for tip in tree.get_terminals():
    parts = tip.name.split("|")
    if len(parts) < 3:
        continue

    kju = parts[0]
    organism = parts[2].replace("_", " ")
    genus = organism.split()[0] if organism else ""

    tip_map[kju] = {
        "tip": tip,
        "genus": genus,
        "organism": organism
    }

rows = []

for _, r in qc.iterrows():

    kju = r["kju_id"]
    genus = r["genus"]
    organism = r["organism"]

    if kju not in tip_map:
        continue

    tip = tip_map[kju]["tip"]

    same_genus = []

    for other_kju, info in tip_map.items():

        if other_kju == kju:
            continue

        if info["genus"] != genus:
            continue

        d = tree.distance(
            tip,
            info["tip"]
        )

        same_genus.append(
            (other_kju, d)
        )

    same_genus.sort(
        key=lambda x: x[1]
    )

    nearest_same_genus = (
        same_genus[0][0]
        if same_genus
        else ""
    )

    nearest_same_distance = (
        same_genus[0][1]
        if same_genus
        else ""
    )

    # Find immediate parent clade
    path = tree.get_path(tip)

    genus_clade_size = 0

    if path:

        # Examine ancestors from closest to furthest
        for ancestor in reversed(path):

            terminals = ancestor.get_terminals()

            taxa_genera = []

            for terminal in terminals:

                name = terminal.name
                p = name.split("|")

                if len(p) >= 3:
                    g = p[2].replace("_", " ").split()[0]
                    taxa_genera.append(g)

            if taxa_genera and all(
                g == genus
                for g in taxa_genera
            ):
                genus_clade_size = len(terminals)
                break

    if genus_clade_size > 0:
        placement = "GENUS_CLADE_SUPPORTED"
    elif same_genus:
        placement = "SAME_GENUS_PRESENT_BUT_NOT_MONOPHYLETIC"
    else:
        placement = "NO_OTHER_GENUS_TAXA"

    rows.append({
        "kju_id": kju,
        "organism": organism,
        "genus": genus,
        "same_genus_taxa_count": len(same_genus),
        "nearest_same_genus_kju": nearest_same_genus,
        "nearest_same_genus_distance":
            nearest_same_distance,
        "genus_clade_size": genus_clade_size,
        "placement_class": placement
    })

result = pd.DataFrame(rows)

result.to_csv(
    OUT,
    index=False
)

summary = [
    "DODDAGGUBBI - GENUS CLADE CONSISTENCY",
    "=" * 100,
    "",
    f"Taxa analyzed: {len(result)}",
    "",
    "PLACEMENT COUNTS",
    "-" * 100,
]

for key, value in result["placement_class"].value_counts().items():
    summary.append(
        f"{key}: {value}"
    )

summary += [
    "",
    "IMPORTANT",
    "-" * 100,
    "This analysis is a screening layer.",
    "It does not automatically establish misidentification.",
    "Final taxonomic decisions should consider clade support,",
    "reference taxon sampling, and additional sequence evidence.",
]

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("=" * 100)
print("GENUS CLADE CONSISTENCY COMPLETE")
print("=" * 100)

print(result["placement_class"].value_counts())

print("")
print("Output:", OUT)
print("Summary:", SUMMARY)
