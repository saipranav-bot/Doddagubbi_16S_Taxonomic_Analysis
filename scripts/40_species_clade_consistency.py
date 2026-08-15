from pathlib import Path
from Bio import Phylo
import pandas as pd

# ============================================================
# DODDAGGUBBI
# SPECIES-LEVEL PHYLOGENETIC CONSISTENCY
#
# SAFE / NON-DESTRUCTIVE
# No source files are modified.
# ============================================================

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

TREE = (
    META
    / "FINAL_109_PHYLOGENY_V1"
    / "KJU109_16S.treefile"
)

QC = (
    META
    / "FINAL_109_SEQUENCE_QC_V2"
    / "FINAL_109_SEQUENCE_QC_V2.csv"
)

OUT = (
    META
    / "FINAL_109_PHYLOGENY_V1"
    / "SPECIES_CLADE_CONSISTENCY_V1.csv"
)

SUMMARY = (
    META
    / "FINAL_109_PHYLOGENY_V1"
    / "SPECIES_CLADE_CONSISTENCY_SUMMARY_V1.txt"
)

print("=" * 110)
print(" DODDAGGUBBI - SPECIES-LEVEL PHYLOGENETIC CONSISTENCY")
print("=" * 110)

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

if not TREE.exists():
    raise SystemExit(f"ERROR: Missing tree:\n{TREE}")

if not QC.exists():
    raise SystemExit(f"ERROR: Missing QC table:\n{QC}")

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------

tree = Phylo.read(TREE, "newick")
qc = pd.read_csv(QC, dtype=str).fillna("")

if "kju_id" not in qc.columns:
    raise SystemExit("ERROR: QC file has no kju_id column.")

if "organism" not in qc.columns:
    raise SystemExit("ERROR: QC file has no organism column.")

# ------------------------------------------------------------
# Parse species
# ------------------------------------------------------------

def species_name(organism):
    parts = str(organism).strip().split()

    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"

    if len(parts) == 1:
        return parts[0]

    return ""

qc["annotated_species"] = qc["organism"].map(species_name)
qc["annotated_genus"] = (
    qc["organism"]
    .astype(str)
    .str.strip()
    .str.split()
    .str[0]
)

# ------------------------------------------------------------
# Build tree tip map
# ------------------------------------------------------------

tip_map = {}

for tip in tree.get_terminals():

    parts = str(tip.name).split("|")

    if len(parts) < 3:
        continue

    kju_id = parts[0].strip()
    accession = parts[1].strip()
    organism = parts[2].replace("_", " ").strip()

    tip_species = species_name(organism)

    tip_genus = (
        organism.split()[0]
        if organism.split()
        else ""
    )

    tip_map[kju_id] = {
        "tip": tip,
        "accession": accession,
        "organism": organism,
        "species": tip_species,
        "genus": tip_genus,
    }

print(f"Tree tips: {len(tree.get_terminals())}")
print(f"QC records: {len(qc)}")

# ------------------------------------------------------------
# Analyze each KJU
# ------------------------------------------------------------

rows = []

for _, row in qc.iterrows():

    kju = row["kju_id"].strip()

    if kju not in tip_map:
        continue

    info = tip_map[kju]
    tip = info["tip"]

    annotated_species = row["annotated_species"]
    annotated_genus = row["annotated_genus"]

    same_species = []
    same_genus = []

    # --------------------------------------------------------
    # Calculate distances
    # --------------------------------------------------------

    for other_kju, other_info in tip_map.items():

        if other_kju == kju:
            continue

        try:
            distance = tree.distance(
                tip,
                other_info["tip"]
            )
        except Exception:
            continue

        if other_info["species"] == annotated_species:
            same_species.append(
                (
                    other_kju,
                    other_info["organism"],
                    distance
                )
            )

        if other_info["genus"] == annotated_genus:
            same_genus.append(
                (
                    other_kju,
                    other_info["organism"],
                    distance
                )
            )

    same_species.sort(key=lambda x: x[2])
    same_genus.sort(key=lambda x: x[2])

    # --------------------------------------------------------
    # Species group
    # --------------------------------------------------------

    species_members = [
        info2["tip"]
        for info2 in tip_map.values()
        if info2["species"] == annotated_species
    ]

    species_group_count = len(species_members)

    # --------------------------------------------------------
    # Determine whether all same-species taxa form
    # an exclusive clade with this sequence.
    # --------------------------------------------------------

    species_clade_supported = False
    species_clade_size = 0

    if species_group_count >= 2:

        try:

            ancestor = tree.common_ancestor(
                species_members
            )

            terminals = ancestor.get_terminals()

            terminal_species = []

            for terminal in terminals:

                p = str(terminal.name).split("|")

                if len(p) < 3:
                    continue

                organism = (
                    p[2]
                    .replace("_", " ")
                    .strip()
                )

                terminal_species.append(
                    species_name(organism)
                )

            if terminal_species:

                if all(
                    x == annotated_species
                    for x in terminal_species
                ):

                    species_clade_supported = True
                    species_clade_size = len(terminals)

        except Exception:
            pass

    # --------------------------------------------------------
    # Classify
    # --------------------------------------------------------

    if species_group_count == 1:

        classification = "SINGLETON_SPECIES"

    elif species_clade_supported:

        classification = "SPECIES_CLADE_SUPPORTED"

    else:

        classification = (
            "SAME_SPECIES_PRESENT_NOT_EXCLUSIVE_CLADE"
        )

    # --------------------------------------------------------
    # Nearest same-species
    # --------------------------------------------------------

    nearest_species_kju = ""
    nearest_species_org = ""
    nearest_species_distance = ""

    if same_species:

        nearest_species_kju = same_species[0][0]
        nearest_species_org = same_species[0][1]
        nearest_species_distance = same_species[0][2]

    # --------------------------------------------------------
    # Nearest same-genus
    # --------------------------------------------------------

    nearest_genus_kju = ""
    nearest_genus_distance = ""

    if same_genus:

        nearest_genus_kju = same_genus[0][0]
        nearest_genus_distance = same_genus[0][2]

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    rows.append({

        "kju_id": kju,

        "accession": row.get(
            "accession",
            ""
        ),

        "organism": row["organism"],

        "annotated_genus": annotated_genus,

        "annotated_species": annotated_species,

        "same_species_taxa_count":
            species_group_count,

        "same_genus_taxa_count":
            len(same_genus),

        "nearest_same_species_kju":
            nearest_species_kju,

        "nearest_same_species_organism":
            nearest_species_org,

        "nearest_same_species_distance":
            nearest_species_distance,

        "nearest_same_genus_kju":
            nearest_genus_kju,

        "nearest_same_genus_distance":
            nearest_genus_distance,

        "species_clade_size":
            species_clade_size,

        "classification":
            classification
    })

# ------------------------------------------------------------
# Create result table
# ------------------------------------------------------------

result = pd.DataFrame(rows)

# Sort by KJU number
result["_num"] = (
    result["kju_id"]
    .str.extract(r"KJU(\d+)")
    [0]
    .astype(int)
)

result = (
    result
    .sort_values("_num")
    .drop(columns="_num")
)

result.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

counts = (
    result["classification"]
    .value_counts()
)

summary = []

summary.append(
    "DODDAGGUBBI - SPECIES-LEVEL PHYLOGENETIC CONSISTENCY"
)

summary.append("=" * 110)
summary.append("")

summary.append(
    f"Tree tips analyzed: {len(result)}"
)

summary.append(
    f"Species groups represented: "
    f"{result['annotated_species'].nunique()}"
)

summary.append("")
summary.append("CLASSIFICATION")
summary.append("-" * 110)

for label in [
    "SPECIES_CLADE_SUPPORTED",
    "SAME_SPECIES_PRESENT_NOT_EXCLUSIVE_CLADE",
    "SINGLETON_SPECIES",
]:

    summary.append(
        f"{label}: {counts.get(label, 0)}"
    )

summary.append("")
summary.append("INTERPRETATION")
summary.append("-" * 110)

summary.append(
    "SPECIES_CLADE_SUPPORTED means all observed members "
    "of the annotated species form an exclusive clade "
    "in the current tree."
)

summary.append(
    "SAME_SPECIES_PRESENT_NOT_EXCLUSIVE_CLADE means the "
    "same species is present but its observed members "
    "do not form an exclusive clade."
)

summary.append(
    "SINGLETON_SPECIES means only one KJU in this dataset "
    "has that species annotation, so species-level "
    "monophyly cannot be tested within this dataset."
)

summary.append("")
summary.append("IMPORTANT")
summary.append("-" * 110)

summary.append(
    "This is a screening analysis, not proof of "
    "species misidentification."
)

summary.append(
    "16S rRNA may have insufficient resolution for "
    "closely related species."
)

summary.append(
    "Any species-level conflict requires additional "
    "taxonomic/reference evidence before changing metadata."
)

# ------------------------------------------------------------
# Candidate list
# ------------------------------------------------------------

summary.append("")
summary.append("REVIEW CANDIDATES")
summary.append("-" * 110)

candidates = result[
    result["classification"] !=
    "SPECIES_CLADE_SUPPORTED"
]

if candidates.empty:

    summary.append("NONE")

else:

    for _, r in candidates.iterrows():

        summary.append(
            f"{r['kju_id']} | "
            f"{r['organism']} | "
            f"{r['classification']} | "
            f"same_species_count="
            f"{r['same_species_taxa_count']} | "
            f"species_clade_size="
            f"{r['species_clade_size']}"
        )

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

print("")
print("=" * 110)
print(" SPECIES-LEVEL CONSISTENCY COMPLETE")
print("=" * 110)

print("")
print("Classification:")
print(counts)

print("")
print("Output:")
print(OUT)

print("")
print("Summary:")
print(SUMMARY)

print("")
print("NO SOURCE DATA WAS MODIFIED.")
