from pathlib import Path
from Bio import Phylo
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata" / "FINAL_109_PHYLOGENY_V1"

TREE = (
    META /
    "REFERENCE_PHYLOGENY_V2/IQTREE/"
    "KJU109_NCBI65_INDEPENDENT.treefile"
)

CLASSIFICATION = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_KJU_CLASSIFICATION_V1.csv"
)

OUT = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_CLade_BOOTSTRAP_REVIEW_V1.csv"
)

SUMMARY = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_CLADE_BOOTSTRAP_REVIEW_SUMMARY_V1.txt"
)

TARGETS = {
    "KJU0007",
    "KJU0017",
    "KJU0088",
    "KJU0099",
    "KJU0105",
    "KJU0106",
    "KJU0208",
    "KJU0217",
    "KJU0221",
}

tree = Phylo.read(TREE, "newick")
classification = pd.read_csv(CLASSIFICATION, dtype=str).fillna("")


def parse_tip(name):

    parts = name.split("|")

    if len(parts) < 3:
        return "", "", ""

    return (
        parts[0],
        parts[1],
        parts[2].replace("_", " ")
    )


# ------------------------------------------------------------
# Build tip metadata
# ------------------------------------------------------------

tips = []

for terminal in tree.get_terminals():

    kju, accession, organism = parse_tip(
        terminal.name
    )

    if not organism:
        continue

    tips.append({
        "terminal": terminal,
        "name": terminal.name,
        "kju": kju,
        "accession": accession,
        "organism": organism,
        "species": organism.lower(),
        "is_kju": kju.startswith("KJU"),
        "is_reference": "|NCBI_REF|" in terminal.name
    })


rows = []


# ------------------------------------------------------------
# Examine target KJUs
# ------------------------------------------------------------

for target in sorted(TARGETS):

    target_info = next(
        x for x in tips
        if x["kju"] == target
    )

    annotation = classification.loc[
        classification["kju_id"] == target,
        "current_annotation"
    ].iloc[0]

    alternative = classification.loc[
        classification["kju_id"] == target,
        "best_alternative_species"
    ].iloc[0]

    # Find terminal
    terminal = target_info["terminal"]

    # --------------------------------------------------------
    # Collect path towards root
    # --------------------------------------------------------

    path = tree.get_path(terminal)

    candidate_clades = []

    # Inspect ancestors from nearest to farthest
    current = tree.root

    for node in path:

        terminals = node.get_terminals()

        names = [
            t.name
            for t in terminals
        ]

        if terminal.name in names:

            organisms = []

            for t in terminals:

                parts = parse_tip(t.name)

                if len(parts) >= 3:
                    organisms.append(
                        parts[2].lower()
                    )

            annotation_count = sum(
                1
                for x in organisms
                if annotation.lower() in x
            )

            alternative_count = sum(
                1
                for x in organisms
                if alternative.lower() in x
            )

            # IQ-TREE uses SH-aLRT / UFBoot annotation
            # on internal branches in the Newick file.
            support = node.confidence

            candidate_clades.append({
                "node": node,
                "size": len(terminals),
                "annotation_count": annotation_count,
                "alternative_count": alternative_count,
                "support": support
            })

    # --------------------------------------------------------
    # Find smallest relevant clade containing target
    # and alternative-species members
    # --------------------------------------------------------

    relevant = [
        x for x in candidate_clades
        if x["alternative_count"] > 0
    ]

    if relevant:

        best = min(
            relevant,
            key=lambda x: x["size"]
        )

        clade_size = best["size"]
        annotation_count = best["annotation_count"]
        alternative_count = best["alternative_count"]
        support = best["support"]

    else:

        clade_size = 1
        annotation_count = 1
        alternative_count = 0
        support = None

    # --------------------------------------------------------
    # Classification
    # --------------------------------------------------------

    if alternative_count == 0:

        interpretation = (
            "ALTERNATIVE_NOT_IN_TARGET_CLADE"
        )

    elif support is None:

        interpretation = (
            "ALTERNATIVE_CLADE_NO_SUPPORT_VALUE"
        )

    elif support >= 95:

        interpretation = (
            "STRONG_SUPPORTED_ALTERNATIVE_CLADE"
        )

    elif support >= 80:

        interpretation = (
            "MODERATELY_SUPPORTED_ALTERNATIVE_CLADE"
        )

    elif support >= 70:

        interpretation = (
            "WEAK_ALTERNATIVE_CLADE"
        )

    else:

        interpretation = (
            "UNSUPPORTED_OR_LOW_SUPPORT_ALTERNATIVE"
        )

    rows.append({

        "kju_id": target,

        "current_annotation":
            annotation,

        "alternative_species":
            alternative,

        "target_annotation_count_in_relevant_clade":
            annotation_count,

        "alternative_species_count_in_relevant_clade":
            alternative_count,

        "relevant_clade_size":
            clade_size,

        "clade_support_percent":
            support
            if support is not None
            else "",

        "interpretation":
            interpretation
    })


out = pd.DataFrame(rows)

out.to_csv(
    OUT,
    index=False
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

with open(SUMMARY, "w") as f:

    f.write("=" * 110 + "\n")
    f.write(
        "DODDAGGUBBI - BOOTSTRAP-SUPPORTED "
        "DISCORDANT CLADE REVIEW\n"
    )
    f.write("=" * 110 + "\n\n")

    f.write(
        f"Discordant KJUs reviewed: {len(out)}\n\n"
    )

    f.write("CLASSIFICATION COUNTS\n")
    f.write("-" * 110 + "\n")

    for cls, count in (
        out["interpretation"]
        .value_counts()
        .items()
    ):

        f.write(
            f"{cls}: {count}\n"
        )

    f.write("\n\n")

    f.write("DETAILED RESULTS\n")
    f.write("-" * 110 + "\n")

    for _, r in out.iterrows():

        f.write(
            f"{r['kju_id']} | "
            f"{r['current_annotation']} | "
            f"alternative={r['alternative_species']} | "
            f"clade_size={r['relevant_clade_size']} | "
            f"support={r['clade_support_percent']} | "
            f"{r['interpretation']}\n"
        )


print("=" * 110)
print(
    "BOOTSTRAP-SUPPORTED DISCORDANT CLADE REVIEW COMPLETE"
)
print("=" * 110)

print(
    out.to_string(index=False)
)

print("\nOutput:", OUT)
print("Summary:", SUMMARY)
