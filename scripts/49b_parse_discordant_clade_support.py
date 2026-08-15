from pathlib import Path
from Bio import Phylo
import pandas as pd
import re

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
    "DISCORDANT_CLADE_SUPPORT_REVIEW_V2.csv"
)

SUMMARY = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_CLADE_SUPPORT_REVIEW_SUMMARY_V2.txt"
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


def parse_tip(name):
    parts = name.split("|")

    if len(parts) < 3:
        return "", "", ""

    return (
        parts[0],
        parts[1],
        parts[2].replace("_", " ")
    )


def extract_support(clade):
    """
    IQ-TREE may store support as either:
        clade.confidence
    or:
        clade.name = '98/100'
        clade.name = '98.5/100'
        clade.name = '98'
        etc.

    Returns:
        sh_alrt, ufboot
    """

    # Case 1: Biopython parsed a numeric confidence
    if clade.confidence is not None:

        try:
            value = float(clade.confidence)

            return value, None

        except Exception:
            pass

    # Case 2: support is stored in internal node name
    label = str(clade.name or "").strip()

    if not label:
        return None, None

    # Common IQ-TREE form:
    # 98/100
    m = re.match(
        r"^\s*([0-9]+(?:\.[0-9]+)?)\s*/\s*"
        r"([0-9]+(?:\.[0-9]+)?)\s*$",
        label
    )

    if m:

        sh_alrt = float(m.group(1))
        ufboot = float(m.group(2))

        return sh_alrt, ufboot

    # Single support number
    m = re.match(
        r"^\s*([0-9]+(?:\.[0-9]+)?)\s*$",
        label
    )

    if m:

        return float(m.group(1)), None

    return None, None


tree = Phylo.read(TREE, "newick")
classification = pd.read_csv(
    CLASSIFICATION,
    dtype=str
).fillna("")

rows = []

for target in sorted(TARGETS):

    target_terminal = next(
        t for t in tree.get_terminals()
        if t.name.startswith(target + "|")
    )

    annotation = classification.loc[
        classification["kju_id"] == target,
        "current_annotation"
    ].iloc[0]

    alternative = classification.loc[
        classification["kju_id"] == target,
        "best_alternative_species"
    ].iloc[0]

    target_genus = annotation.split()[0]
    alternative_genus = alternative.split()[0]

    # --------------------------------------------------------
    # Find clade containing target + alternative species
    # --------------------------------------------------------

    path = tree.get_path(target_terminal)

    candidates = []

    for node in path:

        terminals = node.get_terminals()

        organisms = []

        for t in terminals:

            _, _, organism = parse_tip(t.name)

            if organism:
                organisms.append(organism)

        alternative_count = sum(
            1
            for x in organisms
            if alternative.lower() == x.lower()
        )

        annotation_count = sum(
            1
            for x in organisms
            if annotation.lower() == x.lower()
        )

        if alternative_count > 0:

            sh_alrt, ufboot = extract_support(node)

            candidates.append({
                "node": node,
                "size": len(terminals),
                "annotation_count": annotation_count,
                "alternative_count": alternative_count,
                "sh_alrt": sh_alrt,
                "ufboot": ufboot,
                "label": str(node.name or "")
            })

    if not candidates:

        rows.append({
            "kju_id": target,
            "current_annotation": annotation,
            "alternative_species": alternative,
            "clade_size": "",
            "annotation_count": "",
            "alternative_count": "",
            "node_label": "",
            "sh_alrt_percent": "",
            "ufboot_percent": "",
            "classification":
                "ALTERNATIVE_NOT_IN_RELEVANT_CLADE"
        })

        continue

    # Smallest clade containing the alternative species
    best = min(
        candidates,
        key=lambda x: x["size"]
    )

    sh = best["sh_alrt"]
    ub = best["ufboot"]

    # Prefer UFBoot when available
    support = ub if ub is not None else sh

    if support is None:

        classification_result = (
            "SUPPORT_NOT_PARSED"
        )

    elif support >= 95:

        classification_result = (
            "STRONG_SUPPORTED_ALTERNATIVE_CLADE"
        )

    elif support >= 80:

        classification_result = (
            "MODERATELY_SUPPORTED_ALTERNATIVE_CLADE"
        )

    elif support >= 70:

        classification_result = (
            "WEAK_SUPPORTED_ALTERNATIVE_CLADE"
        )

    else:

        classification_result = (
            "LOW_SUPPORT_ALTERNATIVE_CLADE"
        )

    rows.append({
        "kju_id": target,
        "current_annotation": annotation,
        "alternative_species": alternative,
        "clade_size": best["size"],
        "annotation_count": best["annotation_count"],
        "alternative_count": best["alternative_count"],
        "node_label": best["label"],
        "sh_alrt_percent":
            "" if sh is None else sh,
        "ufboot_percent":
            "" if ub is None else ub,
        "classification":
            classification_result
    })


out = pd.DataFrame(rows)

out.to_csv(
    OUT,
    index=False
)

with open(SUMMARY, "w") as f:

    f.write("=" * 110 + "\n")
    f.write(
        "DODDAGGUBBI - DISCORDANT CLADE SUPPORT REVIEW V2\n"
    )
    f.write("=" * 110 + "\n\n")

    f.write(
        f"Targets reviewed: {len(out)}\n\n"
    )

    f.write("CLASSIFICATION COUNTS\n")
    f.write("-" * 110 + "\n")

    for cls, count in out["classification"].value_counts().items():

        f.write(
            f"{cls}: {count}\n"
        )

    f.write("\nDETAILED RESULTS\n")
    f.write("-" * 110 + "\n")

    for _, r in out.iterrows():

        f.write(
            f"{r['kju_id']} | "
            f"{r['current_annotation']} | "
            f"alternative={r['alternative_species']} | "
            f"clade_size={r['clade_size']} | "
            f"SH-aLRT={r['sh_alrt_percent']} | "
            f"UFBoot={r['ufboot_percent']} | "
            f"{r['classification']}\n"
        )


print("=" * 110)
print("DISCORDANT CLADE SUPPORT REVIEW V2 COMPLETE")
print("=" * 110)

print(out.to_string(index=False))

print("\nOutput:", OUT)
print("Summary:", SUMMARY)
