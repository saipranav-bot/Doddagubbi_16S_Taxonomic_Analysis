from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
BASE = MASTER / "07_Master_Metadata" / "FINAL_109_PHYLOGENY_V1"
REF = BASE / "REFERENCE_PHYLOGENY_V2"
OUT = BASE / "FINAL_REPORT_V1"

OUT.mkdir(parents=True, exist_ok=True)

MAIN_TREE = BASE / "KJU109_16S.treefile"

REF_TREE = (
    REF / "IQTREE" /
    "KJU109_NCBI65_INDEPENDENT.treefile"
)

ASSESSMENT = (
    REF /
    "FINAL_INTEGRATED_9KJU_TAXONOMIC_ASSESSMENT_V1.csv"
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

# ------------------------------------------------------------
# Taxonomic helpers
# ------------------------------------------------------------

GENUS_ORDER = [
    "Aeromonas",
    "Pseudomonas",
    "Comamonas",
    "Delftia",
    "Alcaligenes",
    "Chromobacterium",
    "Bacillus",
    "Priestia",
    "Peribacillus",
    "Paenibacillus",
    "Staphylococcus",
    "Mammaliicoccus",
    "Listeria",
    "Vagococcus",
    "Exiguobacterium",
    "Citrobacter",
    "Escherichia",
    "Edwardsiella",
    "Pantoea",
    "Shewanella",
    "Pseudoalteromonas",
    "Microbacterium",
    "Leifsonia",
    "Cellulomonas",
    "Kocuria",
    "Glutamicibacter",
    "Micrococcus",
    "Hymenobacter",
    "Limosilactobacillus",
]

def parse_tip(name):
    parts = name.split("|")

    kju = parts[0] if len(parts) >= 1 else ""
    accession = parts[1] if len(parts) >= 2 else ""
    organism = parts[2].replace("_", " ") if len(parts) >= 3 else ""

    genus = organism.split()[0] if organism else ""

    is_kju = kju.startswith("KJU")
    is_ref = "|NCBI_REF|" in name

    return {
        "kju": kju,
        "accession": accession,
        "organism": organism,
        "genus": genus,
        "is_kju": is_kju,
        "is_ref": is_ref,
    }


def safe_name(name):
    return name.replace("|", "_").replace(" ", "_")


# ------------------------------------------------------------
# Tree plotting
# ------------------------------------------------------------

def prepare_tree(tree):
    # Ensure branch lengths exist.
    for clade in tree.find_clades():
        if clade.branch_length is None:
            clade.branch_length = 0.0

    # Remove unnecessary root branch if present.
    tree.root.branch_length = 0.0

    return tree


def get_depths(tree):
    depths = tree.depths()

    max_depth = max(depths.values()) if depths else 1.0
    if max_depth == 0:
        max_depth = 1.0

    return depths, max_depth


def assign_y_positions(tree):
    terminals = tree.get_terminals()
    y = {}

    for i, tip in enumerate(terminals):
        y[tip] = i

    def recurse(clade):
        if clade in y:
            return y[clade]

        child_positions = [recurse(c) for c in clade.clades]
        y[clade] = sum(child_positions) / len(child_positions)
        return y[clade]

    recurse(tree.root)
    return y


def annotate_support(ax, tree, depths, y_positions):
    for clade in tree.find_clades(order="level"):
        if clade.is_terminal():
            continue

        support = clade.confidence

        if support is None:
            continue

        if support < 50:
            continue

        x = depths[clade]
        y = y_positions[clade]

        ax.text(
            x,
            y,
            f"{support:.0f}",
            fontsize=5.5,
            va="center",
            ha="left",
            alpha=0.85,
        )


def plot_tree(
    tree,
    outfile,
    title,
    focus_targets=None,
    reference_mode=False,
    max_label_chars=42,
):
    tree = prepare_tree(tree)

    depths, max_depth = get_depths(tree)
    y_positions = assign_y_positions(tree)

    terminals = tree.get_terminals()
    n = len(terminals)

    height = max(12, n * 0.16)
    width = 20 if n <= 120 else 24

    fig, ax = plt.subplots(
        figsize=(width, height),
        dpi=300,
    )

    # --------------------------------------------------------
    # Draw branches
    # --------------------------------------------------------

    for clade in tree.find_clades(order="preorder"):
        if clade == tree.root:
            continue

        parent = None
        for candidate in tree.find_clades(order="preorder"):
            if clade in candidate.clades:
                parent = candidate
                break

        if parent is None:
            continue

        x0 = depths[parent]
        x1 = depths[clade]
        y = y_positions[clade]
        yp = y_positions[parent]

        ax.plot(
            [x0, x1],
            [y, y],
            linewidth=0.65,
            solid_capstyle="round",
            zorder=1,
        )

    # vertical connectors
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue

        child_y = [y_positions[c] for c in clade.clades]

        if not child_y:
            continue

        x = depths[clade]

        ax.plot(
            [x, x],
            [min(child_y), max(child_y)],
            linewidth=0.65,
            zorder=1,
        )

    # --------------------------------------------------------
    # Draw tips / labels
    # --------------------------------------------------------

    for tip in terminals:
        info = parse_tip(tip.name)

        x = depths[tip]
        y = y_positions[tip]

        label = ""

        # Decide display color
        if focus_targets and info["kju"] in focus_targets:
            label = (
                f"{info['kju']} | "
                f"{info['organism']}"
            )

            color = "crimson"
            weight = "bold"
            fontsize = 8

            ax.scatter(
                [x],
                [y],
                s=13,
                zorder=4,
            )

        elif reference_mode and info["is_ref"]:
            label = (
                f"{info['accession']} | "
                f"{info['organism']}"
            )

            color = "dimgray"
            weight = "normal"
            fontsize = 4.8

        else:
            if info["is_kju"]:
                label = (
                    f"{info['kju']} | "
                    f"{info['organism']}"
                )
            else:
                label = info["organism"]

            color = "black"
            weight = "normal"
            fontsize = 5.5 if n > 120 else 6.5

        if len(label) > max_label_chars:
            label = label[:max_label_chars - 3] + "..."

        ax.text(
            x + max_depth * 0.002,
            y,
            label,
            fontsize=fontsize,
            color=color,
            fontweight=weight,
            va="center",
            ha="left",
            zorder=5,
        )

    # --------------------------------------------------------
    # Support values
    # --------------------------------------------------------

    annotate_support(
        ax,
        tree,
        depths,
        y_positions,
    )

    ax.set_ylim(-1, n)
    ax.set_xlim(
        0,
        max_depth * 1.28,
    )

    ax.invert_yaxis()

    ax.set_xlabel(
        "Evolutionary distance (substitutions/site)",
        fontsize=9,
    )

    ax.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        pad=18,
    )

    ax.set_yticks([])

    # clean frame
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.grid(False)

    plt.tight_layout()
    fig.savefig(
        outfile,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    print("Created:", outfile)


# ------------------------------------------------------------
# Figure 1
# ------------------------------------------------------------

print("=" * 100)
print("FIGURE 1 - MAIN 109-KJU PHYLOGENY")
print("=" * 100)

main_tree = Phylo.read(
    MAIN_TREE,
    "newick",
)

plot_tree(
    main_tree,
    OUT / "FIGURE_1_MAIN_109_KJU_PHYLOGENY.png",
    title=(
        "Doddagubbi 109-KJU 16S rRNA Maximum-Likelihood Phylogeny\n"
        "GTR+F+R4 | 1,000 SH-aLRT + 1,000 UFBoot"
    ),
)


# ------------------------------------------------------------
# Figure 2
# ------------------------------------------------------------

print("=" * 100)
print("FIGURE 2 - INDEPENDENT REFERENCE PHYLOGENY")
print("=" * 100)

ref_tree = Phylo.read(
    REF_TREE,
    "newick",
)

plot_tree(
    ref_tree,
    OUT / "FIGURE_2_109_KJU_PLUS_65_INDEPENDENT_NCBI_REFERENCES.png",
    title=(
        "Doddagubbi KJU Sequences with 65 Independent NCBI 16S References\n"
        "GTR+F+R3 | 13 Target Species | 1,000 SH-aLRT + 1,000 UFBoot"
    ),
    focus_targets=TARGETS,
    reference_mode=True,
    max_label_chars=55,
)


# ------------------------------------------------------------
# Figure 3
#
# Keep the 9 targets + all their nearest 5 reference taxa.
# This gives a compact taxonomic-review tree.
# ------------------------------------------------------------

print("=" * 100)
print("FIGURE 3 - FOCUSED 9-KJU TAXONOMIC REVIEW")
print("=" * 100)

# Build distance dictionary based on the reference tree.
all_terminals = ref_tree.get_terminals()

target_tips = []
reference_tips = []

for tip in all_terminals:
    info = parse_tip(tip.name)

    if info["kju"] in TARGETS:
        target_tips.append(tip)

    if info["is_ref"]:
        reference_tips.append(tip)

# Select nearest reference taxa for each target.
selected_names = set()

for target in target_tips:
    distances = []

    for ref_tip in reference_tips:
        try:
            d = ref_tree.distance(target, ref_tip)
            distances.append((d, ref_tip))
        except Exception:
            pass

    distances.sort(key=lambda x: x[0])

    for _, ref_tip in distances[:5]:
        selected_names.add(ref_tip.name)

selected_names.update(
    tip.name for tip in target_tips
)

# Prune a copy
focus_tree = Phylo.read(
    REF_TREE,
    "newick",
)

for tip in list(focus_tree.get_terminals()):
    if tip.name not in selected_names:
        focus_tree.prune(tip)

plot_tree(
    focus_tree,
    OUT / "FIGURE_3_FOCUSED_9_KJU_TAXONOMIC_REVIEW.png",
    title=(
        "Focused Taxonomic Review of the 9 Discordant KJUs\n"
        "Highlighted KJUs compared with nearest independent NCBI references"
    ),
    focus_targets=TARGETS,
    reference_mode=True,
    max_label_chars=60,
)


# ------------------------------------------------------------
# Generate figure manifest
# ------------------------------------------------------------

manifest = OUT / "FINAL_PHYLOGENY_FIGURE_MANIFEST.txt"

manifest.write_text(
    "\n".join([
        "DODDAGGUBBI FINAL PHYLOGENY FIGURE PACKAGE",
        "=" * 80,
        "",
        "FIGURE 1",
        "FIGURE_1_MAIN_109_KJU_PHYLOGENY.png",
        "109 KJUs; main maximum-likelihood phylogeny",
        "Model: GTR+F+R4",
        "Support: 1000 SH-aLRT + 1000 UFBoot",
        "",
        "FIGURE 2",
        "FIGURE_2_109_KJU_PLUS_65_INDEPENDENT_NCBI_REFERENCES.png",
        "109 KJUs + 65 independent NCBI references",
        "13 target species",
        "Model: GTR+F+R3",
        "Support: 1000 SH-aLRT + 1000 UFBoot",
        "",
        "FIGURE 3",
        "FIGURE_3_FOCUSED_9_KJU_TAXONOMIC_REVIEW.png",
        "9 discordant KJUs + nearest independent reference taxa",
        "",
        "NO SOURCE DATA WERE MODIFIED.",
        "",
    ]) + "\n"
)

print("")
print("=" * 100)
print("FINAL FIGURE GENERATION COMPLETE")
print("=" * 100)
print("Output directory:", OUT)
print("Manifest:", manifest)
print("")
print("Created files:")
for f in sorted(OUT.glob("FIGURE_*.png")):
    print(" ", f)
print("")
print("SOURCE DATA WAS NOT MODIFIED.")
