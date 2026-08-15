from pathlib import Path
from Bio import Phylo
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata" / "FINAL_109_PHYLOGENY_V1"

MAIN = META / "REFERENCE_PHYLOGENY_V2/IQTREE/KJU109_NCBI65_INDEPENDENT.treefile"
SENS = META / "REFERENCE_PHYLOGENY_V2/SENSITIVITY/KJU170_no_bias.treefile"

OUT = META / "REFERENCE_PHYLOGENY_V2/TARGET_SPECIES_MAIN_VS_SENSITIVITY_V1.csv"
SUMMARY = META / "REFERENCE_PHYLOGENY_V2/TARGET_SPECIES_MAIN_VS_SENSITIVITY_SUMMARY_V1.txt"

TARGETS = {
    "Aeromonas caviae",
    "Aeromonas hydrophila",
    "Aeromonas veronii",
    "Bacillus cereus",
    "Citrobacter freundii",
    "Shewanella putrefaciens",
}

def parse_tip(name):
    parts = name.split("|")
    accession = parts[1] if len(parts) > 1 else ""
    organism = parts[2].replace("_", " ") if len(parts) > 2 else ""

    if parts[0].startswith("KJU"):
        kind = "KJU"
        kju = parts[0]
    else:
        kind = "NCBI_REF"
        kju = ""

    return {
        "name": name,
        "kind": kind,
        "kju": kju,
        "accession": accession,
        "organism": organism,
        "species": organism,
    }

def build_tree_info(tree):
    info = {}
    for tip in tree.get_terminals():
        p = parse_tip(tip.name)
        p["tip"] = tip
        info[tip.name] = p
    return info

def target_kjus(info):
    return [
        p for p in info.values()
        if p["kind"] == "KJU" and p["species"] in TARGETS
    ]

def nearest_reference(tree, tip):
    best = None

    for ref in tree.get_terminals():
        p = parse_tip(ref.name)

        if p["kind"] != "NCBI_REF":
            continue

        try:
            d = tree.distance(tip, ref)
        except Exception:
            continue

        if best is None or d < best["distance"]:
            best = {
                "reference": ref.name,
                "reference_species": p["species"],
                "distance": d,
            }

    return best

def nearest_any_species(tree, tip, exclude=None):
    best = None

    for other in tree.get_terminals():
        if other.name == tip.name:
            continue

        p = parse_tip(other.name)

        if p["kind"] != "KJU":
            continue

        if exclude and p["species"] == exclude:
            continue

        try:
            d = tree.distance(tip, other)
        except Exception:
            continue

        if best is None or d < best["distance"]:
            best = {
                "kju": p["kju"],
                "species": p["species"],
                "distance": d,
            }

    return best

def same_species_reference_distance(tree, tip, species):
    distances = []

    for ref in tree.get_terminals():
        p = parse_tip(ref.name)

        if p["kind"] != "NCBI_REF":
            continue

        if p["species"] != species:
            continue

        try:
            distances.append(tree.distance(tip, ref))
        except Exception:
            pass

    if not distances:
        return None

    return min(distances)

def summarize_tree(tree, label):
    info = build_tree_info(tree)
    rows = []

    for p in target_kjus(info):
        tip = p["tip"]

        same_ref = same_species_reference_distance(
            tree,
            tip,
            p["species"]
        )

        nearest_ref = nearest_reference(tree, tip)

        nearest_other = nearest_any_species(
            tree,
            tip,
            exclude=p["species"]
        )

        rows.append({
            "analysis": label,
            "kju_id": p["kju"],
            "accession": p["accession"],
            "annotated_species": p["species"],
            "nearest_independent_reference":
                nearest_ref["reference"] if nearest_ref else "",
            "nearest_reference_species":
                nearest_ref["reference_species"] if nearest_ref else "",
            "nearest_reference_distance":
                nearest_ref["distance"] if nearest_ref else "",
            "nearest_other_kju":
                nearest_other["kju"] if nearest_other else "",
            "nearest_other_species":
                nearest_other["species"] if nearest_other else "",
            "nearest_other_distance":
                nearest_other["distance"] if nearest_other else "",
            "nearest_same_species_reference_distance":
                same_ref if same_ref is not None else "",
            "reference_supported":
                (
                    nearest_ref["reference_species"] == p["species"]
                    if nearest_ref else False
                )
        })

    return rows

main_tree = Phylo.read(MAIN, "newick")
sens_tree = Phylo.read(SENS, "newick")

rows = []
rows.extend(summarize_tree(main_tree, "MAIN"))
rows.extend(summarize_tree(sens_tree, "SENSITIVITY"))

df = pd.DataFrame(rows)

# ------------------------------------------------------------
# Direct comparison
# ------------------------------------------------------------

main_df = df[df["analysis"] == "MAIN"].copy()
sens_df = df[df["analysis"] == "SENSITIVITY"].copy()

main_df = main_df.set_index("kju_id")
sens_df = sens_df.set_index("kju_id")

comparison = []

for kju in sorted(set(main_df.index) & set(sens_df.index)):

    m = main_df.loc[kju]
    s = sens_df.loc[kju]

    comparison.append({
        "kju_id": kju,
        "annotated_species": m["annotated_species"],

        "main_nearest_reference":
            m["nearest_reference_species"],

        "sensitivity_nearest_reference":
            s["nearest_reference_species"],

        "main_nearest_reference_distance":
            m["nearest_reference_distance"],

        "sensitivity_nearest_reference_distance":
            s["nearest_reference_distance"],

        "main_nearest_other_species":
            m["nearest_other_species"],

        "sensitivity_nearest_other_species":
            s["nearest_other_species"],

        "main_reference_supported":
            m["reference_supported"],

        "sensitivity_reference_supported":
            s["reference_supported"],

        "placement_stable":
            (
                m["nearest_reference_species"]
                == s["nearest_reference_species"]
            ),

        "annotation_supported_in_both":
            (
                m["reference_supported"]
                and s["reference_supported"]
            )
    })

comp = pd.DataFrame(comparison)

OUT.parent.mkdir(parents=True, exist_ok=True)

comp.to_csv(OUT, index=False)

with open(SUMMARY, "w") as fh:

    fh.write("=" * 100 + "\n")
    fh.write("DODDAGGUBBI - TARGET SPECIES MAIN VS SENSITIVITY\n")
    fh.write("=" * 100 + "\n\n")

    fh.write(f"Target species: {len(TARGETS)}\n")
    fh.write(f"KJUs compared: {len(comp)}\n\n")

    fh.write("PLACEMENT STABILITY\n")
    fh.write("-" * 100 + "\n")

    stable = int(comp["placement_stable"].sum())
    unstable = len(comp) - stable

    fh.write(f"Stable placements: {stable}\n")
    fh.write(f"Changed placements: {unstable}\n\n")

    fh.write("ANNOTATION SUPPORT\n")
    fh.write("-" * 100 + "\n")

    supported = int(comp["annotation_supported_in_both"].sum())
    not_supported = len(comp) - supported

    fh.write(
        f"Supported by nearest independent reference in both trees: "
        f"{supported}\n"
    )
    fh.write(
        f"Not supported in both trees: {not_supported}\n\n"
    )

    fh.write("DETAILS\n")
    fh.write("-" * 100 + "\n")

    fh.write(comp.to_string(index=False))

print("=" * 100)
print("TARGET SPECIES MAIN VS SENSITIVITY COMPLETE")
print("=" * 100)

print("KJUs compared:", len(comp))
print(
    "Stable nearest-reference placements:",
    int(comp["placement_stable"].sum())
)
print(
    "Changed placements:",
    len(comp) - int(comp["placement_stable"].sum())
)

print("")
print(comp.to_string(index=False))

print("")
print("Output:", OUT)
print("Summary:", SUMMARY)
