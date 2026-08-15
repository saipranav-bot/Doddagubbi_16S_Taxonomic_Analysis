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

OUT = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_KJU_REFERENCE_REVIEW_V1.csv"
)

SUMMARY = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_KJU_REFERENCE_REVIEW_SUMMARY_V1.txt"
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

def parse(name):
    parts = name.split("|")

    kju = parts[0] if parts and parts[0].startswith("KJU") else ""
    accession = parts[1] if len(parts) > 1 else ""
    organism = parts[2].replace("_", " ") if len(parts) > 2 else ""

    return {
        "name": name,
        "kju": kju,
        "accession": accession,
        "organism": organism
    }

tips = tree.get_terminals()

parsed = {
    tip.name: parse(tip.name)
    for tip in tips
}

kju_tips = {
    p["kju"]: tip
    for tip, p in zip(tips, parsed.values())
    if p["kju"]
}

reference_tips = [
    tip for tip in tips
    if "|NCBI_REF|" in tip.name
]

rows = []

for kju in sorted(TARGETS):

    if kju not in kju_tips:
        continue

    tip = kju_tips[kju]
    current = parsed[tip.name]

    distances = []

    for ref in reference_tips:

        ref_info = parsed[ref.name]

        d = tree.distance(tip, ref)

        distances.append({
            "kju_id": kju,
            "current_annotation": current["organism"],
            "reference_accession": ref_info["accession"],
            "reference_species": ref_info["organism"],
            "distance": d
        })

    distances.sort(key=lambda x: x["distance"])

    # Keep closest 5 independent references
    for rank, x in enumerate(distances[:5], start=1):
        x["reference_rank"] = rank
        rows.append(x)

df = pd.DataFrame(rows)

OUT.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(OUT, index=False)

with open(SUMMARY, "w") as fh:

    fh.write("=" * 100 + "\n")
    fh.write("DODDAGGUBBI - DISCORDANT KJU INDEPENDENT REFERENCE REVIEW\n")
    fh.write("=" * 100 + "\n\n")

    for kju in sorted(TARGETS):

        sub = df[df["kju_id"] == kju]

        if sub.empty:
            continue

        annotation = sub["current_annotation"].iloc[0]

        fh.write(f"{kju} | {annotation}\n")
        fh.write("-" * 100 + "\n")

        for _, r in sub.iterrows():

            fh.write(
                f"Rank {int(r['reference_rank'])} | "
                f"{r['reference_species']} | "
                f"{r['reference_accession']} | "
                f"distance={r['distance']:.8f}\n"
            )

        fh.write("\n")

print("=" * 100)
print("DISCORDANT KJU REFERENCE REVIEW COMPLETE")
print("=" * 100)

for kju in sorted(TARGETS):

    sub = df[df["kju_id"] == kju]

    if sub.empty:
        continue

    print(f"\n{kju} | annotation={sub['current_annotation'].iloc[0]}")

    for _, r in sub.iterrows():

        print(
            f"  {int(r['reference_rank'])}. "
            f"{r['reference_species']} "
            f"({r['reference_accession']}) "
            f"distance={r['distance']:.8f}"
        )

print("\nOutput:", OUT)
print("Summary:", SUMMARY)
