#!/usr/bin/env python3
"""STEP 107: Phylogenetic + Taxonomic Diversity Integration"""

import hashlib
import re
from pathlib import Path
import pandas as pd

BASE = Path.home() / "Doddagubbi_Master"
OUT = BASE / "downstream_analysis/STEP107_PHYLOGENETIC_TAXONOMIC_DIVERSITY"
OUT.mkdir(parents=True, exist_ok=True)

TREE = BASE / "downstream_analysis/STEP106_REPRODUCIBLE_PHYLOGENY/STEP106_PHYLOGENY.treefile"
META = BASE / "downstream_analysis/STEP102_ANALYSIS_METADATA/STEP102_ANALYSIS_METADATA.csv"

def sha256(p): return hashlib.sha256(p.read_bytes()).hexdigest()

print("STEP 107: PHYLOGENETIC + TAXONOMIC DIVERSITY")
print("=" * 60)

# Load tree as text
tree_text = TREE.read_text()
print(f"Tree length: {len(tree_text)} chars")

# Extract tip labels
tips = re.findall(r"([^,();]+?)(?:[,();])", tree_text)
tips = [t.strip() for t in tips if "|" in t and not t.startswith("(")]
print(f"Tree tips: {len(tips)}")

# Parse tip labels
tip_data = []
for tip in tips:
    parts = tip.split("|")
    if len(parts) >= 3:
        tip_data.append({
            "tip_label": tip,
            "KJU_ID": parts[0],
            "accession": parts[1],
            "taxonomy": parts[2],
        })

tip_df = pd.DataFrame(tip_data)
print(f"Parsed tips: {len(tip_df)}")

# Load metadata
meta = pd.read_csv(META, dtype=str)

# Merge tree tips with metadata
merged = tip_df.merge(meta, on="KJU_ID", how="left", suffixes=("_tree", "_meta"))
merged.to_csv(OUT / "STEP107_TREE_TIP_METADATA.csv", index=False)

# Genus clustering
genus_map = merged.groupby("Genus").agg(
    Count=("KJU_ID", "count"),
    KJU_IDs=("KJU_ID", lambda x: "; ".join(x)),
    Taxonomy=("taxonomy", lambda x: "; ".join(set(x))),
).reset_index()
genus_map = genus_map.sort_values("Count", ascending=False)
genus_map.to_csv(OUT / "STEP107_GENUS_CLUSTERS.csv", index=False)

# Taxonomic diversity metrics
n_genera = merged["Genus"].nunique()
n_species = merged["Species"].nunique()
shannon_genus = -sum(
    (c / len(merged)) * (c / len(merged) > 0 and __import__("math").log(c / len(merged)))
    for c in merged["Genus"].value_counts()
)
shannon_species = -sum(
    (c / len(merged)) * (c / len(merged) > 0 and __import__("math").log(c / len(merged)))
    for c in merged["Species"].value_counts()
)

print(f"Genera: {n_genera}")
print(f"Species: {n_species}")
print(f"Shannon index (genus): {shannon_genus:.3f}")
print(f"Shannon index (species): {shannon_species:.3f}")

# Evenness
evenness_genus = shannon_genus / __import__("math").log(n_genera) if n_genera > 1 else 0
evenness_species = shannon_species / __import__("math").log(n_species) if n_species > 1 else 0
print(f"Evenness (genus): {evenness_genus:.3f}")
print(f"Evenness (species): {evenness_species:.3f}")

# Dominance
dominant_genus = merged["Genus"].value_counts().iloc[0]
dominant_genus_name = merged["Genus"].value_counts().index[0]
dominance = dominant_genus / len(merged)
print(f"Dominant genus: {dominant_genus_name} ({dominant_genus}, {dominance:.1%})")

# Save diversity summary
diversity = pd.DataFrame([
    {"Metric": "Total_KJUs", "Value": len(merged)},
    {"Metric": "Unique_genera", "Value": n_genera},
    {"Metric": "Unique_species", "Value": n_species},
    {"Metric": "Shannon_index_genus", "Value": round(shannon_genus, 3)},
    {"Metric": "Shannon_index_species", "Value": round(shannon_species, 3)},
    {"Metric": "Evenness_genus", "Value": round(evenness_genus, 3)},
    {"Metric": "Evenness_species", "Value": round(evenness_species, 3)},
    {"Metric": "Dominant_genus", "Value": dominant_genus_name},
    {"Metric": "Dominant_genus_count", "Value": dominant_genus},
    {"Metric": "Dominance_percent", "Value": round(dominance * 100, 1)},
])
diversity.to_csv(OUT / "STEP107_DIVERSITY_SUMMARY.csv", index=False)

# Singletons (genera with only 1 species)
singletons = genus_map[genus_map["Count"] == 1]
print(f"Singleton genera (1 species): {len(singletons)}")

# Checksums
for f in sorted(OUT.glob("STEP107_*")):
    if f.is_file():
        with open(OUT / "STEP107_CHECKSUMS.txt", "a") as fh:
            fh.write(f"{sha256(f)}  {f.name}\n")

with open(OUT / "STEP107_REPORT.txt", "w") as fh:
    fh.write("STEP 107: PHYLOGENETIC + TAXONOMIC DIVERSITY\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Input tree: {TREE}\n")
    fh.write(f"Input metadata: {META}\n\n")
    fh.write("TAXONOMIC DIVERSITY:\n")
    fh.write(f"  Total KJUs: {len(merged)}\n")
    fh.write(f"  Unique genera: {n_genera}\n")
    fh.write(f"  Unique species: {n_species}\n\n")
    fh.write("DIVERSITY INDICES:\n")
    fh.write(f"  Shannon index (genus): {shannon_genus:.3f}\n")
    fh.write(f"  Shannon index (species): {shannon_species:.3f}\n")
    fh.write(f"  Evenness (genus): {evenness_genus:.3f}\n")
    fh.write(f"  Evenness (species): {evenness_species:.3f}\n\n")
    fh.write("DOMINANCE:\n")
    fh.write(f"  Dominant genus: {dominant_genus_name} ({dominant_genus}, {dominance:.1%})\n")
    fh.write(f"  Singleton genera: {len(singletons)}\n\n")
    fh.write("GENUS CLUSTERS:\n")
    for _, row in genus_map.iterrows():
        fh.write(f"  {row['Genus']:25s} {row['Count']:3d} species\n")

print("\nOUTPUTS:")
for f in sorted(OUT.glob("STEP107_*")):
    print(f"  {f.name}")
print("\nSTEP 107 COMPLETE")
