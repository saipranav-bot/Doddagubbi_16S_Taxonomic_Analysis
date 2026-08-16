#!/usr/bin/env python3
"""STEP 103: Taxonomic Composition Analysis"""

import hashlib
from pathlib import Path
import pandas as pd

BASE = Path.home() / "Doddagubbi_Master"
OUT = BASE / "downstream_analysis/STEP103_TAXONOMIC_COMPOSITION"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(BASE / "downstream_analysis/STEP102_ANALYSIS_METADATA/STEP102_ANALYSIS_METADATA.csv", dtype=str)

def sha256(p): return hashlib.sha256(p.read_bytes()).hexdigest()

print("STEP 103: TAXONOMIC COMPOSITION")
print("=" * 60)

# Genus counts
genus_counts = meta["Genus"].value_counts().reset_index()
genus_counts.columns = ["Genus", "Count"]
genus_counts["Percent"] = (genus_counts["Count"] / genus_counts["Count"].sum() * 100).round(2)
genus_counts.to_csv(OUT / "STEP103_GENUS_COUNTS.csv", index=False)

# Species counts
species_counts = meta["Species"].value_counts().reset_index()
species_counts.columns = ["Species", "Count"]
species_counts["Percent"] = (species_counts["Count"] / species_counts["Count"].sum() * 100).round(2)
species_counts.to_csv(OUT / "STEP103_SPECIES_COUNTS.csv", index=False)

# Full taxonomy table
tax_table = meta[["KJU_ID", "Accepted_accession", "Accepted_taxonomy", "Genus", "Species"]].copy()
tax_table.to_csv(OUT / "STEP103_FULL_TAXONOMY.csv", index=False)

# Summary
print(f"Unique genera: {meta['Genus'].nunique()}")
print(f"Unique species: {meta['Species'].nunique()}")
print()
print("GENUS DISTRIBUTION:")
for _, row in genus_counts.iterrows():
    bar = "█" * int(row["Percent"] / 2)
    print(f"  {row['Genus']:25s} {row['Count']:3d} ({row['Percent']:5.1f}%) {bar}")
print()
print("SPECIES DISTRIBUTION:")
for _, row in species_counts.iterrows():
    print(f"  {row['Species']:40s} {row['Count']:3d} ({row['Percent']:5.1f}%)")

# Checksums + report
for f in sorted(OUT.glob("STEP103_*")):
    if f.is_file():
        with open(OUT / "STEP103_CHECKSUMS.txt", "a") as fh:
            fh.write(f"{sha256(f)}  {f.name}\n")

with open(OUT / "STEP103_REPORT.txt", "w") as fh:
    fh.write("STEP 103: TAXONOMIC COMPOSITION\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Total canonical KJUs: {len(meta)}\n")
    fh.write(f"Unique genera: {meta['Genus'].nunique()}\n")
    fh.write(f"Unique species: {meta['Species'].nunique()}\n\n")
    fh.write("GENUS COUNTS:\n")
    for _, row in genus_counts.iterrows():
        fh.write(f"  {row['Genus']:25s} {row['Count']:3d} ({row['Percent']:5.1f}%)\n")
    fh.write("\nSPECIES COUNTS:\n")
    for _, row in species_counts.iterrows():
        fh.write(f"  {row['Species']:40s} {row['Count']:3d} ({row['Percent']:5.1f}%)\n")

print("\nOUTPUTS:")
for f in sorted(OUT.glob("STEP103_*")):
    print(f"  {f.name}")
print("\nSTEP 103 COMPLETE")
