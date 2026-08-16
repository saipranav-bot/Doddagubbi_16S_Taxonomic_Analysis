#!/usr/bin/env python3
"""
STEP 102: Build Analysis-Ready Metadata

Merges the V2 frozen master with genbank sample metadata for the 111
canonical KJUs. Explicitly preserves the 22 KJUs with missing metadata.

Downstream of Step 101. Does NOT modify V2 or V1 freezes.
"""

import hashlib
from pathlib import Path

import pandas as pd

# ──────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────
BASE = Path.home() / "Doddagubbi_Master"
RESULTS = BASE / "github_results/reconciliation_111_vs_109/final_integrated_master"

V2_CSV = RESULTS / "STEP98_FINAL_FREEZE_V2/FINAL_KJU_MASTER_FREEZE_V2/FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv"
V2_FASTA = RESULTS / "STEP98_FINAL_FREEZE_V2/FINAL_KJU_MASTER_FREEZE_V2/FINAL_CANONICAL_KJU_16S_V2.fasta"
GENBANK_META = BASE / "07_Master_Metadata/direct_genbank_kju_metadata.csv"

OUT = BASE / "downstream_analysis/STEP102_ANALYSIS_METADATA"
OUT.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────────────────────
def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calc_gc(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100, 2) if seq else 0.0


print("STEP 102: BUILD ANALYSIS-READY METADATA")
print("=" * 70)
print()

# ──────────────────────────────────────────────────────────
# LOAD DATA
# ──────────────────────────────────────────────────────────
v2 = pd.read_csv(V2_CSV, dtype=str)
genbank = pd.read_csv(GENBANK_META, dtype=str)

print(f"V2 master: {len(v2)} rows")
print(f"Genbank metadata: {len(genbank)} rows")
print()

# ──────────────────────────────────────────────────────────
# FILTER CANONICAL KJUs
# ──────────────────────────────────────────────────────────
canon = v2[v2["Final_status"] == "ADJUDICATED_CANONICAL"].copy()
print(f"Canonical KJUs: {len(canon)}")

# ──────────────────────────────────────────────────────────
# MERGE WITH GENBANK METADATA
# ──────────────────────────────────────────────────────────
# Genbank uses 'kju_id' column; V2 uses 'KJU_ID'
genbank_renamed = genbank.rename(columns={"kju_id": "KJU_ID"})

# Merge: left join to keep all 111 canonical KJUs
merged = canon.merge(genbank_renamed, on="KJU_ID", how="left", suffixes=("", "_genbank"))

# Explicitly track metadata availability
merged["metadata_available"] = merged["accession"].notna() & (merged["accession"] != "")
n_with_meta = merged["metadata_available"].sum()
n_without_meta = (~merged["metadata_available"]).sum()

print(f"With metadata: {n_with_meta}")
print(f"Without metadata: {n_without_meta}")
print()

# ──────────────────────────────────────────────────────────
# CALCULATE GC% FROM CANONICAL SEQUENCES
# ──────────────────────────────────────────────────────────
merged["GC_percent"] = merged["Canonical_sequence"].apply(
    lambda x: calc_gc(str(x)) if pd.notna(x) and str(x).strip() not in {"", "nan", "none"} else None
)

# ──────────────────────────────────────────────────────────
# EXTRACT GENUS FROM SPECIES
# ──────────────────────────────────────────────────────────
merged["Genus"] = merged["Accepted_taxonomy"].apply(
    lambda x: str(x).strip().split()[0] if pd.notna(x) and str(x).strip() else ""
)

merged["Species"] = merged["Accepted_taxonomy"].apply(
    lambda x: " ".join(str(x).strip().split()[:2]) if pd.notna(x) and str(x).strip() else ""
)

# ──────────────────────────────────────────────────────────
# BUILD CLEAN ANALYSIS TABLE
# ──────────────────────────────────────────────────────────
analysis_cols = [
    "KJU_ID",
    "Accepted_accession",
    "Accepted_taxonomy",
    "Genus",
    "Species",
    "Accepted_length",
    "Accepted_sequence_type",
    "Canonical_sequence_length",
    "GC_percent",
    "Canonical_accession",
    "Biological_identity",
    "Final_biological_identity",
    "Confidence",
    "Adjudication_type",
    "Adjudication_reason",
    "Selection_status",
    "Sequence_type_validated",
    "BLAST_taxonomy_agreement",
    "In_phylogenetic_tree",
    "metadata_available",
    # Genbank sample metadata (if available)
    "isolation_source",
    "geo_loc_name",
    "collection_date",
    "collected_by",
    "strain",
]

# Only include columns that exist
analysis_cols = [c for c in analysis_cols if c in merged.columns]

analysis = merged[analysis_cols].copy()

# Clean up NaN strings
for col in analysis.columns:
    analysis[col] = analysis[col].apply(
        lambda x: "" if pd.isna(x) or str(x).strip() in {"nan", "none", "None"} else str(x).strip()
    )

print(f"Analysis table: {len(analysis)} rows × {len(analysis.columns)} columns")
print()

# ──────────────────────────────────────────────────────────
# SUMMARY STATISTICS
# ──────────────────────────────────────────────────────────
print("TAXONOMIC SUMMARY")
print("-" * 70)
genus_counts = analysis["Genus"].value_counts()
print(f"Unique genera: {analysis['Genus'].nunique()}")
print(f"Unique species: {analysis['Species'].nunique()}")
print()
print("Top 10 genera:")
for genus, count in genus_counts.head(10).items():
    print(f"  {genus}: {count}")
print()

print("SEQUENCE STATISTICS")
print("-" * 70)
gc_vals = analysis["GC_percent"].apply(lambda x: float(x) if x else None).dropna()
len_vals = analysis["Canonical_sequence_length"].apply(lambda x: int(float(x)) if x else None).dropna()
print(f"GC% range: {gc_vals.min():.1f}% – {gc_vals.max():.1f}%")
print(f"GC% mean: {gc_vals.mean():.1f}%")
print(f"Length range: {int(len_vals.min())} – {int(len_vals.max())} bp")
print(f"Length mean: {len_vals.mean():.0f} bp")
print()

print("METADATA AVAILABILITY")
print("-" * 70)
print(f"With sample metadata: {n_with_meta}")
print(f"Without sample metadata: {n_without_meta}")
if n_without_meta > 0:
    missing_kjus = analysis[analysis["metadata_available"] == "False"]["KJU_ID"].tolist()
    print(f"Missing KJUs: {missing_kjus[:10]}{'...' if len(missing_kjus) > 10 else ''}")
print()

# ──────────────────────────────────────────────────────────
# KJUs WITHOUT METADATA — EXPLICIT LIST
# ──────────────────────────────────────────────────────────
if n_without_meta > 0:
    print("KJUs WITHOUT SAMPLE METADATA (explicit list)")
    print("-" * 70)
    missing = analysis[analysis["metadata_available"] == "False"]
    for _, row in missing.iterrows():
        print(f"  {row['KJU_ID']}: {row['Accepted_taxonomy']}")
    print()

# ──────────────────────────────────────────────────────────
# SAVE OUTPUTS
# ──────────────────────────────────────────────────────────
analysis.to_csv(OUT / "STEP102_ANALYSIS_METADATA.csv", index=False)

# Also save a version with only the 89 KJUs that have full metadata
analysis_complete = analysis[analysis["metadata_available"] == "True"].copy()
analysis_complete.to_csv(OUT / "STEP102_ANALYSIS_METADATA_COMPLETE.csv", index=False)

# Save missing metadata report
missing_report = analysis[analysis["metadata_available"] == "False"].copy()
missing_report.to_csv(OUT / "STEP102_MISSING_METADATA_REPORT.csv", index=False)

# ──────────────────────────────────────────────────────────
# CHECKSUMS
# ──────────────────────────────────────────────────────────
sha256_map = {}
for f in sorted(OUT.glob("STEP102_*")):
    if f.is_file():
        sha = sha256(f)
        sha256_map[f.name] = sha

with open(OUT / "STEP102_CHECKSUMS.txt", "w") as fh:
    fh.write("STEP 102 CHECKSUMS\n")
    fh.write("=" * 60 + "\n\n")
    for fname, sha in sha256_map.items():
        fh.write(f"{sha}  {fname}\n")
    fh.write(f"\nV2 master SHA256: {sha256(V2_CSV)}\n")

# ──────────────────────────────────────────────────────────
# QC REPORT
# ──────────────────────────────────────────────────────────
with open(OUT / "STEP102_REPORT.txt", "w") as fh:
    fh.write("STEP 102: BUILD ANALYSIS-READY METADATA\n")
    fh.write("=" * 70 + "\n\n")

    fh.write("INPUTS\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  V2 master: {V2_CSV}\n")
    fh.write(f"  Genbank metadata: {GENBANK_META}\n")
    fh.write(f"  V2 SHA256: {sha256(V2_CSV)}\n\n")

    fh.write("RESULTS\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  Canonical KJUs: {len(canon)}\n")
    fh.write(f"  With sample metadata: {n_with_meta}\n")
    fh.write(f"  Without sample metadata: {n_without_meta}\n")
    fh.write(f"  Analysis columns: {len(analysis.columns)}\n\n")

    fh.write("TAXONOMIC COMPOSITION\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  Unique genera: {analysis['Genus'].nunique()}\n")
    fh.write(f"  Unique species: {analysis['Species'].nunique()}\n\n")

    fh.write("Top 10 genera:\n")
    for genus, count in genus_counts.head(10).items():
        fh.write(f"    {genus}: {count}\n")
    fh.write("\n")

    fh.write("SEQUENCE STATISTICS\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  GC% range: {gc_vals.min():.1f}% – {gc_vals.max():.1f}%\n")
    fh.write(f"  GC% mean: {gc_vals.mean():.1f}%\n")
    fh.write(f"  Length range: {int(len_vals.min())} – {int(len_vals.max())} bp\n")
    fh.write(f"  Length mean: {len_vals.mean():.0f} bp\n\n")

    if n_without_meta > 0:
        fh.write("KJUs WITHOUT SAMPLE METADATA\n")
        fh.write("-" * 70 + "\n")
        for _, row in missing.iterrows():
            fh.write(f"    {row['KJU_ID']}: {row['Accepted_taxonomy']}\n")
        fh.write("\n")

    fh.write("OUTPUTS\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  STEP102_ANALYSIS_METADATA.csv (all 111)\n")
    fh.write(f"  STEP102_ANALYSIS_METADATA_COMPLETE.csv ({n_with_meta} with metadata)\n")
    fh.write(f"  STEP102_MISSING_METADATA_REPORT.csv ({n_without_meta} without metadata)\n")
    fh.write(f"  STEP102_REPORT.txt\n")
    fh.write(f"  STEP102_CHECKSUMS.txt\n\n")

    fh.write("DATA LIMITATIONS\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  {n_without_meta} KJUs have no sample metadata\n")
    fh.write(f"  No abundance/read-count data available\n")
    fh.write(f"  No multiple-sample design\n")
    fh.write(f"  All sequences from Doddagubbi Lake area\n")

# ──────────────────────────────────────────────────────────
# FINAL SUMMARY
# ──────────────────────────────────────────────────────────
print("OUTPUTS:")
for f in sorted(OUT.glob("STEP102_*")):
    print(f"  {f.name}")
print()

print("ANALYSIS TABLE COLUMNS:")
for i, col in enumerate(analysis.columns, 1):
    print(f"  {i:02d}. {col}")
print()

print("=" * 70)
print(f"  Total canonical: {len(analysis)}")
print(f"  With metadata: {n_with_meta}")
print(f"  Without metadata: {n_without_meta}")
print(f"  Genera: {analysis['Genus'].nunique()}")
print(f"  Species: {analysis['Species'].nunique()}")
print(f"  GC%: {gc_vals.mean():.1f}% mean")
print(f"  Length: {len_vals.mean():.0f} bp mean")
print()
print("STEP 102 COMPLETE")
