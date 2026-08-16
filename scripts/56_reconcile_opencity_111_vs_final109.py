#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

# ============================================================
# INPUTS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

OPENCITY = BASE / "731b3fb2-c763-4491-9621-e04f9e1a9740.csv"

FINAL_METADATA = (
    BASE
    / "07_Master_Metadata"
    / "FINAL_KJU_16S_DATASET_VALIDATED_V1.csv"
)

OUTDIR = BASE / "github_results" / "reconciliation_111_vs_109"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD DATA
# ============================================================

print("=" * 70)
print("DODDAGUBBI 111 vs FINAL 109 RECONCILIATION")
print("=" * 70)

print("\nLoading OpenCity dataset...")
open_df = pd.read_csv(OPENCITY)

print("Loading final KJU metadata...")
final_df = pd.read_csv(FINAL_METADATA)

print("\nOpenCity shape:", open_df.shape)
print("Final metadata shape:", final_df.shape)

# ============================================================
# CLEAN COLUMN NAMES
# ============================================================

open_df.columns = open_df.columns.str.strip()
final_df.columns = final_df.columns.str.strip()

print("\nOpenCity columns:")
print(open_df.columns.tolist())

print("\nFinal metadata columns:")
print(final_df.columns.tolist())

# ============================================================
# IDENTIFY KJU COLUMN
# ============================================================

def find_column(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

open_kju_col = find_column(
    open_df,
    ["Strain ID", "KJU_ID", "KJU ID", "strain_id", "StrainID"]
)

final_kju_col = find_column(
    final_df,
    ["Strain ID", "KJU_ID", "KJU ID", "strain_id", "StrainID", "KJU", "kju_id"]
)

if open_kju_col is None:
    raise ValueError("Could not identify KJU column in OpenCity dataset.")

if final_kju_col is None:
    raise ValueError("Could not identify KJU column in final metadata.")

print("\nOpenCity KJU column :", open_kju_col)
print("Final KJU column    :", final_kju_col)

# ============================================================
# STANDARDIZE KJU IDS
# ============================================================

open_df["_KJU"] = (
    open_df[open_kju_col]
    .astype(str)
    .str.strip()
    .str.upper()
)

final_df["_KJU"] = (
    final_df[final_kju_col]
    .astype(str)
    .str.strip()
    .str.upper()
)

# ============================================================
# BASIC COUNTS
# ============================================================

open_kjus = set(open_df["_KJU"].dropna())
final_kjus = set(final_df["_KJU"].dropna())

both = sorted(open_kjus & final_kjus)
open_only = sorted(open_kjus - final_kjus)
final_only = sorted(final_kjus - open_kjus)

print("\n" + "=" * 70)
print("KJU INVENTORY")
print("=" * 70)

print("OpenCity unique KJU IDs :", len(open_kjus))
print("Final unique KJU IDs    :", len(final_kjus))
print("Present in both         :", len(both))
print("OpenCity only           :", len(open_only))
print("Final only              :", len(final_only))

# ============================================================
# SAVE KJU SET DIFFERENCES
# ============================================================

pd.DataFrame({"KJU_ID": both}).to_csv(
    OUTDIR / "KJU_PRESENT_IN_BOTH.csv",
    index=False
)

pd.DataFrame({"KJU_ID": open_only}).to_csv(
    OUTDIR / "KJU_IN_OPENCITY_NOT_FINAL.csv",
    index=False
)

pd.DataFrame({"KJU_ID": final_only}).to_csv(
    OUTDIR / "KJU_IN_FINAL_NOT_OPENCITY.csv",
    index=False
)

# ============================================================
# MERGE ORIGINAL AND FINAL METADATA
# ============================================================

open_small = open_df.copy()
final_small = final_df.copy()

# Prefix columns to prevent confusion
open_small = open_small.rename(
    columns={
        c: f"OpenCity_{c}"
        for c in open_small.columns
        if c != "_KJU"
    }
)

final_small = final_small.rename(
    columns={
        c: f"Final_{c}"
        for c in final_small.columns
        if c != "_KJU"
    }
)

merged = pd.merge(
    open_small,
    final_small,
    on="_KJU",
    how="outer",
    indicator=True
)

# ============================================================
# ACCESSION COMPARISON
# ============================================================

open_acc_col = "OpenCity_GenBank Acc. No."

final_acc_candidates = [
    "Final_GenBank Acc. No.",
    "Final_GenBank_Accession",
    "Final_Accession",
    "Final_ACCESSION",
    "Final_ACCESSION_VERSION",
    "Final_Accession_Version"
]

final_acc_col = None

for c in final_acc_candidates:
    if c in merged.columns:
        final_acc_col = c
        break

print("\nOpenCity accession column:", open_acc_col)
print("Final accession column   :", final_acc_col)

if open_acc_col in merged.columns and final_acc_col:

    merged["ACCESSION_MATCH"] = (
        merged[open_acc_col]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
        ==
        merged[final_acc_col]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

else:
    merged["ACCESSION_MATCH"] = "NOT_COMPARABLE"

# ============================================================
# SPECIES COMPARISON
# ============================================================

open_species_col = "OpenCity_Species"

final_species_candidates = [
    "Final_Species",
    "Final_species",
    "Final_Organism",
    "Final_Taxonomy",
    "Final_Taxon"
]

final_species_col = None

for c in final_species_candidates:
    if c in merged.columns:
        final_species_col = c
        break

print("OpenCity species column :", open_species_col)
print("Final species column    :", final_species_col)

if open_species_col in merged.columns and final_species_col:

    merged["SPECIES_EXACT_MATCH"] = (
        merged[open_species_col]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        ==
        merged[final_species_col]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )

else:
    merged["SPECIES_EXACT_MATCH"] = "NOT_COMPARABLE"

# ============================================================
# STATUS CLASSIFICATION
# ============================================================

def classify(row):

    if row["_merge"] == "left_only":
        return "OPENCITY_ONLY"

    if row["_merge"] == "right_only":
        return "FINAL_ONLY"

    accession = row["ACCESSION_MATCH"]

    if accession is False:
        return "ACCESSION_DIFFERENCE"

    species = row["SPECIES_EXACT_MATCH"]

    if species is False:
        return "SPECIES_DIFFERENCE"

    return "MATCHED"

merged["RECONCILIATION_STATUS"] = merged.apply(
    classify,
    axis=1
)

# ============================================================
# SAVE COMPLETE RECONCILIATION
# ============================================================

merged.to_csv(
    OUTDIR / "OPENCITY_111_VS_FINAL_109_COMPLETE_RECONCILIATION.csv",
    index=False
)

# ============================================================
# SAVE PROBLEMATIC RECORDS
# ============================================================

problems = merged[
    merged["RECONCILIATION_STATUS"] != "MATCHED"
].copy()

problems.to_csv(
    OUTDIR / "OPENCITY_111_VS_FINAL_109_REVIEW_REQUIRED.csv",
    index=False
)

# ============================================================
# SUMMARY
# ============================================================

summary = (
    merged["RECONCILIATION_STATUS"]
    .value_counts(dropna=False)
    .rename_axis("STATUS")
    .reset_index(name="COUNT")
)

summary.to_csv(
    OUTDIR / "OPENCITY_111_VS_FINAL_109_RECONCILIATION_SUMMARY.csv",
    index=False
)

# ============================================================
# REPORT
# ============================================================

report = []

report.append("DODDAGUBBI OPENCITY 111 vs FINAL 109 RECONCILIATION")
report.append("=" * 70)
report.append("")

report.append(f"Original OpenCity records : {len(open_df)}")
report.append(f"Original OpenCity unique KJU IDs : {len(open_kjus)}")
report.append("")

report.append(f"Final metadata records : {len(final_df)}")
report.append(f"Final metadata unique KJU IDs : {len(final_kjus)}")
report.append("")

report.append(f"KJU IDs present in both : {len(both)}")
report.append(f"KJU IDs only in OpenCity : {len(open_only)}")
report.append(f"KJU IDs only in final metadata : {len(final_only)}")
report.append("")

report.append("RECONCILIATION STATUS")
report.append("-" * 70)

for _, row in summary.iterrows():
    report.append(
        f"{row['STATUS']}: {row['COUNT']}"
    )

report.append("")
report.append("OpenCity-only KJU IDs:")
report.append(", ".join(open_only) if open_only else "None")

report.append("")
report.append("Final-only KJU IDs:")
report.append(", ".join(final_only) if final_only else "None")

report.append("")
report.append("IMPORTANT:")
report.append(
    "This report compares the original OpenCity inventory against "
    "the final KJU metadata. It does not by itself prove sequence "
    "identity or taxonomic correctness."
)

with open(
    OUTDIR / "OPENCITY_111_VS_FINAL_109_RECONCILIATION_REPORT.txt",
    "w",
    encoding="utf-8"
) as f:
    f.write("\n".join(report))

# ============================================================
# PRINT RESULTS
# ============================================================

print("\n" + "=" * 70)
print("RECONCILIATION RESULTS")
print("=" * 70)

print(summary.to_string(index=False))

print("\nOpenCity-only KJU IDs:")
print(open_only if open_only else "None")

print("\nFinal-only KJU IDs:")
print(final_only if final_only else "None")

print("\nOutput directory:")
print(OUTDIR)

print("\nGenerated files:")

for f in sorted(OUTDIR.iterdir()):
    print(" -", f.name)

print("\nDONE.")
