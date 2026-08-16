#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
)

INPUT = ROOT / "OPENCITY_111_VS_FINAL_109_COMPLETE_RECONCILIATION.csv"

OUT = ROOT / "biological_reconciliation"

OUT.mkdir(parents=True, exist_ok=True)

OUTPUT = OUT / "SHARED_KJU_BIOLOGICAL_RECONCILIATION.csv"
SUMMARY = OUT / "SHARED_KJU_BIOLOGICAL_RECONCILIATION_SUMMARY.txt"


# ============================================================
# HELPERS
# ============================================================

def clean(value):
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_accession(value):
    """
    Normalize accession for comparison.

    Example:
        PZ666011
        PZ666011.1

    Both become:
        PZ666011
    """
    value = clean(value)

    if not value:
        return ""

    value = value.upper()

    return re.sub(r"\.\d+$", "", value)


def normalize_species(value):
    """
    Basic species normalization:
    - strip whitespace
    - collapse repeated spaces
    - lowercase
    """
    value = clean(value)

    if not value:
        return ""

    value = re.sub(r"\s+", " ", value)

    return value.lower()


# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INPUT)

# Only records where KJU exists in BOTH datasets
shared = df[df["_merge"] == "both"].copy()

print("=" * 90)
print("STEP 71: SHARED KJU BIOLOGICAL RECONCILIATION")
print("=" * 90)

print(f"Total rows: {len(df)}")
print(f"Shared KJU records: {len(shared)}")


# ============================================================
# NORMALIZED ACCESSIONS
# ============================================================

shared["OpenCity_accession_normalized"] = (
    shared["OpenCity_GenBank Acc. No."]
    .apply(normalize_accession)
)

shared["Final_accession_normalized"] = (
    shared["Final_accession"]
    .apply(normalize_accession)
)


# ============================================================
# ACCESSION COMPARISON
# ============================================================

def accession_status(row):

    oc = row["OpenCity_accession_normalized"]
    final = row["Final_accession_normalized"]

    if not oc and not final:
        return "BOTH_MISSING"

    if not oc:
        return "MISSING_IN_OPENCITY"

    if not final:
        return "MISSING_IN_FINAL"

    if oc == final:
        return "ACCESSION_MATCH"

    return "ACCESSION_MISMATCH"


shared["ACCESSION_RECONCILIATION"] = (
    shared.apply(accession_status, axis=1)
)


# ============================================================
# NORMALIZED TAXONOMY
# ============================================================

shared["OpenCity_species_normalized"] = (
    shared["OpenCity_Species"]
    .apply(normalize_species)
)

shared["Final_organism_normalized"] = (
    shared["Final_organism"]
    .apply(normalize_species)
)


# ============================================================
# SPECIES COMPARISON
# ============================================================

def species_status(row):

    oc = row["OpenCity_species_normalized"]
    final = row["Final_organism_normalized"]

    if not oc and not final:
        return "BOTH_MISSING"

    if not oc:
        return "MISSING_IN_OPENCITY"

    if not final:
        return "MISSING_IN_FINAL"

    if oc == final:
        return "SPECIES_EXACT_MATCH"

    # Compare genus + species tokens
    oc_tokens = oc.split()
    final_tokens = final.split()

    if len(oc_tokens) >= 2 and len(final_tokens) >= 2:

        if (
            oc_tokens[-1] == final_tokens[-1]
            and oc_tokens[0] == final_tokens[0]
        ):
            return "SPECIES_FORMAT_VARIATION"

    return "SPECIES_MISMATCH"


shared["SPECIES_RECONCILIATION"] = (
    shared.apply(species_status, axis=1)
)


# ============================================================
# OVERALL BIOLOGICAL STATUS
# ============================================================

def overall_status(row):

    accession = row["ACCESSION_RECONCILIATION"]
    species = row["SPECIES_RECONCILIATION"]

    if (
        accession == "ACCESSION_MATCH"
        and species == "SPECIES_EXACT_MATCH"
    ):
        return "FULL_BIOLOGICAL_MATCH"

    if (
        accession == "ACCESSION_MATCH"
        and species == "SPECIES_FORMAT_VARIATION"
    ):
        return "MATCH_WITH_TAXONOMY_FORMAT_VARIATION"

    if accession == "ACCESSION_MISMATCH":
        return "ACCESSION_CONFLICT"

    if species == "SPECIES_MISMATCH":
        return "SPECIES_CONFLICT"

    return "REQUIRES_REVIEW"


shared["OVERALL_BIOLOGICAL_STATUS"] = (
    shared.apply(overall_status, axis=1)
)


# ============================================================
# SAVE
# ============================================================

shared.to_csv(OUTPUT, index=False)


# ============================================================
# SUMMARY
# ============================================================

with open(SUMMARY, "w") as f:

    f.write("=" * 90 + "\n")
    f.write("DODDAGUBBI SHARED KJU BIOLOGICAL RECONCILIATION\n")
    f.write("=" * 90 + "\n\n")

    f.write(f"Total shared KJU records: {len(shared)}\n\n")

    f.write("ACCESSION RECONCILIATION\n")
    f.write("-" * 90 + "\n")

    for status, count in (
        shared["ACCESSION_RECONCILIATION"]
        .value_counts()
        .items()
    ):
        f.write(f"{status}: {count}\n")

    f.write("\n")

    f.write("SPECIES RECONCILIATION\n")
    f.write("-" * 90 + "\n")

    for status, count in (
        shared["SPECIES_RECONCILIATION"]
        .value_counts()
        .items()
    ):
        f.write(f"{status}: {count}\n")

    f.write("\n")

    f.write("OVERALL BIOLOGICAL STATUS\n")
    f.write("-" * 90 + "\n")

    for status, count in (
        shared["OVERALL_BIOLOGICAL_STATUS"]
        .value_counts()
        .items()
    ):
        f.write(f"{status}: {count}\n")

    f.write("\n")
    f.write("=" * 90 + "\n")
    f.write("STEP 71 COMPLETE\n")
    f.write("=" * 90 + "\n")

print("\nOUTPUT:")
print(OUTPUT)

print("\nSUMMARY:")
print(SUMMARY)

print("\nOVERALL STATUS:")
print(shared["OVERALL_BIOLOGICAL_STATUS"].value_counts())

