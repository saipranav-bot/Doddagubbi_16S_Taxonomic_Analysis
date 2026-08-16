#!/usr/bin/env python3

from pathlib import Path
from Bio import Entrez, SeqIO
import pandas as pd
import time

# ============================================================
# SETTINGS
# ============================================================

Entrez.email = "saipranavkar@gmail.com"

BASE = Path.home() / "Doddagubbi_Master"

INPUT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
    / "KJU_TAXONOMIC_IDENTITY_RESULTS.csv"
)

OUTDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

OUTFILE = OUTDIR / "BEST_REFERENCE_TAXONOMY.csv"

# ============================================================
# LOAD RESULTS
# ============================================================

df = pd.read_csv(INPUT)

# Get best reference for each KJU
best = (
    df.sort_values(
        ["query", "identity_percent", "query_coverage_percent"],
        ascending=[True, False, False]
    )
    .groupby("query")
    .head(5)
    .copy()
)

accessions = best["reference"].astype(str).tolist()

print("=" * 80)
print("STEP 64 — IDENTIFY TAXONOMY OF BEST NCBI REFERENCES")
print("=" * 80)

print("\nReferences to inspect:")
for acc in accessions:
    print(" -", acc)

# ============================================================
# QUERY NCBI
# ============================================================

taxonomy = {}

for acc in accessions:

    print("\nChecking:", acc)

    try:

        handle = Entrez.efetch(
            db="nuccore",
            id=acc,
            rettype="gb",
            retmode="text"
        )

        text = handle.read()
        handle.close()

        record = SeqIO.read(
            __import__("io").StringIO(text),
            "genbank"
        )

        organism = record.annotations.get(
            "organism",
            "UNKNOWN"
        )

        description = record.description

        taxonomy[acc] = {
            "organism": organism,
            "description": description
        }

        print("Organism   :", organism)
        print("Description:", description)

        time.sleep(0.34)

    except Exception as e:

        print("ERROR:", e)

        taxonomy[acc] = {
            "organism": "ERROR",
            "description": str(e)
        }

# ============================================================
# ADD TAXONOMY
# ============================================================

best["reference_organism"] = best["reference"].map(
    lambda x: taxonomy.get(
        x,
        {}
    ).get(
        "organism",
        "UNKNOWN"
    )
)

best["reference_description"] = best["reference"].map(
    lambda x: taxonomy.get(
        x,
        {}
    ).get(
        "description",
        "UNKNOWN"
    )
)

# ============================================================
# SAVE
# ============================================================

best.to_csv(
    OUTFILE,
    index=False
)

# ============================================================
# PRINT SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("REFERENCE TAXONOMY SUMMARY")
print("=" * 80)

for query in ["KJU0224", "KJU0225"]:

    print("\n" + "-" * 80)
    print(query)
    print("-" * 80)

    subset = best[
        best["query"] == query
    ]

    for _, row in subset.iterrows():

        print(
            f"\nReference : {row['reference']}"
        )

        print(
            f"Organism  : {row['reference_organism']}"
        )

        print(
            f"Identity  : "
            f"{row['identity_percent']:.2f}%"
        )

        print(
            f"Coverage  : "
            f"{row['query_coverage_percent']:.2f}%"
        )

# ============================================================
# DONE
# ============================================================

print("\n" + "=" * 80)
print("STEP 64 COMPLETE")
print("=" * 80)

print("\nSaved:")
print(OUTFILE)

