#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

ALIGNMENT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
    / "NEW_KJUS_PLUS_NCBI_REFERENCES_MAFFT.fasta"
)

OUTDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "KJU_TAXONOMIC_IDENTITY_RESULTS.csv"

# ============================================================
# LOAD ALIGNMENT
# ============================================================

records = list(SeqIO.parse(ALIGNMENT, "fasta"))

new_records = [
    r for r in records
    if r.id in ["KJU0224", "KJU0225"]
]

reference_records = [
    r for r in records
    if r.id not in ["KJU0224", "KJU0225"]
]

print("=" * 80)
print("STEP 63 — TAXONOMIC IDENTITY VALIDATION")
print("=" * 80)

print(f"\nAligned sequences : {len(records)}")
print(f"New KJU sequences : {len(new_records)}")
print(f"NCBI references   : {len(reference_records)}")

# ============================================================
# PAIRWISE ALIGNMENT SETTINGS
# ============================================================

aligner = PairwiseAligner()

aligner.mode = "global"

aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = 0
aligner.extend_gap_score = 0

results = []

# ============================================================
# CALCULATE IDENTITY
# ============================================================

for query in new_records:

    query_seq = str(query.seq).upper()

    print("\n" + "=" * 80)
    print("QUERY:", query.id)
    print("=" * 80)

    query_results = []

    for ref in reference_records:

        ref_seq = str(ref.seq).upper()

        # ----------------------------------------------------
        # Remove alignment gaps for comparison
        # ----------------------------------------------------

        q = query_seq.replace("-", "")
        r = ref_seq.replace("-", "")

        alignment = aligner.align(q, r)[0]

        aligned_q = alignment[0]
        aligned_r = alignment[1]

        matches = 0
        aligned_positions = 0

        for a, b in zip(aligned_q, aligned_r):

            if a == "-" or b == "-":
                continue

            aligned_positions += 1

            if a == b:
                matches += 1

        if aligned_positions == 0:
            continue

        identity = (
            matches / aligned_positions
        ) * 100

        coverage = (
            aligned_positions / len(q)
        ) * 100

        query_results.append({
            "query": query.id,
            "reference": ref.id,
            "query_length": len(q),
            "reference_length": len(r),
            "matches": matches,
            "aligned_positions": aligned_positions,
            "identity_percent": identity,
            "query_coverage_percent": coverage
        })

    # --------------------------------------------------------
    # SORT BY IDENTITY
    # --------------------------------------------------------

    query_results.sort(
        key=lambda x: (
            x["identity_percent"],
            x["query_coverage_percent"]
        ),
        reverse=True
    )

    results.extend(query_results)

    print("\nTOP 10 REFERENCES")
    print("-" * 80)

    for i, row in enumerate(query_results[:10], 1):

        print(
            f"{i:2d}. "
            f"{row['reference']:<60} "
            f"identity={row['identity_percent']:.2f}% "
            f"coverage={row['query_coverage_percent']:.2f}%"
        )

# ============================================================
# SAVE RESULTS
# ============================================================

df = pd.DataFrame(results)

df.to_csv(
    OUTFILE,
    index=False
)

# ============================================================
# BEST MATCHES
# ============================================================

print("\n" + "=" * 80)
print("BEST TAXONOMIC MATCHES")
print("=" * 80)

for query in ["KJU0224", "KJU0225"]:

    subset = df[df["query"] == query]

    if subset.empty:
        continue

    best = subset.iloc[0]

    print(f"\n{query}")

    print(
        "Best reference :",
        best["reference"]
    )

    print(
        "Identity       :",
        f"{best['identity_percent']:.2f}%"
    )

    print(
        "Coverage       :",
        f"{best['query_coverage_percent']:.2f}%"
    )

    print(
        "Matches        :",
        int(best["matches"])
    )

    print(
        "Aligned sites  :",
        int(best["aligned_positions"])
    )

# ============================================================
# DONE
# ============================================================

print("\n" + "=" * 80)
print("STEP 63 COMPLETE")
print("=" * 80)

print("\nResults saved to:")
print(OUTFILE)

