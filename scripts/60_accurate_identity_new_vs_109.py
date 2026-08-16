#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pathlib import Path
import pandas as pd

BASE = Path.home() / "Doddagubbi_Master"

EXISTING = (
    BASE
    / "07_Master_Metadata"
    / "FINAL_109_ALIGNMENT_V1"
    / "FINAL_109_KJU_16S_MAFFT.fasta"
)

NEW = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "ncbi_check"
    / "OPENCITY_MISSING_KJU0224_KJU0225.fasta"
)

OUTDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "sequence_comparison"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "ACCURATE_NEW_VS_109_IDENTITY.csv"

# ============================================================
# LOAD SEQUENCES
# ============================================================

old_records = list(SeqIO.parse(EXISTING, "fasta"))
new_records = list(SeqIO.parse(NEW, "fasta"))

print("=" * 80)
print("ACCURATE PAIRWISE IDENTITY ANALYSIS")
print("=" * 80)

print(f"\nExisting sequences : {len(old_records)}")
print(f"New sequences      : {len(new_records)}")

# ============================================================
# ALIGNER
# ============================================================

aligner = PairwiseAligner()

# Global alignment
aligner.mode = "global"

# Match / mismatch
aligner.match_score = 1
aligner.mismatch_score = -1

# Gap penalties
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

# ============================================================
# FUNCTION
# ============================================================

def calculate_identity(seq1, seq2):

    alignment = aligner.align(seq1, seq2)[0]

    # Coordinates describe aligned blocks
    coords = alignment.coordinates

    matches = 0
    mismatches = 0
    gaps = 0
    aligned_positions = 0

    for i in range(0, len(coords[0]) - 1):

        start1 = coords[0][i]
        end1   = coords[0][i + 1]

        start2 = coords[1][i]
        end2   = coords[1][i + 1]

        len1 = end1 - start1
        len2 = end2 - start2

        # Both sequences contain bases
        if len1 > 0 and len2 > 0:

            block_len = min(len1, len2)

            s1 = seq1[start1:start1 + block_len]
            s2 = seq2[start2:start2 + block_len]

            for a, b in zip(s1, s2):

                aligned_positions += 1

                if a == b:
                    matches += 1
                else:
                    mismatches += 1

        # Gap in sequence 2
        elif len1 > 0 and len2 == 0:

            gaps += len1
            aligned_positions += len1

        # Gap in sequence 1
        elif len2 > 0 and len1 == 0:

            gaps += len2
            aligned_positions += len2

    identity = (
        matches / aligned_positions * 100
        if aligned_positions
        else 0
    )

    # Query coverage
    query_length = len(seq1)

    covered_query = 0

    for i in range(0, len(coords[0]) - 1):

        start1 = coords[0][i]
        end1 = coords[0][i + 1]

        start2 = coords[1][i]
        end2 = coords[1][i + 1]

        if end1 > start1 and end2 > start2:
            covered_query += min(
                end1 - start1,
                end2 - start2
            )

    coverage = covered_query / query_length * 100

    return {
        "alignment_score": alignment.score,
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "aligned_positions": aligned_positions,
        "identity_percent": identity,
        "query_coverage_percent": coverage,
        "alignment_length": aligned_positions,
    }


# ============================================================
# RUN COMPARISONS
# ============================================================

results = []

for new_record in new_records:

    query_id = new_record.id
    query_seq = str(new_record.seq).upper()

    print("\n" + "=" * 80)
    print("QUERY:", query_id)
    print("Length:", len(query_seq))
    print("=" * 80)

    query_results = []

    for old_record in old_records:

        target_id = old_record.id
        target_seq = str(old_record.seq).upper()

        metrics = calculate_identity(
            query_seq,
            target_seq
        )

        row = {
            "query": query_id,
            "target": target_id,
            "query_length": len(query_seq),
            "target_length": len(target_seq),
            **metrics
        }

        query_results.append(row)

    # Sort by identity
    query_results = sorted(
        query_results,
        key=lambda x: (
            x["identity_percent"],
            x["query_coverage_percent"]
        ),
        reverse=True
    )

    # Keep all results
    results.extend(query_results)

    print("\nTOP 10 BY TRUE PAIRWISE IDENTITY")

    for rank, row in enumerate(query_results[:10], 1):

        print(
            f"{rank:2d}. "
            f"{row['target']} "
            f"identity={row['identity_percent']:.2f}% "
            f"coverage={row['query_coverage_percent']:.2f}% "
            f"matches={row['matches']} "
            f"mismatches={row['mismatches']} "
            f"gaps={row['gaps']}"
        )

# ============================================================
# SAVE ALL RESULTS
# ============================================================

df = pd.DataFrame(results)

df.to_csv(
    OUTFILE,
    index=False
)

print("\n" + "=" * 80)
print("BEST MATCHES")
print("=" * 80)

for query in df["query"].unique():

    sub = df[df["query"] == query].sort_values(
        ["identity_percent", "query_coverage_percent"],
        ascending=False
    )

    best = sub.iloc[0]

    print("\n", query)

    print(
        "Best existing sequence:",
        best["target"]
    )

    print(
        f"Identity: {best['identity_percent']:.2f}%"
    )

    print(
        f"Coverage: {best['query_coverage_percent']:.2f}%"
    )

    print(
        f"Matches: {best['matches']}"
    )

    print(
        f"Mismatches: {best['mismatches']}"
    )

    print(
        f"Gaps: {best['gaps']}"
    )

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)

print("\nFull results:")
print(OUTFILE)

