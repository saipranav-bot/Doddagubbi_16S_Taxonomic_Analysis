#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from pathlib import Path

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

# ------------------------------------------------------------
# LOAD
# ------------------------------------------------------------

old_records = list(SeqIO.parse(EXISTING, "fasta"))
new_records = list(SeqIO.parse(NEW, "fasta"))

print("=" * 80)
print("NEW SEQUENCES vs EXISTING 109 KJU SEQUENCES")
print("=" * 80)

print("\nExisting sequences:", len(old_records))
print("New sequences     :", len(new_records))

# ------------------------------------------------------------
# ALIGNER
# ------------------------------------------------------------

aligner = PairwiseAligner()

aligner.mode = "global"

aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.1

results = []

for new_record in new_records:

    new_seq = str(new_record.seq).upper()

    print("\n" + "=" * 80)
    print("QUERY:", new_record.id)
    print("Length:", len(new_seq))
    print("=" * 80)

    comparisons = []

    for old_record in old_records:

        old_seq = str(old_record.seq).upper()

        score = aligner.score(new_seq, old_seq)

        max_len = max(len(new_seq), len(old_seq))

        identity_estimate = (
            score / max_len * 100
            if max_len
            else 0
        )

        comparisons.append(
            (
                old_record.id,
                len(old_seq),
                score,
                identity_estimate
            )
        )

    comparisons.sort(
        key=lambda x: x[3],
        reverse=True
    )

    print("\nTop 10 matches:")

    for rank, item in enumerate(comparisons[:10], 1):

        kju, length, score, identity = item

        print(
            f"{rank:2d}. "
            f"{kju:12s} "
            f"length={length:4d} "
            f"score={score:8.2f} "
            f"identity_est={identity:6.2f}%"
        )

        results.append({
            "query": new_record.id,
            "existing_kju": kju,
            "existing_length": length,
            "query_length": len(new_seq),
            "alignment_score": score,
            "identity_estimate_percent": identity
        })

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------

import pandas as pd

df = pd.DataFrame(results)

outfile = OUTDIR / "NEW_VS_109_SEQUENCE_COMPARISON.csv"

df.to_csv(outfile, index=False)

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)

print("\nResults:")
print(outfile)

print("\nHighest similarity for each new sequence:")

best = (
    df.sort_values(
        ["query", "identity_estimate_percent"],
        ascending=[True, False]
    )
    .groupby("query")
    .head(1)
)

print(
    best[
        [
            "query",
            "existing_kju",
            "query_length",
            "existing_length",
            "alignment_score",
            "identity_estimate_percent"
        ]
    ].to_string(index=False)
)
