from pathlib import Path
from Bio import AlignIO
from collections import Counter
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
ALIGN = MASTER / "07_Master_Metadata/FINAL_109_ALIGNMENT_V1/FINAL_109_KJU_16S_MAFFT.fasta"

OUT = MASTER / "07_Master_Metadata/FINAL_109_ALIGNMENT_V1/ALIGNMENT_QC_V1.txt"

aln = AlignIO.read(ALIGN, "fasta")

n = len(aln)
L = aln.get_alignment_length()

print("=" * 90)
print(" DODDAGGUBBI - MAFFT ALIGNMENT QC")
print("=" * 90)

print("Sequences:", n)
print("Alignment length:", L)

# ------------------------------------------------------------
# Column statistics
# ------------------------------------------------------------

gap_chars = set("-.")

column_gap_percent = []
variable_sites = 0
parsimony_informative = 0
constant_sites = 0

for i in range(L):

    column = [
        str(record.seq[i]).upper()
        for record in aln
    ]

    nongap = [
        x for x in column
        if x not in gap_chars
    ]

    if n == 0:
        gap_pct = 100
    else:
        gap_pct = (
            (n - len(nongap)) / n * 100
        )

    column_gap_percent.append(gap_pct)

    states = Counter(nongap)

    if len(states) == 1:
        constant_sites += 1

    elif len(states) > 1:
        variable_sites += 1

        # At least two states occur in >=2 sequences
        informative_states = sum(
            1 for count in states.values()
            if count >= 2
        )

        if informative_states >= 2:
            parsimony_informative += 1


# ------------------------------------------------------------
# Sequence gap statistics
# ------------------------------------------------------------

sequence_stats = []

for record in aln:

    seq = str(record.seq).upper()

    gaps = sum(
        1 for x in seq
        if x in gap_chars
    )

    nongap = L - gaps

    gap_pct = (
        gaps / L * 100
    )

    sequence_stats.append({
        "id": record.id,
        "alignment_length": L,
        "gap_count": gaps,
        "nongap_count": nongap,
        "gap_percent": gap_pct,
    })

seq_df = pd.DataFrame(sequence_stats)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

max_col_gap = max(column_gap_percent)
mean_col_gap = sum(column_gap_percent) / L

high_gap_50 = sum(
    1 for x in column_gap_percent
    if x >= 50
)

high_gap_80 = sum(
    1 for x in column_gap_percent
    if x >= 80
)

high_gap_95 = sum(
    1 for x in column_gap_percent
    if x >= 95
)

print("\nSITE SUMMARY")
print("-" * 90)

print("Constant sites:", constant_sites)
print("Variable sites:", variable_sites)
print("Parsimony-informative sites:", parsimony_informative)

print("\nGAP SUMMARY")
print("-" * 90)

print("Mean column gap %:", round(mean_col_gap, 2))
print("Maximum column gap %:", round(max_col_gap, 2))
print("Columns >=50% gaps:", high_gap_50)
print("Columns >=80% gaps:", high_gap_80)
print("Columns >=95% gaps:", high_gap_95)

print("\nSEQUENCE GAP SUMMARY")
print("-" * 90)

print(
    "Minimum sequence gap %:",
    round(seq_df["gap_percent"].min(), 2)
)

print(
    "Maximum sequence gap %:",
    round(seq_df["gap_percent"].max(), 2)
)

print(
    "Mean sequence gap %:",
    round(seq_df["gap_percent"].mean(), 2)
)

print("\nMost gappy sequences:")
print(
    seq_df.sort_values(
        "gap_percent",
        ascending=False
    ).head(10).to_string(index=False)
)

# ------------------------------------------------------------
# Write report
# ------------------------------------------------------------

lines = []

lines.append("DODDAGGUBBI - MAFFT ALIGNMENT QC")
lines.append("=" * 90)
lines.append("")
lines.append(f"Sequences: {n}")
lines.append(f"Alignment length: {L}")
lines.append("")
lines.append(f"Constant sites: {constant_sites}")
lines.append(f"Variable sites: {variable_sites}")
lines.append(
    f"Parsimony-informative sites: {parsimony_informative}"
)
lines.append("")
lines.append(
    f"Mean column gap %: {mean_col_gap:.2f}"
)
lines.append(
    f"Maximum column gap %: {max_col_gap:.2f}"
)
lines.append(
    f"Columns >=50% gaps: {high_gap_50}"
)
lines.append(
    f"Columns >=80% gaps: {high_gap_80}"
)
lines.append(
    f"Columns >=95% gaps: {high_gap_95}"
)
lines.append("")
lines.append(
    f"Minimum sequence gap %: "
    f"{seq_df['gap_percent'].min():.2f}"
)
lines.append(
    f"Maximum sequence gap %: "
    f"{seq_df['gap_percent'].max():.2f}"
)
lines.append(
    f"Mean sequence gap %: "
    f"{seq_df['gap_percent'].mean():.2f}"
)
lines.append("")
lines.append("No source sequences were modified.")

OUT.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8"
)

print("\nReport:")
print(OUT)

print("\nALIGNMENT QC COMPLETE")
