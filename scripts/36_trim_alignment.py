from pathlib import Path
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter

MASTER = Path.home() / "Doddagubbi_Master"

INPUT = (
    MASTER
    / "07_Master_Metadata/FINAL_109_ALIGNMENT_V1/FINAL_109_KJU_16S_MAFFT.fasta"
)

OUTDIR = (
    MASTER
    / "07_Master_Metadata/FINAL_109_ALIGNMENT_V1"
)

OUTPUT = OUTDIR / "FINAL_109_KJU_16S_MAFFT_TRIMMED_50GAP.fasta"
REPORT = OUTDIR / "ALIGNMENT_TRIMMING_REPORT_V1.txt"

aln = AlignIO.read(INPUT, "fasta")

n = len(aln)
L = aln.get_alignment_length()

gap_chars = set("-.")

keep_columns = []
removed_columns = []

for i in range(L):

    column = [
        str(record.seq[i]).upper()
        for record in aln
    ]

    gaps = sum(
        1 for x in column
        if x in gap_chars
    )

    gap_fraction = gaps / n

    # Remove columns with >=50% gaps
    if gap_fraction < 0.50:
        keep_columns.append(i)
    else:
        removed_columns.append(i)


trimmed_records = []

for record in aln:

    seq = "".join(
        str(record.seq[i])
        for i in keep_columns
    )

    trimmed_records.append(
        record[:0] + seq
    )


trimmed = MultipleSeqAlignment(
    trimmed_records
)

# ------------------------------------------------------------
# Site statistics
# ------------------------------------------------------------

constant = 0
variable = 0
informative = 0

for i in range(trimmed.get_alignment_length()):

    column = [
        str(record.seq[i]).upper()
        for record in trimmed
        if str(record.seq[i]).upper()
        not in gap_chars
    ]

    states = Counter(column)

    if len(states) == 1:
        constant += 1

    elif len(states) > 1:

        variable += 1

        if sum(
            1 for count in states.values()
            if count >= 2
        ) >= 2:

            informative += 1


# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

AlignIO.write(
    trimmed,
    OUTPUT,
    "fasta"
)

report = []

report.append(
    "DODDAGGUBBI - ALIGNMENT TRIMMING REPORT"
)

report.append("=" * 90)
report.append("")

report.append(
    f"Original sequences: {n}"
)

report.append(
    f"Original alignment length: {L}"
)

report.append(
    "Trimming rule: remove columns with >=50% gaps"
)

report.append(
    f"Columns removed: {len(removed_columns)}"
)

report.append(
    f"Columns retained: {len(keep_columns)}"
)

report.append("")

report.append(
    f"Trimmed constant sites: {constant}"
)

report.append(
    f"Trimmed variable sites: {variable}"
)

report.append(
    f"Trimmed parsimony-informative sites: "
    f"{informative}"
)

report.append("")

report.append(
    f"Removed fraction: "
    f"{len(removed_columns)/L*100:.2f}%"
)

report.append("")

report.append(
    "Original MAFFT alignment was NOT modified."
)

REPORT.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8"
)

print("=" * 90)
print(" ALIGNMENT TRIMMING COMPLETE")
print("=" * 90)

print("Sequences:", len(trimmed))
print("Original length:", L)
print("Trimmed length:", trimmed.get_alignment_length())
print("Columns removed:", len(removed_columns))
print("Columns retained:", len(keep_columns))
print("Variable sites:", variable)
print("Parsimony-informative sites:", informative)

print("")
print("Output:", OUTPUT)
print("Report:", REPORT)

print("")
print("ORIGINAL ALIGNMENT WAS NOT MODIFIED.")
