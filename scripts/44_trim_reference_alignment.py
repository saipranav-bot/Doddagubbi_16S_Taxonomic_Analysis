from pathlib import Path
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

MASTER = Path.home() / "Doddagubbi_Master"

INPUT = (
    MASTER
    / "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
      "REFERENCE_PHYLOGENY_V1/KJU109_NCBI65_MAFFT.fasta"
)

OUTPUT = (
    MASTER
    / "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
      "REFERENCE_PHYLOGENY_V1/"
      "KJU109_NCBI65_MAFFT_TRIMMED_50GAP.fasta"
)

REPORT = (
    MASTER
    / "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
      "REFERENCE_PHYLOGENY_V1/"
      "REFERENCE_ALIGNMENT_TRIMMING_REPORT_V1.txt"
)

aln = AlignIO.read(INPUT, "fasta")

n = len(aln)
L = aln.get_alignment_length()

gap_chars = set("-.")

keep = []
removed = []

for i in range(L):

    column = [
        str(record.seq[i]).upper()
        for record in aln
    ]

    gap_count = sum(
        1 for x in column
        if x in gap_chars
    )

    gap_percent = (
        gap_count / n
    ) * 100

    if gap_percent >= 50:
        removed.append(i)
    else:
        keep.append(i)

trimmed_records = []

for record in aln:

    sequence = "".join(
        str(record.seq[i])
        for i in keep
    )

    record_copy = record[:]
    record_copy.seq = record.seq.__class__(sequence)

    trimmed_records.append(record_copy)

trimmed = MultipleSeqAlignment(
    trimmed_records
)

AlignIO.write(
    trimmed,
    OUTPUT,
    "fasta"
)

variable = 0
informative = 0
constant = 0

for i in range(
    trimmed.get_alignment_length()
):

    states = set(
        str(record.seq[i]).upper()
        for record in trimmed
        if str(record.seq[i]).upper()
        not in gap_chars
    )

    if len(states) <= 1:
        constant += 1
    else:
        variable += 1

    counts = {}

    for record in trimmed:

        base = str(
            record.seq[i]
        ).upper()

        if base in gap_chars:
            continue

        counts[base] = (
            counts.get(base, 0) + 1
        )

    informative_states = sum(
        1 for count in counts.values()
        if count >= 2
    )

    if (
        len(counts) >= 2
        and
        informative_states >= 2
    ):
        informative += 1

report = [
    "DODDAGGUBBI - COMBINED REFERENCE ALIGNMENT TRIMMING",
    "=" * 100,
    "",
    f"Sequences: {n}",
    f"Original alignment length: {L}",
    f"Trimmed alignment length: "
    f"{trimmed.get_alignment_length()}",
    f"Columns removed: {len(removed)}",
    f"Columns retained: {len(keep)}",
    "",
    f"Constant sites: {constant}",
    f"Variable sites: {variable}",
    f"Parsimony-informative sites: {informative}",
    "",
    "Rule:",
    "Columns with >=50% gaps were removed.",
    "",
    f"Output: {OUTPUT}",
]

REPORT.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8"
)

print("=" * 100)
print("COMBINED ALIGNMENT TRIMMING COMPLETE")
print("=" * 100)
print("Sequences:", n)
print("Original length:", L)
print(
    "Trimmed length:",
    trimmed.get_alignment_length()
)
print("Columns removed:", len(removed))
print("Columns retained:", len(keep))
print("Variable sites:", variable)
print(
    "Parsimony-informative sites:",
    informative
)
print("Output:", OUTPUT)
print("Report:", REPORT)
