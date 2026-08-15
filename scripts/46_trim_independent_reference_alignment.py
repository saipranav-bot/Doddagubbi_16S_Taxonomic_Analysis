from pathlib import Path
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

MASTER = Path.home() / "Doddagubbi_Master"

INPUT = (
    MASTER
    / "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
      "REFERENCE_PHYLOGENY_V2/"
      "KJU109_NCBI65_INDEPENDENT_MAFFT.fasta"
)

OUTPUT = (
    MASTER
    / "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
      "REFERENCE_PHYLOGENY_V2/"
      "KJU109_NCBI65_INDEPENDENT_MAFFT_TRIMMED_50GAP.fasta"
)

REPORT = (
    MASTER
    / "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
      "REFERENCE_PHYLOGENY_V2/"
      "INDEPENDENT_REFERENCE_ALIGNMENT_TRIMMING_REPORT_V1.txt"
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

    gaps = sum(
        1 for x in column
        if x in gap_chars
    )

    gap_percent = (gaps / n) * 100

    if gap_percent >= 50:
        removed.append(i)
    else:
        keep.append(i)

trimmed_records = []

for record in aln:

    seq = "".join(
        str(record.seq[i])
        for i in keep
    )

    new_record = record[:]
    new_record.seq = record.seq.__class__(seq)

    trimmed_records.append(new_record)

trimmed = MultipleSeqAlignment(trimmed_records)

AlignIO.write(
    trimmed,
    OUTPUT,
    "fasta"
)

variable = 0
informative = 0
constant = 0

for i in range(trimmed.get_alignment_length()):

    states = []

    for record in trimmed:

        base = str(
            record.seq[i]
        ).upper()

        if base not in gap_chars:
            states.append(base)

    unique = set(states)

    if len(unique) <= 1:
        constant += 1
    else:
        variable += 1

    counts = {}

    for base in states:
        counts[base] = counts.get(base, 0) + 1

    repeated_states = sum(
        1 for count in counts.values()
        if count >= 2
    )

    if len(unique) >= 2 and repeated_states >= 2:
        informative += 1

report = [
    "DODDAGGUBBI - INDEPENDENT REFERENCE ALIGNMENT TRIMMING",
    "=" * 110,
    "",
    f"Sequences: {n}",
    f"Original alignment length: {L}",
    f"Trimmed alignment length: {trimmed.get_alignment_length()}",
    f"Columns removed: {len(removed)}",
    f"Columns retained: {len(keep)}",
    "",
    f"Constant sites: {constant}",
    f"Variable sites: {variable}",
    f"Parsimony-informative sites: {informative}",
    "",
    "Trimming rule:",
    "Columns containing >=50% gaps were removed.",
    "",
    f"Output: {OUTPUT}",
]

REPORT.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8"
)

print("=" * 110)
print(" INDEPENDENT REFERENCE ALIGNMENT TRIMMING COMPLETE")
print("=" * 110)

print("Sequences:", n)
print("Original length:", L)
print(
    "Trimmed length:",
    trimmed.get_alignment_length()
)
print("Columns removed:", len(removed))
print("Columns retained:", len(keep))
print("Constant sites:", constant)
print("Variable sites:", variable)
print(
    "Parsimony-informative sites:",
    informative
)

print("")
print("Output:", OUTPUT)
print("Report:", REPORT)

print("")
print("SOURCE DATA WAS NOT MODIFIED.")
