from pathlib import Path
from Bio import SeqIO

MASTER = Path.home() / "Doddagubbi_Master"

FASTA = (
    MASTER /
    "07_Master_Metadata/FINAL_109_SEQUENCE_QC_V2/"
    "FINAL_109_KJU_16S_QC_V2.fasta"
)

OUT = (
    MASTER /
    "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
    "REFERENCE_PHYLOGENY_V2/BLAST_CONFIRMATION_V1/"
    "DISCORDANT_9KJU_QUERIES.fasta"
)

TARGETS = {
    "KJU0007",
    "KJU0017",
    "KJU0088",
    "KJU0099",
    "KJU0105",
    "KJU0106",
    "KJU0208",
    "KJU0217",
    "KJU0221",
}

records = list(SeqIO.parse(FASTA, "fasta"))

selected = []

for r in records:
    kju = r.id.split("|")[0]
    if kju in TARGETS:
        selected.append(r)

if len(selected) != 9:
    raise SystemExit(
        f"ERROR: Expected 9 KJUs, found {len(selected)}"
    )

OUT.parent.mkdir(parents=True, exist_ok=True)

SeqIO.write(selected, OUT, "fasta")

print("=" * 90)
print("9-KJU BLAST QUERY PREPARATION")
print("=" * 90)

print("Queries:", len(selected))

for r in selected:
    print(r.id, len(r.seq))

print("")
print("Output:", OUT)
print("QUERY PREPARATION: PASS")
