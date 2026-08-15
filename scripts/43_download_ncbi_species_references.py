from pathlib import Path
import requests
import time
import re
import pandas as pd
from Bio import Entrez, SeqIO

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

OUTDIR = META / "FINAL_109_PHYLOGENY_V1" / "NCBI_SPECIES_REFERENCES_V1"
OUTDIR.mkdir(parents=True, exist_ok=True)

FASTA = OUTDIR / "NCBI_REFERENCE_16S.fasta"
TABLE = OUTDIR / "NCBI_REFERENCE_METADATA.csv"
SUMMARY = OUTDIR / "NCBI_REFERENCE_SUMMARY.txt"

# ------------------------------------------------------------
# NCBI identification
# ------------------------------------------------------------

Entrez.email = "pranav@example.com"

# ------------------------------------------------------------
# Target species
# ------------------------------------------------------------

TARGET_SPECIES = {
    "Aeromonas caviae",
    "Aeromonas hydrophila",
    "Aeromonas veronii",
    "Aeromonas salmonicida",

    "Bacillus cereus",
    "Bacillus wiedmannii",
    "Bacillus subtilis",
    "Bacillus velezensis",

    "Citrobacter freundii",
    "Citrobacter braakii",
    "Citrobacter koseri",

    "Shewanella putrefaciens",
    "Shewanella xiamenensis",
}

# Maximum references per species
MAX_PER_SPECIES = 5

print("=" * 110)
print(" DODDAGGUBBI - NCBI 16S REFERENCE SEQUENCE COLLECTION")
print("=" * 110)

all_records = []
seen_accessions = set()

# ------------------------------------------------------------
# Search NCBI
# ------------------------------------------------------------

for species in sorted(TARGET_SPECIES):

    query = (
        f'"{species}"[Organism] '
        f'AND ("16S ribosomal RNA"[Title] '
        f'OR "16S rRNA"[All Fields]) '
        f'AND 1000:1700[SLEN]'
    )

    print(f"\nSEARCHING: {species}")
    print(f"Query: {query}")

    try:
        handle = Entrez.esearch(
            db="nuccore",
            term=query,
            retmode="xml",
            retmax=30
        )

        result = Entrez.read(handle)
        handle.close()

        ids = result.get("IdList", [])

    except Exception as e:
        print("SEARCH ERROR:", e)
        continue

    print("Hits:", len(ids))

    if not ids:
        continue

    # --------------------------------------------------------
    # Fetch candidate records
    # --------------------------------------------------------

    try:
        handle = Entrez.efetch(
            db="nuccore",
            id=",".join(ids),
            rettype="gb",
            retmode="text"
        )

        records = list(
            SeqIO.parse(handle, "genbank")
        )

        handle.close()

    except Exception as e:
        print("FETCH ERROR:", e)
        continue

    selected = 0

    for record in records:

        accession = str(
            record.id
        ).strip()

        if accession in seen_accessions:
            continue

        sequence = str(
            record.seq
        ).upper()

        # ----------------------------------------------------
        # Basic sequence checks
        # ----------------------------------------------------

        if len(sequence) < 1000:
            continue

        if len(sequence) > 1700:
            continue

        if not re.fullmatch(
            r"[ACGTRYSWKMBDHVN]+",
            sequence
        ):
            continue

        # ----------------------------------------------------
        # Extract organism
        # ----------------------------------------------------

        organism = ""

        for feature in record.features:

            if feature.type != "source":
                continue

            qualifiers = feature.qualifiers

            if "organism" in qualifiers:
                organism = (
                    qualifiers["organism"][0]
                    .strip()
                )

            break

        if not organism:
            continue

        # Must match target species exactly
        if organism.lower() != species.lower():
            continue

        # ----------------------------------------------------
        # Extract 16S annotation evidence
        # ----------------------------------------------------

        description = str(
            record.description
        )

        text = (
            description
            + " "
            + str(record.annotations)
        ).lower()

        if (
            "16s" not in text
            and
            "ribosomal rna" not in text
        ):
            continue

        # ----------------------------------------------------
        # Avoid duplicated sequences from same species
        # ----------------------------------------------------

        seq_hash = __import__("hashlib").sha256(
            sequence.encode()
        ).hexdigest()

        # Don't take multiple identical reference sequences
        if any(
            x["sequence_sha256"] == seq_hash
            for x in all_records
        ):
            continue

        all_records.append({
            "accession": accession,
            "species": species,
            "organism": organism,
            "length": len(sequence),
            "description": description,
            "sequence_sha256": seq_hash,
            "sequence": sequence
        })

        seen_accessions.add(accession)

        selected += 1

        if selected >= MAX_PER_SPECIES:
            break

    print(
        "Selected:",
        selected
    )

    time.sleep(0.34)

# ------------------------------------------------------------
# Save FASTA
# ------------------------------------------------------------

with open(FASTA, "w") as fh:

    for r in all_records:

        safe_species = (
            r["species"]
            .replace(" ", "_")
        )

        fh.write(
            f">{r['accession']}|NCBI_REF|"
            f"{safe_species}\n"
        )

        seq = r["sequence"]

        for i in range(0, len(seq), 80):
            fh.write(
                seq[i:i+80] + "\n"
            )

# ------------------------------------------------------------
# Save metadata
# ------------------------------------------------------------

metadata = []

for r in all_records:

    metadata.append({
        "accession": r["accession"],
        "species": r["species"],
        "organism": r["organism"],
        "length": r["length"],
        "description": r["description"],
        "sequence_sha256":
            r["sequence_sha256"]
    })

pd.DataFrame(metadata).to_csv(
    TABLE,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = []

summary.append(
    "DODDAGGUBBI - NCBI REFERENCE 16S COLLECTION"
)

summary.append("=" * 110)
summary.append("")

summary.append(
    f"Target species: {len(TARGET_SPECIES)}"
)

summary.append(
    f"Reference sequences collected: {len(all_records)}"
)

summary.append(
    f"Unique species represented: "
    f"{len(set(r['species'] for r in all_records))}"
)

summary.append("")
summary.append("REFERENCES PER SPECIES")
summary.append("-" * 110)

counts = {}

for r in all_records:
    counts[r["species"]] = (
        counts.get(r["species"], 0) + 1
    )

for species in sorted(TARGET_SPECIES):

    summary.append(
        f"{species}: {counts.get(species, 0)}"
    )

summary.append("")
summary.append("QC")
summary.append("-" * 110)

summary.append(
    "All retained sequences were between 1000 and 1700 bp."
)

summary.append(
    "All retained sequences passed the basic IUPAC DNA check."
)

summary.append(
    "Only records annotated to the requested organism/species "
    "were retained."
)

summary.append("")
summary.append("IMPORTANT")
summary.append("-" * 110)

summary.append(
    "These are independent NCBI reference sequences."
)

summary.append(
    "They were not added to the 109-KJU final dataset."
)

summary.append(
    "The original project metadata and sequence files were not modified."
)

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("")
print("=" * 110)
print("NCBI REFERENCE COLLECTION COMPLETE")
print("=" * 110)

print("Reference sequences:", len(all_records))
print("FASTA:", FASTA)
print("Metadata:", TABLE)
print("Summary:", SUMMARY)

print("")
print("NO DODDAGGUBBI SOURCE DATA WAS MODIFIED.")
