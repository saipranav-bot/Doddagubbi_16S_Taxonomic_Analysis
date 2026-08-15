from pathlib import Path
import pandas as pd
import re
import time
import hashlib
from Bio import Entrez, SeqIO
from Bio.Seq import UndefinedSequenceError

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

KJU_FILE = META / "FINAL_KJU_16S_DATASET_VALIDATED_V1.csv"

OUTDIR = (
    META
    / "FINAL_109_PHYLOGENY_V1"
    / "INDEPENDENT_NCBI_REFERENCES_V2"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

FASTA = OUTDIR / "INDEPENDENT_NCBI_16S_REFERENCES.fasta"
TABLE = OUTDIR / "INDEPENDENT_NCBI_REFERENCE_METADATA.csv"
SUMMARY = OUTDIR / "INDEPENDENT_NCBI_REFERENCE_SUMMARY.txt"

Entrez.email = "pranav@example.com"

TARGET_SPECIES = {
    "Aeromonas caviae",
    "Aeromonas hydrophila",
    "Aeromonas salmonicida",
    "Aeromonas veronii",
    "Bacillus cereus",
    "Bacillus subtilis",
    "Bacillus velezensis",
    "Bacillus wiedmannii",
    "Citrobacter braakii",
    "Citrobacter freundii",
    "Citrobacter koseri",
    "Shewanella putrefaciens",
    "Shewanella xiamenensis",
}

MAX_PER_SPECIES = 5

print("=" * 110)
print(" DODDAGGUBBI - TRULY INDEPENDENT NCBI REFERENCE COLLECTION")
print("=" * 110)

# ------------------------------------------------------------
# Load KJU accessions
# ------------------------------------------------------------

kju = pd.read_csv(
    KJU_FILE,
    dtype=str
).fillna("")

kju_accessions = set(
    kju["accession"]
    .astype(str)
    .str.strip()
    .str.upper()
)

print(
    "KJU accessions protected:",
    len(kju_accessions)
)

records_out = []
seen_accessions = set()
seen_hashes = set()

skipped_undefined = []
skipped_kju = []
skipped_duplicate = []
skipped_invalid = []

# ------------------------------------------------------------
# Search NCBI
# ------------------------------------------------------------

for species in sorted(TARGET_SPECIES):

    print("")
    print("=" * 80)
    print("SEARCH:", species)
    print("=" * 80)

    query = (
        f'"{species}"[Organism] '
        f'AND ("16S ribosomal RNA"[Title] '
        f'OR "16S rRNA"[All Fields]) '
        f'AND 1000:1700[SLEN]'
    )

    try:
        h = Entrez.esearch(
            db="nuccore",
            term=query,
            retmode="xml",
            retmax=100
        )

        result = Entrez.read(h)
        h.close()

        ids = result.get("IdList", [])

    except Exception as e:
        print("SEARCH ERROR:", e)
        continue

    print("NCBI hits:", len(ids))

    selected = 0

    # --------------------------------------------------------
    # Fetch candidates in batches
    # --------------------------------------------------------

    for start in range(0, len(ids), 25):

        batch = ids[start:start + 25]

        try:

            h = Entrez.efetch(
                db="nuccore",
                id=",".join(batch),
                rettype="gb",
                retmode="text"
            )

            gb_records = list(
                SeqIO.parse(h, "genbank")
            )

            h.close()

        except Exception as e:

            print("FETCH ERROR:", e)
            continue

        for record in gb_records:

            accession = str(record.id).strip()
            accession_upper = accession.upper()

            # ------------------------------------------------
            # Exclude KJU-owned accessions
            # ------------------------------------------------

            if accession_upper in kju_accessions:

                skipped_kju.append(accession)
                continue

            # ------------------------------------------------
            # Exclude duplicate accession
            # ------------------------------------------------

            if accession_upper in seen_accessions:

                skipped_duplicate.append(accession)
                continue

            # ------------------------------------------------
            # Safely obtain sequence
            # ------------------------------------------------

            try:

                sequence = str(
                    record.seq
                ).upper()

            except (
                UndefinedSequenceError,
                ValueError,
                AttributeError
            ):

                skipped_undefined.append(accession)
                continue

            # ------------------------------------------------
            # Basic sequence validation
            # ------------------------------------------------

            if not (
                1000 <= len(sequence) <= 1700
            ):

                skipped_invalid.append(
                    f"{accession}:length={len(sequence)}"
                )
                continue

            if not re.fullmatch(
                r"[ACGTRYSWKMBDHVN]+",
                sequence
            ):

                skipped_invalid.append(
                    f"{accession}:invalid_DNA"
                )
                continue

            # ------------------------------------------------
            # Extract organism
            # ------------------------------------------------

            organism = ""

            for feature in record.features:

                if feature.type == "source":

                    if "organism" in feature.qualifiers:

                        organism = (
                            feature.qualifiers["organism"][0]
                            .strip()
                        )

                    break

            if not organism:
                skipped_invalid.append(
                    f"{accession}:missing_organism"
                )
                continue

            if organism.lower() != species.lower():
                skipped_invalid.append(
                    f"{accession}:organism_mismatch"
                )
                continue

            # ------------------------------------------------
            # Confirm 16S annotation
            # ------------------------------------------------

            description = str(
                record.description
            )

            feature_text = " ".join(
                str(feature.qualifiers)
                for feature in record.features
            )

            combined_text = (
                description
                + " "
                + feature_text
            ).lower()

            if (
                "16s" not in combined_text
                and
                "ribosomal rna" not in combined_text
            ):

                skipped_invalid.append(
                    f"{accession}:not_16S"
                )
                continue

            # ------------------------------------------------
            # Sequence duplicate check
            # ------------------------------------------------

            sha = hashlib.sha256(
                sequence.encode()
            ).hexdigest()

            if sha in seen_hashes:

                skipped_duplicate.append(
                    f"{accession}:same_sequence"
                )
                continue

            # ------------------------------------------------
            # Store
            # ------------------------------------------------

            records_out.append({

                "accession": accession,
                "species": species,
                "organism": organism,
                "length": len(sequence),
                "description": description,
                "sequence_sha256": sha,
                "sequence": sequence,

            })

            seen_accessions.add(
                accession_upper
            )

            seen_hashes.add(sha)

            selected += 1

            if selected >= MAX_PER_SPECIES:
                break

        if selected >= MAX_PER_SPECIES:
            break

        time.sleep(0.34)

    print(
        "Independent references selected:",
        selected
    )

# ------------------------------------------------------------
# Save FASTA
# ------------------------------------------------------------

with open(FASTA, "w") as fh:

    for r in records_out:

        species_safe = (
            r["species"]
            .replace(" ", "_")
        )

        fh.write(
            f">{r['accession']}|NCBI_REF|"
            f"{species_safe}\n"
        )

        seq = r["sequence"]

        for i in range(0, len(seq), 80):

            fh.write(
                seq[i:i + 80] + "\n"
            )

# ------------------------------------------------------------
# Save metadata
# ------------------------------------------------------------

metadata = []

for r in records_out:

    metadata.append({
        "accession": r["accession"],
        "species": r["species"],
        "organism": r["organism"],
        "length": r["length"],
        "description": r["description"],
        "sequence_sha256": r["sequence_sha256"],
    })

pd.DataFrame(metadata).to_csv(
    TABLE,
    index=False
)

# ------------------------------------------------------------
# Counts
# ------------------------------------------------------------

if metadata:

    metadata_df = pd.DataFrame(metadata)

    counts = (
        metadata_df
        .groupby("species")
        .size()
    )

else:

    counts = pd.Series(
        dtype=int
    )

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = [

    "DODDAGGUBBI - TRULY INDEPENDENT NCBI REFERENCES",

    "=" * 110,
    "",

    f"KJU accessions protected: "
    f"{len(kju_accessions)}",

    f"Independent reference sequences: "
    f"{len(records_out)}",

    f"Species represented: "
    f"{len(counts)}",

    "",
    "REFERENCES PER SPECIES",
    "-" * 110,
]

for species in sorted(TARGET_SPECIES):

    summary.append(
        f"{species}: "
        f"{int(counts.get(species, 0))}"
    )

summary += [

    "",
    "EXCLUSIONS",
    "-" * 110,

    f"Excluded KJU-owned accessions: "
    f"{len(skipped_kju)}",

    f"Skipped undefined-sequence records: "
    f"{len(skipped_undefined)}",

    f"Skipped duplicates: "
    f"{len(skipped_duplicate)}",

    f"Skipped invalid/non-matching records: "
    f"{len(skipped_invalid)}",

    "",
    "QC",
    "-" * 110,

    "Length range: 1000-1700 bp",
    "IUPAC DNA alphabet required",
    "Exact organism/species match required",
    "16S/ribosomal RNA annotation required",
    "Any accession already present in the 109-KJU "
    "dataset was excluded.",

    "",
    "IMPORTANT",
    "-" * 110,

    "These references are intended to be independent "
    "of the Doddagubbi KJU sequence dataset.",

    "The 109-KJU dataset was not modified.",

]

if skipped_undefined:

    summary += [
        "",
        "UNDEFINED-SEQUENCE ACCESSIONS",
        "-" * 110,
    ]

    summary.extend(
        sorted(set(skipped_undefined))
    )

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("")
print("=" * 110)
print("INDEPENDENT REFERENCE COLLECTION COMPLETE")
print("=" * 110)

print(
    "Independent references:",
    len(records_out)
)

print(
    "Species represented:",
    len(counts)
)

print("")
print("FASTA:", FASTA)
print("Metadata:", TABLE)
print("Summary:", SUMMARY)

print("")
print("NO KJU SOURCE FILES WERE MODIFIED.")
