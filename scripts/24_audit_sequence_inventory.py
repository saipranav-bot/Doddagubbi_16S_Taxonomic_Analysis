from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

CANONICAL = META / "CANONICAL_KJU_METADATA_V1.csv"
FASTA_DIR = MASTER / "01_NCBI" / "fasta"
GB_DIR = MASTER / "01_NCBI" / "genbank"

OUT = META / "SEQUENCE_INVENTORY_AUDIT_V1.csv"
SUMMARY = META / "SEQUENCE_INVENTORY_AUDIT_SUMMARY_V1.txt"

print("=" * 90)
print(" DODDAGGUBBI - SEQUENCE INVENTORY AUDIT")
print("=" * 90)

canonical = pd.read_csv(CANONICAL, dtype=str).fillna("")

# ------------------------------------------------------------
# FASTA inventory
# ------------------------------------------------------------

fasta_records = {}

for f in sorted(FASTA_DIR.glob("*")):

    if not f.is_file():
        continue

    accession = f.name

    # Read FASTA headers
    headers = []

    try:
        with open(f, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                line = line.strip()

                if line.startswith(">"):
                    headers.append(line[1:])

    except Exception as e:
        headers.append(f"ERROR: {e}")

    fasta_records[accession] = {
        "file": f.name,
        "headers": " || ".join(headers),
        "header_count": len(headers),
    }

# ------------------------------------------------------------
# GenBank inventory
# ------------------------------------------------------------

gb_records = {}

for f in sorted(GB_DIR.glob("*.gb")):

    accession = f.name

    text = f.read_text(
        encoding="utf-8",
        errors="ignore"
    )

    version_match = re.search(
        r"VERSION\s+([A-Z0-9_.]+)",
        text
    )

    accession_match = re.search(
        r"ACCESSION\s+([A-Z0-9]+)",
        text
    )

    organism_match = re.search(
        r"^\s+ORGANISM\s+(.+)$",
        text,
        re.MULTILINE
    )

    kju_matches = sorted(
        set(
            re.findall(
                r"KJU\d{4}",
                text,
                flags=re.IGNORECASE
            )
        )
    )

    gb_records[accession] = {
        "file": f.name,
        "genbank_accession": (
            version_match.group(1)
            if version_match
            else (
                accession_match.group(1)
                if accession_match
                else ""
            )
        ),
        "organism": (
            organism_match.group(1).strip()
            if organism_match
            else ""
        ),
        "kju_ids_in_record": ";".join(kju_matches),
    }

# ------------------------------------------------------------
# Build accession-level inventory
# ------------------------------------------------------------

rows = []

all_accessions = sorted(
    set(
        canonical["canonical_accession"]
    )
    |
    set(fasta_records.keys())
    |
    set(gb_records.keys())
)

for accession in all_accessions:

    c = canonical[
        canonical["canonical_accession"] == accession
    ]

    if len(c) > 0:

        kju_id = ";".join(
            sorted(c["kju_id"].tolist())
        )

        canonical_organism = ";".join(
            sorted(
                set(
                    x for x in c["canonical_organism"]
                    if x
                )
            )
        )

        sequence_type = ";".join(
            sorted(
                set(
                    x for x in c["canonical_sequence_type"]
                    if x
                )
            )
        )

    else:
        kju_id = ""
        canonical_organism = ""
        sequence_type = ""

    fasta = fasta_records.get(
        accession,
        {}
    )

    gb = gb_records.get(
        accession,
        {}
    )

    rows.append({

        "accession": accession,

        "canonical_kju_id": kju_id,

        "canonical_organism": canonical_organism,

        "canonical_sequence_type": sequence_type,

        "fasta_file": fasta.get(
            "file",
            ""
        ),

        "fasta_header_count": fasta.get(
            "header_count",
            ""
        ),

        "fasta_headers": fasta.get(
            "headers",
            ""
        ),

        "genbank_file": gb.get(
            "file",
            ""
        ),

        "genbank_accession": gb.get(
            "genbank_accession",
            ""
        ),

        "genbank_organism": gb.get(
            "organism",
            ""
        ),

        "kju_ids_in_genbank": gb.get(
            "kju_ids_in_record",
            ""
        ),

    })

inventory = pd.DataFrame(rows)

# ------------------------------------------------------------
# Classification
# ------------------------------------------------------------

def classify(row):

    has_canonical = bool(
        row["canonical_kju_id"]
    )

    has_fasta = bool(
        row["fasta_file"]
    )

    has_gb = bool(
        row["genbank_file"]
    )

    if has_canonical and has_fasta and has_gb:
        return "CANONICAL_WITH_FASTA_AND_GENBANK"

    if has_canonical and has_fasta:
        return "CANONICAL_WITH_FASTA_ONLY"

    if has_canonical and has_gb:
        return "CANONICAL_WITH_GENBANK_ONLY"

    if has_canonical:
        return "CANONICAL_NO_SEQUENCE_FILE"

    if has_fasta and has_gb:
        return "ORPHAN_FASTA_AND_GENBANK"

    if has_fasta:
        return "ORPHAN_FASTA_ONLY"

    if has_gb:
        return "ORPHAN_GENBANK_ONLY"

    return "UNCLASSIFIED"


inventory["classification"] = inventory.apply(
    classify,
    axis=1
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

inventory.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = []

summary.append(
    "DODDAGGUBBI - SEQUENCE INVENTORY AUDIT"
)

summary.append(
    "=" * 80
)

summary.append(
    f"Canonical KJU records: {len(canonical)}"
)

summary.append(
    f"Canonical records with accession: "
    f"{(canonical['canonical_accession'] != '').sum()}"
)

summary.append(
    f"FASTA files: {len(fasta_records)}"
)

summary.append(
    f"GenBank files: {len(gb_records)}"
)

summary.append("")

summary.append(
    "CLASSIFICATION"
)

summary.append(
    "-" * 80
)

summary.append(
    inventory["classification"]
    .value_counts()
    .to_string()
)

summary.append("")

summary.append(
    "CANONICAL ACCESSIONS WITH NO FASTA"
)

summary.append(
    "-" * 80
)

missing_fasta = inventory[
    (inventory["canonical_kju_id"] != "")
    &
    (inventory["fasta_file"] == "")
]

if len(missing_fasta) > 0:
    summary.extend(
        missing_fasta[
            [
                "canonical_kju_id",
                "accession",
                "canonical_organism",
            ]
        ]
        .astype(str)
        .apply(
            lambda x:
            " | ".join(x),
            axis=1
        )
        .tolist()
    )
else:
    summary.append("NONE")

summary.append("")

summary.append(
    "CANONICAL ACCESSIONS WITH NO GENBANK"
)

summary.append(
    "-" * 80
)

missing_gb = inventory[
    (inventory["canonical_kju_id"] != "")
    &
    (inventory["genbank_file"] == "")
]

if len(missing_gb) > 0:
    summary.extend(
        missing_gb[
            [
                "canonical_kju_id",
                "accession",
                "canonical_organism",
            ]
        ]
        .astype(str)
        .apply(
            lambda x:
            " | ".join(x),
            axis=1
        )
        .tolist()
    )
else:
    summary.append("NONE")

summary.append("")

summary.append(
    "ORPHAN FILES"
)

summary.append(
    "-" * 80
)

orphans = inventory[
    inventory["canonical_kju_id"] == ""
]

if len(orphans) > 0:
    summary.extend(
        orphans[
            [
                "accession",
                "fasta_file",
                "genbank_file",
            ]
        ]
        .astype(str)
        .apply(
            lambda x:
            " | ".join(x),
            axis=1
        )
        .tolist()
    )
else:
    summary.append("NONE")

summary.append("")

summary.append(
    "NO FILES WERE MODIFIED."
)

SUMMARY.write_text(
    "\n".join(summary),
    encoding="utf-8"
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

print(
    "\nCanonical records:",
    len(canonical)
)

print(
    "Canonical accessions:",
    (canonical["canonical_accession"] != "").sum()
)

print(
    "FASTA files:",
    len(fasta_records)
)

print(
    "GenBank files:",
    len(gb_records)
)

print("\nClassification:")
print(
    inventory["classification"]
    .value_counts()
)

print("\nOutput:")
print(OUT)

print("\nSummary:")
print(SUMMARY)

print("\nMASTER DATA WAS NOT MODIFIED.")
