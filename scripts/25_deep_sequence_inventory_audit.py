from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

CANONICAL_FILE = META / "CANONICAL_KJU_METADATA_V1.csv"

FASTA_DIR = MASTER / "01_NCBI" / "fasta"
GENBANK_DIR = MASTER / "01_NCBI" / "genbank"

OUT = META / "DEEP_SEQUENCE_INVENTORY_AUDIT_V1.csv"
SUMMARY = META / "DEEP_SEQUENCE_INVENTORY_AUDIT_SUMMARY_V1.txt"

print("=" * 90)
print(" DODDAGGUBBI - DEEP SEQUENCE INVENTORY AUDIT")
print("=" * 90)

canonical = pd.read_csv(CANONICAL_FILE, dtype=str).fillna("")

# ------------------------------------------------------------
# Canonical lookup
# ------------------------------------------------------------

canonical_by_acc = {}
canonical_by_kju = {}

for _, row in canonical.iterrows():

    kju = row["kju_id"].strip()
    acc = row.get("canonical_accession", "").strip()

    if kju:
        canonical_by_kju[kju] = row.to_dict()

    if acc:
        canonical_by_acc[acc] = row.to_dict()

# ------------------------------------------------------------
# Accession normalization
# ------------------------------------------------------------

def normalize_accession(x):

    x = x.strip()

    if not x:
        return ""

    # Remove common sequence-file extensions
    x = re.sub(r"\.(fasta|fa|fna|gb|gbk|genbank)$", "", x,
               flags=re.IGNORECASE)

    return x


def extract_accessions(text):

    """
    Extract accession-like identifiers from headers.
    Handles examples such as:

    PZ666011
    PZ666011.1
    """

    pattern = r"\b[A-Z]{1,5}\d{5,}(?:\.\d+)?\b"

    return sorted(set(re.findall(pattern, text)))


# ------------------------------------------------------------
# FASTA audit
# ------------------------------------------------------------

records = []

if FASTA_DIR.exists():

    for f in sorted(FASTA_DIR.iterdir()):

        if not f.is_file():
            continue

        text = f.read_text(
            encoding="utf-8",
            errors="ignore"
        )

        headers = [
            line.strip()[1:]
            for line in text.splitlines()
            if line.startswith(">")
        ]

        header_text = " | ".join(headers)

        header_accessions = extract_accessions(header_text)

        filename_acc = normalize_accession(f.name)

        matched_acc = ""
        matched_kju = ""
        match_type = ""

        # filename accession
        for acc in [filename_acc] + header_accessions:

            acc_base = normalize_accession(acc)

            for canonical_acc, row in canonical_by_acc.items():

                canonical_base = normalize_accession(canonical_acc)

                if (
                    acc_base == canonical_base
                    or acc_base.split(".")[0] == canonical_base.split(".")[0]
                ):

                    matched_acc = canonical_acc
                    matched_kju = row["kju_id"]

                    if acc == filename_acc:
                        match_type = "FILENAME_ACCESSION_MATCH"
                    else:
                        match_type = "HEADER_ACCESSION_MATCH"

                    break

            if matched_kju:
                break

        if matched_kju:
            classification = "MATCHED_CANONICAL"
        elif header_accessions:
            classification = "ACCESSION_NOT_IN_CANONICAL"
        else:
            classification = "NO_ACCESSION_FOUND"

        records.append({
            "file_type": "FASTA",
            "file_name": f.name,
            "file_path": str(f),
            "filename_accession": filename_acc,
            "header_accessions": ";".join(header_accessions),
            "header_text": header_text,
            "matched_canonical_accession": matched_acc,
            "matched_kju_id": matched_kju,
            "match_type": match_type,
            "classification": classification,
        })

# ------------------------------------------------------------
# GenBank audit
# ------------------------------------------------------------

if GENBANK_DIR.exists():

    for f in sorted(GENBANK_DIR.iterdir()):

        if not f.is_file():
            continue

        text = f.read_text(
            encoding="utf-8",
            errors="ignore"
        )

        # GenBank ACCESSION
        accession_matches = re.findall(
            r"^ACCESSION\s+(\S+)",
            text,
            flags=re.MULTILINE
        )

        # GenBank VERSION
        version_matches = re.findall(
            r"^VERSION\s+(\S+)",
            text,
            flags=re.MULTILINE
        )

        # Definition
        definition_matches = re.findall(
            r"^DEFINITION\s+(.+)",
            text,
            flags=re.MULTILINE
        )

        filename_acc = normalize_accession(f.name)

        accessions = []

        for x in accession_matches + version_matches:
            if x not in accessions:
                accessions.append(x)

        matched_acc = ""
        matched_kju = ""
        match_type = ""

        candidates = [filename_acc] + accessions

        for acc in candidates:

            acc_base = normalize_accession(acc)

            for canonical_acc, row in canonical_by_acc.items():

                canonical_base = normalize_accession(canonical_acc)

                if (
                    acc_base == canonical_base
                    or acc_base.split(".")[0] == canonical_base.split(".")[0]
                ):

                    matched_acc = canonical_acc
                    matched_kju = row["kju_id"]

                    if acc == filename_acc:
                        match_type = "FILENAME_ACCESSION_MATCH"
                    else:
                        match_type = "GENBANK_ACCESSION_MATCH"

                    break

            if matched_kju:
                break

        definition = (
            definition_matches[0].strip()
            if definition_matches
            else ""
        )

        if matched_kju:
            classification = "MATCHED_CANONICAL"
        elif accessions:
            classification = "ACCESSION_NOT_IN_CANONICAL"
        else:
            classification = "NO_ACCESSION_FOUND"

        records.append({
            "file_type": "GENBANK",
            "file_name": f.name,
            "file_path": str(f),
            "filename_accession": filename_acc,
            "header_accessions": ";".join(accessions),
            "header_text": definition,
            "matched_canonical_accession": matched_acc,
            "matched_kju_id": matched_kju,
            "match_type": match_type,
            "classification": classification,
        })

# ------------------------------------------------------------
# Save audit
# ------------------------------------------------------------

audit = pd.DataFrame(records)

audit.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = []

summary.append("DODDAGGUBBI - DEEP SEQUENCE INVENTORY AUDIT")
summary.append("=" * 90)
summary.append("")

summary.append(f"Canonical records: {len(canonical)}")
summary.append(
    f"Canonical records with accession: "
    f"{(canonical['canonical_accession'] != '').sum()}"
)
summary.append(
    f"Canonical records without accession: "
    f"{(canonical['canonical_accession'] == '').sum()}"
)

summary.append("")
summary.append("SEQUENCE FILES")
summary.append("-" * 90)

summary.append(
    f"FASTA files: "
    f"{(audit['file_type'] == 'FASTA').sum()}"
)

summary.append(
    f"GenBank files: "
    f"{(audit['file_type'] == 'GENBANK').sum()}"
)

summary.append("")
summary.append("CLASSIFICATION")
summary.append("-" * 90)

if len(audit):
    summary.append(
        audit["classification"]
        .value_counts()
        .to_string()
    )

summary.append("")
summary.append("MATCH TYPE")
summary.append("-" * 90)

if len(audit):
    summary.append(
        audit["match_type"]
        .replace("", "NO_MATCH")
        .value_counts()
        .to_string()
    )

summary.append("")
summary.append("IMPORTANT")
summary.append("-" * 90)
summary.append(
    "NO FILES WERE RENAMED."
)
summary.append(
    "NO SEQUENCES WERE MODIFIED."
)
summary.append(
    "NO MASTER/CANONICAL METADATA WAS MODIFIED."
)
summary.append(
    "This audit must be reviewed before any rename operation."
)

SUMMARY.write_text(
    "\n".join(summary),
    encoding="utf-8"
)

print("")
print("=" * 90)
print(" DEEP AUDIT COMPLETE")
print("=" * 90)

print(f"Audit   : {OUT}")
print(f"Summary : {SUMMARY}")

print("")
print("CLASSIFICATION:")
print(
    audit["classification"]
    .value_counts()
    .to_string()
)

print("")
print("NO FILES WERE MODIFIED.")
