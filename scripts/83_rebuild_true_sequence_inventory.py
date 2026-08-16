#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
import csv


# ============================================================
# STEP 83
# REBUILD TRUE SEQUENCE INVENTORY
#
# PURPOSE
# ============================================================
#
# STEP 82 established that:
#
#     Sequence_file = FILE PATH
#
# and NOT the nucleotide sequence itself.
#
# Therefore Steps 78-81 could not correctly determine sequence
# identity from Sequence_file.
#
# This step:
#
#   1. Reads STEP77 raw sequence records
#   2. Opens the actual FASTA / FASTQ / GenBank files
#   3. Extracts REAL nucleotide sequences
#   4. Calculates a new SHA256 from the nucleotide sequence
#   5. Calculates sequence length
#   6. Extracts sequence accession
#   7. Maps sequence records back to KJU IDs
#   8. Detects exact duplicate sequences
#   9. Detects multiple sequences per KJU
#  10. Detects same accession / different sequence
#  11. Detects different accession / identical sequence
#  12. Produces a clean true sequence inventory
#
# IMPORTANT
# ============================================================
#
# NEVER modifies:
#
#     FINAL_KJU_MASTER_FREEZE_V1
#
# The frozen V1 master remains untouched.
#
# This step only creates new QC files.
#
# ============================================================


from pathlib import Path
import pandas as pd
import hashlib
import re
import gzip


# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "final_integrated_master"
)

STEP77 = ROOT / "STEP77_SEQUENCE_QC"

RAW_INPUT = STEP77 / "STEP77_RAW_SEQUENCE_RECORDS.csv"

OUT = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"
OUT.mkdir(parents=True, exist_ok=True)

MASTER = (
    ROOT
    / "FINAL_KJU_MASTER_FREEZE_V1"
    / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
)


TRUE_INVENTORY = OUT / "STEP83_TRUE_SEQUENCE_INVENTORY.csv"
TRUE_RECORDS = OUT / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
DUPLICATES = OUT / "STEP83_EXACT_DUPLICATE_SEQUENCES.csv"
KJU_AUDIT = OUT / "STEP83_KJU_SEQUENCE_AUDIT.csv"
ACCESSION_AUDIT = OUT / "STEP83_ACCESSION_SEQUENCE_AUDIT.csv"
MISSING = OUT / "STEP83_MISSING_SEQUENCE_FILES.csv"
ERRORS = OUT / "STEP83_SEQUENCE_FILE_ERRORS.csv"
REPORT = OUT / "STEP83_TRUE_SEQUENCE_INVENTORY_REPORT.txt"


# ============================================================
# HELPERS
# ============================================================

def clean(value):
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    value = str(value).strip()

    if value.lower() in {
        "",
        "nan",
        "none",
        "null",
        "na",
        "n/a",
    }:
        return ""

    return value


def normalize_accession(value):
    value = clean(value)

    if not value:
        return ""

    value = value.strip()

    # Remove version suffix for normalized comparison
    # Example:
    #     PZ666011.1 -> PZ666011
    value = re.sub(r"\.\d+$", "", value)

    return value.upper()


def normalize_sequence(sequence):
    """
    Normalize a nucleotide sequence.

    Removes:
      - whitespace
      - FASTA formatting
      - digits

    Keeps:
      A C G T U N and IUPAC ambiguity symbols.
    """

    if sequence is None:
        return ""

    sequence = str(sequence).upper()

    sequence = re.sub(r"\s+", "", sequence)
    sequence = re.sub(r"[0-9]", "", sequence)

    # Remove accidental FASTA marker
    sequence = sequence.replace(">", "")

    # Keep valid nucleotide/IUPAC symbols only
    sequence = re.sub(r"[^ACGTUNRYSWKMBDHV]", "", sequence)

    return sequence


def sequence_sha256(sequence):
    return hashlib.sha256(
        sequence.encode("utf-8")
    ).hexdigest()


def open_text(path):
    """
    Open plain-text or gzipped sequence file.
    """

    path = Path(path)

    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")

    return open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
    )


def parse_fasta(path):
    """
    Simple FASTA parser.

    Returns:
        list of dictionaries:
            Record_ID
            Description
            Sequence
    """

    records = []

    current_id = None
    current_description = ""
    chunks = []

    with open_text(path) as handle:

        for line in handle:

            line = line.rstrip("\n\r")

            if not line:
                continue

            if line.startswith(">"):

                # Save previous record
                if current_id is not None:

                    sequence = normalize_sequence(
                        "".join(chunks)
                    )

                    records.append(
                        {
                            "Record_ID": current_id,
                            "Description": current_description,
                            "Sequence": sequence,
                        }
                    )

                header = line[1:].strip()

                parts = header.split(None, 1)

                current_id = parts[0]

                if len(parts) > 1:
                    current_description = parts[1]
                else:
                    current_description = ""

                chunks = []

            else:
                chunks.append(line.strip())

    # Save final record
    if current_id is not None:

        sequence = normalize_sequence(
            "".join(chunks)
        )

        records.append(
            {
                "Record_ID": current_id,
                "Description": current_description,
                "Sequence": sequence,
            }
        )

    return records


def parse_fastq(path):
    """
    Basic FASTQ parser.
    """

    records = []

    with open_text(path) as handle:

        while True:

            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence:
                break

            header = header.strip()
            sequence = normalize_sequence(sequence)

            if header.startswith("@"):
                header = header[1:]

            parts = header.split(None, 1)

            record_id = parts[0]

            description = (
                parts[1]
                if len(parts) > 1
                else ""
            )

            records.append(
                {
                    "Record_ID": record_id,
                    "Description": description,
                    "Sequence": sequence,
                }
            )

    return records


def parse_sequence_file(path):

    suffixes = [
        suffix.lower()
        for suffix in Path(path).suffixes
    ]

    name = Path(path).name.lower()

    if (
        ".fastq" in suffixes
        or ".fq" in suffixes
        or name.endswith(".fastq.gz")
        or name.endswith(".fq.gz")
    ):
        return parse_fastq(path)

    # Default to FASTA for .fa/.fasta/.fas/.fna
    return parse_fasta(path)


def extract_accession(record_id, description):

    text = (
        clean(record_id)
        + " "
        + clean(description)
    )

    # Common GenBank accession-like identifiers
    patterns = [

        # Explicit accession labels
        r"(?:accession|acc|gb)[=: ]+([A-Z]{1,4}\d{5,9}(?:\.\d+)?)",

        # PZ accessions
        r"\b(PZ\d{5,}(?:\.\d+)?)\b",

        # General GenBank-style accessions
        r"\b([A-Z]{1,4}\d{5,9}(?:\.\d+)?)\b",
    ]

    for pattern in patterns:

        match = re.search(
            pattern,
            text,
            flags=re.IGNORECASE,
        )

        if match:

            return (
                match.group(1).upper()
            )

    return ""


def extract_kju(text):

    text = clean(text)

    match = re.search(
        r"\bKJU\d{4}\b",
        text,
        flags=re.IGNORECASE,
    )

    if match:
        return match.group(0).upper()

    return ""


# ============================================================
# START
# ============================================================

print("=" * 110)
print("STEP 83: REBUILD TRUE SEQUENCE INVENTORY")
print("=" * 110)


# ============================================================
# LOAD RAW STEP 77
# ============================================================

if not RAW_INPUT.exists():
    raise FileNotFoundError(
        f"Missing STEP77 raw sequence file:\n{RAW_INPUT}"
    )


raw = pd.read_csv(
    RAW_INPUT,
    dtype=str,
    keep_default_na=False,
)

raw.columns = [
    str(c).strip()
    for c in raw.columns
]


# ============================================================
# LOAD MASTER
# ============================================================

if not MASTER.exists():
    raise FileNotFoundError(
        f"Missing frozen master:\n{MASTER}"
    )


master = pd.read_csv(
    MASTER,
    dtype=str,
    keep_default_na=False,
)

master.columns = [
    str(c).strip()
    for c in master.columns
]


# ============================================================
# DETECT COLUMNS
# ============================================================

def detect_column(df, candidates, required=True):

    lookup = {
        str(c).strip().lower(): c
        for c in df.columns
    }

    for candidate in candidates:

        key = candidate.lower()

        if key in lookup:
            return lookup[key]

    if required:
        raise ValueError(
            "Could not detect required column. "
            f"Tried: {candidates}\n"
            f"Available: {list(df.columns)}"
        )

    return None


raw_kju_col = detect_column(
    raw,
    [
        "KJU_ID",
        "kju_id",
        "KJU",
        "kju",
    ],
)

raw_file_col = detect_column(
    raw,
    [
        "Sequence_file",
        "sequence_file",
        "File",
        "file",
        "Sequence_path",
    ],
)

raw_record_col = detect_column(
    raw,
    [
        "Record_ID",
        "record_id",
        "Sequence_ID",
        "sequence_id",
    ],
    required=False,
)

raw_acc_col = detect_column(
    raw,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "Accession",
        "accession",
    ],
    required=False,
)


master_kju_col = detect_column(
    master,
    [
        "KJU_ID",
        "kju_id",
    ],
)


print()
print("INPUTS")
print("-" * 110)

print(f"Frozen master : {MASTER}")
print(f"Raw inventory : {RAW_INPUT}")

print()
print("INPUT SIZES")
print("-" * 110)

print(f"Master rows : {len(master)}")
print(f"Raw rows    : {len(raw)}")

print()
print("DETECTED COLUMNS")
print("-" * 110)

print(f"Raw KJU       : {raw_kju_col}")
print(f"Raw file      : {raw_file_col}")
print(f"Raw record ID : {raw_record_col}")
print(f"Raw accession : {raw_acc_col}")
print(f"Master KJU    : {master_kju_col}")


# ============================================================
# UNIQUE FILE PATHS
# ============================================================

file_paths = []

for value in raw[raw_file_col]:

    value = clean(value)

    if value and value not in file_paths:
        file_paths.append(value)


print()
print("SEQUENCE FILE INVENTORY")
print("-" * 110)

print(
    f"Unique sequence file paths : {len(file_paths)}"
)


# ============================================================
# PARSE ACTUAL SEQUENCE FILES
# ============================================================

true_records = []
file_errors = []
missing_files = []

for file_path_string in file_paths:

    file_path = Path(file_path_string)

    if not file_path.exists():

        missing_files.append(
            {
                "Sequence_file": file_path_string,
                "Reason": "FILE_NOT_FOUND",
            }
        )

        continue

    try:

        parsed = parse_sequence_file(
            file_path
        )

        for record in parsed:

            sequence = clean(
                record.get("Sequence", "")
            )

            if not sequence:
                continue

            record_id = clean(
                record.get("Record_ID", "")
            )

            description = clean(
                record.get("Description", "")
            )

            accession = extract_accession(
                record_id,
                description,
            )

            kju_from_record = extract_kju(
                record_id
                + " "
                + description
            )

            sha = sequence_sha256(
                sequence
            )

            true_records.append(
                {
                    "Sequence_file": str(
                        file_path
                    ),
                    "Sequence_file_name": file_path.name,
                    "Record_ID": record_id,
                    "Description": description,
                    "KJU_ID_from_record": kju_from_record,
                    "Sequence_accession_detected": accession,
                    "Normalized_accession": normalize_accession(
                        accession
                    ),
                    "Sequence_length": len(sequence),
                    "Sequence_type": (
                        "16S"
                        if len(sequence) >= 500
                        else "SHORT_OR_OTHER"
                    ),
                    "Sequence": sequence,
                    "Sequence_SHA256": sha,
                }
            )

    except Exception as exc:

        file_errors.append(
            {
                "Sequence_file": str(
                    file_path
                ),
                "Error": repr(exc),
            }
        )


true_df = pd.DataFrame(
    true_records
)


# ============================================================
# MAP RAW KJU TO TRUE RECORDS
# ============================================================

#
# The same FASTA file may contain many records.
#
# STEP 77 knows which KJU generated each raw row.
#
# Therefore we build:
#
#     (file path, record ID)
#
# -> KJU
#
# from STEP77.
#

record_to_kju = {}

for _, row in raw.iterrows():

    kju = clean(
        row[raw_kju_col]
    ).upper()

    file_path = clean(
        row[raw_file_col]
    )

    record_id = ""

    if raw_record_col:
        record_id = clean(
            row[raw_record_col]
        )

    if not kju:
        continue

    if not file_path:
        continue

    key = (
        file_path,
        record_id,
    )

    if record_id:

        record_to_kju[key] = kju


# ============================================================
# ASSIGN KJU
# ============================================================

if not true_df.empty:

    assigned_kjus = []

    for _, row in true_df.iterrows():

        file_path = clean(
            row["Sequence_file"]
        )

        record_id = clean(
            row["Record_ID"]
        )

        description = clean(
            row["Description"]
        )

        kju = extract_kju(
            record_id
            + " "
            + description
        )

        if not kju:

            kju = record_to_kju.get(
                (
                    file_path,
                    record_id,
                ),
                "",
            )

        assigned_kjus.append(
            kju
        )

    true_df["KJU_ID"] = assigned_kjus

else:

    true_df = pd.DataFrame(
        columns=[
            "Sequence_file",
            "Sequence_file_name",
            "Record_ID",
            "Description",
            "KJU_ID_from_record",
            "KJU_ID",
            "Sequence_accession_detected",
            "Normalized_accession",
            "Sequence_length",
            "Sequence_type",
            "Sequence",
            "Sequence_SHA256",
        ]
    )


# ============================================================
# SECONDARY KJU MAPPING
# ============================================================

#
# If record IDs do not contain KJU IDs, use STEP77
# accession/file mapping where possible.
#

raw_map = {}

for _, row in raw.iterrows():

    kju = clean(
        row[raw_kju_col]
    ).upper()

    file_path = clean(
        row[raw_file_col]
    )

    accession = ""

    if raw_acc_col:
        accession = normalize_accession(
            row[raw_acc_col]
        )

    if not kju:
        continue

    if file_path and accession:

        raw_map[
            (
                file_path,
                accession,
            )
        ] = kju


if not true_df.empty:

    for idx in true_df.index:

        if clean(
            true_df.at[idx, "KJU_ID"]
        ):
            continue

        file_path = clean(
            true_df.at[idx, "Sequence_file"]
        )

        accession = normalize_accession(
            true_df.at[
                idx,
                "Sequence_accession_detected",
            ]
        )

        if accession:

            kju = raw_map.get(
                (
                    file_path,
                    accession,
                ),
                "",
            )

            if kju:
                true_df.at[
                    idx,
                    "KJU_ID"
                ] = kju


# ============================================================
# SAVE TRUE RECORDS
# ============================================================

true_df.to_csv(
    TRUE_RECORDS,
    index=False,
)


# ============================================================
# DUPLICATE SEQUENCE ANALYSIS
# ============================================================

if not true_df.empty:

    sha_groups = (
        true_df
        .groupby("Sequence_SHA256")
        .agg(
            Record_count=(
                "Sequence_SHA256",
                "size",
            ),
            KJU_count=(
                "KJU_ID",
                lambda x: len(
                    set(
                        v
                        for v in x
                        if clean(v)
                    )
                ),
            ),
            KJU_IDs=(
                "KJU_ID",
                lambda x: ";".join(
                    sorted(
                        set(
                            v
                            for v in x
                            if clean(v)
                        )
                    )
                ),
            ),
            Accessions=(
                "Normalized_accession",
                lambda x: ";".join(
                    sorted(
                        set(
                            v
                            for v in x
                            if clean(v)
                        )
                    )
                ),
            ),
            Length=(
                "Sequence_length",
                "first",
            ),
        )
        .reset_index()
    )

    duplicate_groups = sha_groups[
        sha_groups["Record_count"] > 1
    ].copy()

else:

    duplicate_groups = pd.DataFrame()


duplicate_groups.to_csv(
    DUPLICATES,
    index=False,
)


# ============================================================
# MASTER KJU SET
# ============================================================

master["KJU_ID"] = (
    master[master_kju_col]
    .map(clean)
    .str.upper()
)

master_kjus = set(
    master["KJU_ID"]
)


# ============================================================
# KJU SEQUENCE AUDIT
# ============================================================

kju_rows = []

for kju in sorted(
    master_kjus
):

    subset = true_df[
        true_df["KJU_ID"]
        == kju
    ]

    sequences = set(
        subset[
            "Sequence_SHA256"
        ].dropna()
    )

    accessions = set(
        a
        for a in subset[
            "Normalized_accession"
        ]
        if clean(a)
    )

    lengths = set(
        subset[
            "Sequence_length"
        ].dropna()
        .astype(int)
    )

    if len(subset) == 0:

        status = "NO_SEQUENCE"

    elif len(sequences) == 1:

        status = "ONE_UNIQUE_SEQUENCE"

    else:

        status = "MULTIPLE_UNIQUE_SEQUENCES"

    kju_rows.append(
        {
            "KJU_ID": kju,
            "Raw_sequence_records": len(subset),
            "Unique_sequences": len(
                sequences
            ),
            "Unique_accessions": len(
                accessions
            ),
            "Unique_lengths": len(
                lengths
            ),
            "Accessions": ";".join(
                sorted(accessions)
            ),
            "Lengths": ";".join(
                str(x)
                for x in sorted(lengths)
            ),
            "Sequence_status": status,
        }
    )


kju_audit = pd.DataFrame(
    kju_rows
)

kju_audit.to_csv(
    KJU_AUDIT,
    index=False,
)


# ============================================================
# ACCESSION AUDIT
# ============================================================

accession_rows = []

for accession, subset in (
    true_df[
        true_df[
            "Normalized_accession"
        ] != ""
    ]
    .groupby(
        "Normalized_accession"
    )
):

    sequences = set(
        subset[
            "Sequence_SHA256"
        ]
    )

    kjus = set(
        k
        for k in subset[
            "KJU_ID"
        ]
        if clean(k)
    )

    lengths = set(
        subset[
            "Sequence_length"
        ]
        .dropna()
        .astype(int)
    )

    accession_rows.append(
        {
            "Normalized_accession": accession,
            "Records": len(subset),
            "KJU_count": len(kjus),
            "KJU_IDs": ";".join(
                sorted(kjus)
            ),
            "Unique_sequences": len(
                sequences
            ),
            "Lengths": ";".join(
                str(x)
                for x in sorted(lengths)
            ),
            "Status": (
                "CONFLICT"
                if len(sequences) > 1
                else "CONSISTENT"
            ),
        }
    )


accession_audit = pd.DataFrame(
    accession_rows
)

accession_audit.to_csv(
    ACCESSION_AUDIT,
    index=False,
)


# ============================================================
# FROZEN MASTER QC
# ============================================================

master_bytes = MASTER.read_bytes()

master_sha256 = hashlib.sha256(
    master_bytes
).hexdigest()


# ============================================================
# STATISTICS
# ============================================================

total_true_records = len(
    true_df
)

kju_with_sequence = (
    kju_audit[
        kju_audit[
            "Raw_sequence_records"
        ]
        > 0
    ]
    .shape[0]
)

kju_without_sequence = (
    kju_audit[
        kju_audit[
            "Raw_sequence_records"
        ]
        == 0
    ]
    .shape[0]
)

one_unique = (
    kju_audit[
        kju_audit[
            "Unique_sequences"
        ]
        == 1
    ]
    .shape[0]
)

multiple_unique = (
    kju_audit[
        kju_audit[
            "Unique_sequences"
        ]
        > 1
    ]
    .shape[0]
)

unresolved_kjus = set(
    master.loc[
        master[
            "Accepted_taxonomy"
        ]
        .map(clean)
        == "",
        "KJU_ID",
    ]
)

unresolved_with_sequence = (
    unresolved_kjus
    & set(
        true_df[
            "KJU_ID"
        ]
        .dropna()
        .map(clean)
    )
)

same_accession_diff_sequence = (
    accession_audit[
        accession_audit[
            "Unique_sequences"
        ]
        > 1
    ]
    .shape[0]
)

different_accession_same_sequence = (
    duplicate_groups[
        duplicate_groups[
            "KJU_count"
        ]
        > 1
    ]
    .shape[0]
)


# ============================================================
# REPORT
# ============================================================

report_lines = []

report_lines.append(
    "=" * 110
)

report_lines.append(
    "STEP 83: TRUE SEQUENCE INVENTORY REPORT"
)

report_lines.append(
    "=" * 110
)

report_lines.append("")

report_lines.append(
    "INPUT"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Frozen master : {MASTER}"
)

report_lines.append(
    f"STEP77 input  : {RAW_INPUT}"
)

report_lines.append("")

report_lines.append(
    "INPUT SIZES"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Master KJU records       : {len(master)}"
)

report_lines.append(
    f"STEP77 raw rows          : {len(raw)}"
)

report_lines.append(
    f"Unique sequence files    : {len(file_paths)}"
)

report_lines.append(
    f"Actual sequence records  : {total_true_records}"
)

report_lines.append("")

report_lines.append(
    "SEQUENCE COVERAGE"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"KJUs with actual sequence    : {kju_with_sequence}"
)

report_lines.append(
    f"KJUs without sequence       : {kju_without_sequence}"
)

report_lines.append(
    f"Exactly one unique sequence : {one_unique}"
)

report_lines.append(
    f"Multiple unique sequences   : {multiple_unique}"
)

report_lines.append("")

report_lines.append(
    "DUPLICATE SEQUENCE QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Exact duplicate SHA groups              : {len(duplicate_groups)}"
)

report_lines.append(
    f"Same accession / different sequence     : {same_accession_diff_sequence}"
)

report_lines.append(
    f"Different KJU / identical sequence      : {different_accession_same_sequence}"
)

report_lines.append("")

report_lines.append(
    "UNRESOLVED QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Unresolved KJUs with actual sequence : {len(unresolved_with_sequence)}"
)

if len(unresolved_with_sequence) == 0:

    report_lines.append(
        "PASS — no unresolved KJU has an actual sequence"
    )

else:

    report_lines.append(
        "WARNING — unresolved KJUs have actual sequences"
    )

    for kju in sorted(
        unresolved_with_sequence
    ):

        report_lines.append(
            f"  {kju}"
        )

report_lines.append("")

report_lines.append(
    "FILE QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Missing sequence files : {len(missing_files)}"
)

report_lines.append(
    f"Files with parse errors: {len(file_errors)}"
)

report_lines.append("")

report_lines.append(
    "FROZEN MASTER QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    "PASS — frozen V1 master was not modified"
)

report_lines.append(
    f"Frozen V1 SHA256: {master_sha256}"
)

report_lines.append("")

report_lines.append(
    "OUTPUT FILES"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"True records     : {TRUE_RECORDS}"
)

report_lines.append(
    f"True inventory   : {TRUE_INVENTORY}"
)

report_lines.append(
    f"Duplicate groups : {DUPLICATES}"
)

report_lines.append(
    f"KJU audit        : {KJU_AUDIT}"
)

report_lines.append(
    f"Accession audit  : {ACCESSION_AUDIT}"
)

report_lines.append(
    f"Missing files    : {MISSING}"
)

report_lines.append(
    f"Parse errors     : {ERRORS}"
)

report_lines.append("")

report_lines.append(
    "=" * 110
)

report_lines.append(
    "STEP 83 COMPLETE"
)

report_lines.append(
    "FROZEN V1 MASTER WAS NOT MODIFIED"
)

report_lines.append(
    "=" * 110
)


REPORT.write_text(
    "\n".join(report_lines),
    encoding="utf-8",
)


# ============================================================
# TRUE INVENTORY
# ============================================================

inventory_rows = []

for kju in sorted(
    master_kjus
):

    subset = true_df[
        true_df["KJU_ID"]
        == kju
    ]

    if len(subset) == 0:

        inventory_rows.append(
            {
                "KJU_ID": kju,
                "Sequence_records": 0,
                "Unique_sequences": 0,
                "Unique_accessions": 0,
                "Sequence_lengths": "",
                "Accessions": "",
                "Status": "NO_SEQUENCE",
            }
        )

        continue

    sequences = set(
        subset[
            "Sequence_SHA256"
        ]
    )

    accessions = set(
        a
        for a in subset[
            "Normalized_accession"
        ]
        if clean(a)
    )

    lengths = set(
        subset[
            "Sequence_length"
        ]
        .dropna()
        .astype(int)
    )

    inventory_rows.append(
        {
            "KJU_ID": kju,
            "Sequence_records": len(
                subset
            ),
            "Unique_sequences": len(
                sequences
            ),
            "Unique_accessions": len(
                accessions
            ),
            "Sequence_lengths": ";".join(
                str(x)
                for x in sorted(lengths)
            ),
            "Accessions": ";".join(
                sorted(accessions)
            ),
            "Status": (
                "ONE_UNIQUE_SEQUENCE"
                if len(sequences) == 1
                else "MULTIPLE_UNIQUE_SEQUENCES"
            ),
        }
    )


pd.DataFrame(
    inventory_rows
).to_csv(
    TRUE_INVENTORY,
    index=False,
)


# ============================================================
# MISSING / ERROR FILES
# ============================================================

pd.DataFrame(
    missing_files,
    columns=[
        "Sequence_file",
        "Reason",
    ],
).to_csv(
    MISSING,
    index=False,
)


pd.DataFrame(
    file_errors,
    columns=[
        "Sequence_file",
        "Error",
    ],
).to_csv(
    ERRORS,
    index=False,
)


# ============================================================
# PRINT REPORT
# ============================================================

print()
print("\n".join(report_lines))

