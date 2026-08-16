#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 82
# SEQUENCE FIELD INTEGRITY DIAGNOSTIC
#
# PURPOSE:
#   Diagnose why STEP 81 reported:
#
#       2320 raw rows
#       2320 unique candidates
#
# despite STEP 77 reporting:
#
#       1946 exact duplicate rows
#
# This step determines:
#
#   1. What Sequence_file actually contains
#   2. Whether Sequence_SHA256 is valid
#   3. Whether duplicate hashes exist
#   4. Whether Sequence_file contains:
#        - actual nucleotide sequence
#        - file path
#        - filename
#        - FASTA text
#        - another identifier
#   5. Whether sequence normalization is destroying
#      duplicate detection
#
# IMPORTANT:
#   READ-ONLY.
#   DOES NOT MODIFY THE FROZEN MASTER.
# ============================================================


from pathlib import Path
import pandas as pd
import hashlib
import re


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

RAW = STEP77 / "STEP77_RAW_SEQUENCE_RECORDS.csv"

OUT = ROOT / "STEP82_SEQUENCE_FIELD_DIAGNOSTIC"

OUT.mkdir(parents=True, exist_ok=True)

REPORT = OUT / "STEP82_SEQUENCE_FIELD_DIAGNOSTIC_REPORT.txt"
SAMPLE = OUT / "STEP82_SEQUENCE_FIELD_SAMPLE.csv"
HASH_AUDIT = OUT / "STEP82_HASH_AUDIT.csv"
KJU_AUDIT = OUT / "STEP82_KJU_DUPLICATE_AUDIT.csv"


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

    return str(value).strip()


def normalize_sequence(value):

    value = clean(value)

    if not value:
        return ""

    # FASTA text
    if ">" in value:

        lines = value.splitlines()

        seq_lines = [
            x.strip()
            for x in lines
            if x.strip()
            and not x.strip().startswith(">")
        ]

        value = "".join(seq_lines)

    # Remove whitespace
    value = re.sub(r"\s+", "", value)

    return value.upper()


def sha256_text(value):

    return hashlib.sha256(
        value.encode("utf-8")
    ).hexdigest()


def is_dna_like(value):

    value = normalize_sequence(value)

    if not value:
        return False

    return bool(
        re.fullmatch(
            r"[ACGTNRYKMSWBDHVU\-]+",
            value
        )
    )


def classify_field(value):

    value = clean(value)

    if not value:
        return "BLANK"

    if value.startswith(">"):
        return "FASTA_TEXT"

    if is_dna_like(value):
        return "NUCLEOTIDE_SEQUENCE"

    lower = value.lower()

    if (
        "/" in value
        or "\\" in value
        or lower.endswith(".fasta")
        or lower.endswith(".fa")
        or lower.endswith(".fas")
        or lower.endswith(".gb")
        or lower.endswith(".gbk")
        or lower.endswith(".genbank")
        or lower.endswith(".txt")
    ):
        return "FILE_PATH_OR_FILENAME"

    if len(value) < 100:
        return "IDENTIFIER_OR_METADATA"

    return "OTHER_TEXT"


# ============================================================
# START
# ============================================================

print("=" * 110)
print("STEP 82: SEQUENCE FIELD INTEGRITY DIAGNOSTIC")
print("=" * 110)


if not RAW.exists():

    raise FileNotFoundError(
        f"Missing Step 77 raw sequence file:\n{RAW}"
    )


raw = pd.read_csv(
    RAW,
    dtype=str,
).fillna("")


print("\nINPUT")
print("-" * 110)
print(f"Raw sequence file : {RAW}")
print(f"Rows              : {len(raw)}")
print(f"Columns           : {list(raw.columns)}")


# ============================================================
# REQUIRED COLUMNS
# ============================================================

def detect_column(df, candidates):

    lower = {
        str(c).lower(): c
        for c in df.columns
    }

    for candidate in candidates:

        if candidate.lower() in lower:

            return lower[candidate.lower()]

    return None


kju_col = detect_column(
    raw,
    [
        "KJU_ID",
        "kju_id",
        "KJU",
    ]
)

seq_col = detect_column(
    raw,
    [
        "Sequence_file",
        "sequence",
        "Sequence",
        "sequence_data",
    ]
)

sha_col = detect_column(
    raw,
    [
        "Sequence_SHA256",
        "sequence_sha256",
        "SHA256",
        "sha256",
    ]
)

acc_col = detect_column(
    raw,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "accession",
    ]
)

length_col = detect_column(
    raw,
    [
        "Sequence_length",
        "sequence_length",
        "length",
    ]
)


print("\nDETECTED COLUMNS")
print("-" * 110)
print(f"KJU column       : {kju_col}")
print(f"Sequence column  : {seq_col}")
print(f"SHA256 column    : {sha_col}")
print(f"Accession column : {acc_col}")
print(f"Length column    : {length_col}")


if kju_col is None:
    raise ValueError("KJU column could not be detected.")

if seq_col is None:
    raise ValueError(
        "Sequence column could not be detected."
    )


# ============================================================
# FIELD CLASSIFICATION
# ============================================================

raw["_KJU"] = raw[kju_col].map(clean)
raw["_SEQ_FIELD"] = raw[seq_col].map(clean)

raw["_FIELD_TYPE"] = raw["_SEQ_FIELD"].map(
    classify_field
)


print("\nSEQUENCE FIELD TYPE")
print("-" * 110)

field_counts = (
    raw["_FIELD_TYPE"]
    .value_counts()
    .sort_index()
)

for field_type, count in field_counts.items():

    print(
        f"{field_type:<30} : {count}"
    )


# ============================================================
# SAMPLE VALUES
# ============================================================

sample_cols = [
    c
    for c in [
        kju_col,
        seq_col,
        sha_col,
        acc_col,
        length_col,
    ]
    if c is not None
]

sample = raw[
    sample_cols
].head(30).copy()

sample.to_csv(
    SAMPLE,
    index=False,
)


print("\nFIRST 10 SEQUENCE-FIELD VALUES")
print("-" * 110)

for i, value in enumerate(
    raw["_SEQ_FIELD"].head(10),
    start=1
):

    preview = value

    if len(preview) > 150:
        preview = preview[:150] + "..."

    print(
        f"{i:02d}: {preview}"
    )


# ============================================================
# SHA256 AUDIT
# ============================================================

if sha_col is not None:

    raw["_REPORTED_SHA"] = (
        raw[sha_col]
        .map(clean)
        .str.lower()
    )

else:

    raw["_REPORTED_SHA"] = ""


# Recalculate hash from actual field
raw["_FIELD_NORMALIZED"] = raw[
    "_SEQ_FIELD"
].map(normalize_sequence)

raw["_CALCULATED_SHA"] = raw[
    "_FIELD_NORMALIZED"
].map(
    lambda x:
        sha256_text(x)
        if x
        else ""
)


raw["_SHA_MATCH"] = (
    (raw["_REPORTED_SHA"] != "")
    &
    (
        raw["_REPORTED_SHA"]
        == raw["_CALCULATED_SHA"]
    )
)


print("\nHASH QC")
print("-" * 110)

reported_hashes = (
    raw.loc[
        raw["_REPORTED_SHA"] != "",
        "_REPORTED_SHA"
    ]
)

calculated_hashes = (
    raw.loc[
        raw["_CALCULATED_SHA"] != "",
        "_CALCULATED_SHA"
    ]
)

print(
    f"Rows with reported SHA256 : "
    f"{len(reported_hashes)}"
)

print(
    f"Unique reported SHA256    : "
    f"{reported_hashes.nunique()}"
)

print(
    f"Rows with calculated SHA   : "
    f"{len(calculated_hashes)}"
)

print(
    f"Unique calculated SHA      : "
    f"{calculated_hashes.nunique()}"
)

print(
    f"Reported/calculated match  : "
    f"{int(raw['_SHA_MATCH'].sum())}"
)


# ============================================================
# DUPLICATE HASH ANALYSIS
# ============================================================

hash_counts = (
    raw.loc[
        raw["_REPORTED_SHA"] != "",
        "_REPORTED_SHA"
    ]
    .value_counts()
)

duplicate_hashes = (
    hash_counts[
        hash_counts > 1
    ]
)


print("\nREPORTED SHA256 DUPLICATES")
print("-" * 110)

print(
    f"Duplicate SHA groups : "
    f"{len(duplicate_hashes)}"
)

print(
    f"Rows participating   : "
    f"{int(duplicate_hashes.sum())}"
)


# ============================================================
# CALCULATED HASH DUPLICATES
# ============================================================

calc_hash_counts = (
    raw.loc[
        raw["_CALCULATED_SHA"] != "",
        "_CALCULATED_SHA"
    ]
    .value_counts()
)

calc_duplicate_hashes = (
    calc_hash_counts[
        calc_hash_counts > 1
    ]
)


print("\nCALCULATED HASH DUPLICATES")
print("-" * 110)

print(
    f"Duplicate calculated groups : "
    f"{len(calc_duplicate_hashes)}"
)

print(
    f"Rows participating          : "
    f"{int(calc_duplicate_hashes.sum())}"
)


# ============================================================
# KJU DUPLICATE ANALYSIS
# ============================================================

kju_rows = []

for kju, group in raw.groupby(
    "_KJU",
    dropna=False,
):

    if not kju:
        continue

    reported_unique = group[
        group["_REPORTED_SHA"] != ""
    ]["_REPORTED_SHA"].nunique()

    calculated_unique = group[
        group["_CALCULATED_SHA"] != ""
    ]["_CALCULATED_SHA"].nunique()

    field_unique = group[
        group["_FIELD_NORMALIZED"] != ""
    ]["_FIELD_NORMALIZED"].nunique()

    kju_rows.append({
        "KJU_ID": kju,
        "Raw_rows": len(group),
        "Reported_SHA_unique":
            reported_unique,
        "Calculated_SHA_unique":
            calculated_unique,
        "Normalized_field_unique":
            field_unique,
    })


kju_audit = pd.DataFrame(
    kju_rows
)

kju_audit.to_csv(
    KJU_AUDIT,
    index=False,
)


# ============================================================
# CRITICAL DIAGNOSTIC
# ============================================================

print("\nCRITICAL INTERPRETATION")
print("-" * 110)

sequence_rows = int(
    (
        raw["_FIELD_TYPE"]
        == "NUCLEOTIDE_SEQUENCE"
    ).sum()
)

path_rows = int(
    (
        raw["_FIELD_TYPE"]
        == "FILE_PATH_OR_FILENAME"
    ).sum()
)

fasta_rows = int(
    (
        raw["_FIELD_TYPE"]
        == "FASTA_TEXT"
    ).sum()
)


if sequence_rows > 0:

    print(
        f"Detected nucleotide-like rows : "
        f"{sequence_rows}"
    )

else:

    print(
        "WARNING — Sequence_file does not "
        "appear to contain nucleotide sequences."
    )


if path_rows > 0:

    print(
        f"Detected file/path rows        : "
        f"{path_rows}"
    )


if fasta_rows > 0:

    print(
        f"Detected FASTA-text rows       : "
        f"{fasta_rows}"
    )


# ============================================================
# REPORT
# ============================================================

with open(
    REPORT,
    "w",
    encoding="utf-8",
) as fh:

    fh.write(
        "=" * 110 + "\n"
    )

    fh.write(
        "STEP 82: SEQUENCE FIELD "
        "INTEGRITY DIAGNOSTIC\n"
    )

    fh.write(
        "=" * 110 + "\n\n"
    )

    fh.write(
        f"Input rows: {len(raw)}\n\n"
    )

    fh.write(
        "DETECTED COLUMNS\n"
    )

    fh.write(
        f"KJU: {kju_col}\n"
    )

    fh.write(
        f"Sequence field: {seq_col}\n"
    )

    fh.write(
        f"SHA256: {sha_col}\n"
    )

    fh.write(
        f"Accession: {acc_col}\n"
    )

    fh.write(
        f"Length: {length_col}\n\n"
    )

    fh.write(
        "FIELD TYPES\n"
    )

    for field_type, count in field_counts.items():

        fh.write(
            f"{field_type}: {count}\n"
        )

    fh.write("\n")

    fh.write(
        "HASH QC\n"
    )

    fh.write(
        f"Reported hashes: "
        f"{len(reported_hashes)}\n"
    )

    fh.write(
        f"Unique reported hashes: "
        f"{reported_hashes.nunique()}\n"
    )

    fh.write(
        f"Calculated hashes: "
        f"{len(calculated_hashes)}\n"
    )

    fh.write(
        f"Unique calculated hashes: "
        f"{calculated_hashes.nunique()}\n"
    )

    fh.write(
        f"Hash matches: "
        f"{int(raw['_SHA_MATCH'].sum())}\n"
    )

    fh.write("\n")

    fh.write(
        "DUPLICATE HASHES\n"
    )

    fh.write(
        f"Reported duplicate groups: "
        f"{len(duplicate_hashes)}\n"
    )

    fh.write(
        f"Calculated duplicate groups: "
        f"{len(calc_duplicate_hashes)}\n"
    )

    fh.write("\n")

    fh.write(
        "IMPORTANT:\n"
    )

    fh.write(
        "This diagnostic does not modify the "
        "frozen V1 master.\n"
    )


# ============================================================
# FINAL
# ============================================================

print("\nOUTPUT FILES")
print("-" * 110)
print(f"Report       : {REPORT}")
print(f"Sample       : {SAMPLE}")
print(f"Hash audit   : {HASH_AUDIT}")
print(f"KJU audit    : {KJU_AUDIT}")

print("\nFROZEN MASTER")
print("-" * 110)
print(
    "PASS — frozen V1 master was not modified."
)

print("\n" + "=" * 110)
print("STEP 82 COMPLETE")
print("=" * 110)
