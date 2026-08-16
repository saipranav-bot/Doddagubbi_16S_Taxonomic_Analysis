#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 79
# SEQUENCE DISCREPANCY INVESTIGATION
#
# INPUT:
#   STEP 76 frozen master
#   STEP 77 raw sequence records
#   STEP 77 sequence inventory
#
# PURPOSE:
#   1. Explain the 101 length mismatches reported by STEP 78
#   2. Determine whether mismatches are caused by:
#        - multiple raw records per KJU
#        - multiple records per accession
#        - different sequence regions
#        - duplicate records
#        - genuinely conflicting sequences
#   3. Determine the number of unique sequences per KJU
#   4. Determine the number of unique sequences per accession
#   5. Identify exact problematic records
#   6. Produce a clean discrepancy table
#
# IMPORTANT:
#   This step DOES NOT modify the frozen V1 master.
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

FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"

MASTER = (
    FREEZE
    / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
)

RAW = (
    ROOT
    / "STEP77_SEQUENCE_QC"
    / "STEP77_RAW_SEQUENCE_RECORDS.csv"
)

INVENTORY = (
    ROOT
    / "STEP77_SEQUENCE_QC"
    / "STEP77_SEQUENCE_INVENTORY.csv"
)

OUT = ROOT / "STEP79_SEQUENCE_DISCREPANCY"
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# OUTPUT FILES
# ============================================================

KJU_AUDIT = OUT / "STEP79_KJU_SEQUENCE_AUDIT.csv"
ACCESSION_AUDIT = OUT / "STEP79_ACCESSION_SEQUENCE_AUDIT.csv"
LENGTH_MISMATCHES = OUT / "STEP79_LENGTH_MISMATCHES.csv"
CONFLICTS = OUT / "STEP79_SEQUENCE_CONFLICTS.csv"
DUPLICATES = OUT / "STEP79_DUPLICATE_SEQUENCE_GROUPS.csv"
REPORT = OUT / "STEP79_SEQUENCE_DISCREPANCY_REPORT.txt"


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

    # Remove whitespace
    value = re.sub(r"\s+", "", value)

    # Remove version suffix
    value = re.sub(r"\.\d+$", "", value)

    return value.upper()


def normalize_sequence(value):
    value = clean(value)

    if not value:
        return ""

    # Remove FASTA whitespace
    value = re.sub(r"\s+", "", value)

    # Remove non-IUPAC nucleotide symbols
    value = re.sub(r"[^ACGTUNRYKMSWBDHV\-]", "", value.upper())

    return value


def sequence_hash(sequence):
    sequence = normalize_sequence(sequence)

    if not sequence:
        return ""

    return hashlib.sha256(
        sequence.encode("utf-8")
    ).hexdigest()


def detect_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None


# ============================================================
# LOAD
# ============================================================

print("=" * 110)
print("STEP 79: SEQUENCE DISCREPANCY INVESTIGATION")
print("=" * 110)

print()
print("INPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"Raw sequences : {RAW}")
print(f"Inventory     : {INVENTORY}")


if not MASTER.exists():
    raise FileNotFoundError(f"Frozen master not found: {MASTER}")

if not RAW.exists():
    raise FileNotFoundError(f"Raw sequence file not found: {RAW}")

if not INVENTORY.exists():
    raise FileNotFoundError(f"Inventory file not found: {INVENTORY}")


master = pd.read_csv(MASTER, dtype=str).fillna("")
raw = pd.read_csv(RAW, dtype=str).fillna("")
inventory = pd.read_csv(INVENTORY, dtype=str).fillna("")


print()
print("INPUT SIZES")
print("-" * 110)
print(f"Master rows        : {len(master)}")
print(f"Raw sequence rows  : {len(raw)}")
print(f"Inventory rows     : {len(inventory)}")


# ============================================================
# DETECT COLUMNS
# ============================================================

master_kju_col = detect_column(
    master,
    ["KJU_ID", "kju_id", "KJU", "kju"]
)

master_acc_col = detect_column(
    master,
    [
        "Accepted_accession",
        "accepted_accession",
        "canonical_accession",
        "accession",
    ]
)

master_length_col = detect_column(
    master,
    [
        "Accepted_length",
        "accepted_length",
        "canonical_length",
        "length",
    ]
)

master_tax_col = detect_column(
    master,
    [
        "Accepted_taxonomy",
        "accepted_taxonomy",
        "canonical_organism",
        "organism",
    ]
)

master_identity_col = detect_column(
    master,
    [
        "Biological_identity",
        "biological_identity",
    ]
)

raw_kju_col = detect_column(
    raw,
    ["KJU_ID", "kju_id", "KJU", "kju"]
)

raw_acc_col = detect_column(
    raw,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "accession",
        "Accession",
    ]
)

raw_length_col = detect_column(
    raw,
    [
        "Sequence_length",
        "sequence_length",
        "length",
        "Length",
    ]
)

raw_seq_col = detect_column(
    raw,
    [
        "Sequence",
        "sequence",
        "Sequence_file",
    ]
)

raw_hash_col = detect_column(
    raw,
    [
        "Sequence_SHA256",
        "sequence_sha256",
        "SHA256",
    ]
)


print()
print("DETECTED COLUMNS")
print("-" * 110)
print(f"Master KJU       : {master_kju_col}")
print(f"Master accession : {master_acc_col}")
print(f"Master length    : {master_length_col}")
print(f"Master taxonomy  : {master_tax_col}")
print(f"Master identity  : {master_identity_col}")
print(f"Raw KJU          : {raw_kju_col}")
print(f"Raw accession    : {raw_acc_col}")
print(f"Raw length       : {raw_length_col}")
print(f"Raw sequence     : {raw_seq_col}")
print(f"Raw SHA256       : {raw_hash_col}")


required = [
    ("master KJU", master_kju_col),
    ("master accession", master_acc_col),
    ("master length", master_length_col),
    ("raw KJU", raw_kju_col),
    ("raw accession", raw_acc_col),
    ("raw length", raw_length_col),
]

missing = [
    name
    for name, col in required
    if col is None
]

if missing:
    raise ValueError(
        "Required columns missing: "
        + ", ".join(missing)
    )


# ============================================================
# NORMALIZE MASTER
# ============================================================

master["_KJU"] = master[master_kju_col].map(clean)
master["_ACCESSION"] = (
    master[master_acc_col]
    .map(normalize_accession)
)

master["_LENGTH"] = pd.to_numeric(
    master[master_length_col]
    .map(clean),
    errors="coerce"
)

master["_TAXONOMY"] = (
    master[master_tax_col].map(clean)
    if master_tax_col
    else ""
)

master["_IDENTITY"] = (
    master[master_identity_col].map(clean)
    if master_identity_col
    else ""
)


# ============================================================
# NORMALIZE RAW
# ============================================================

raw["_KJU"] = raw[raw_kju_col].map(clean)

raw["_ACCESSION"] = (
    raw[raw_acc_col]
    .map(normalize_accession)
)

raw["_LENGTH"] = pd.to_numeric(
    raw[raw_length_col].map(clean),
    errors="coerce"
)


# ============================================================
# SEQUENCE HASH
# ============================================================

if raw_seq_col is not None:

    raw["_SEQUENCE"] = raw[raw_seq_col].map(
        normalize_sequence
    )

    raw["_HASH"] = raw["_SEQUENCE"].map(
        sequence_hash
    )

else:

    raw["_SEQUENCE"] = ""
    raw["_HASH"] = ""

    if raw_hash_col is not None:
        raw["_HASH"] = raw[raw_hash_col].map(clean)


# ============================================================
# ONLY REAL SEQUENCE RECORDS
# ============================================================

sequence_rows = raw[
    (raw["_KJU"] != "")
    &
    (
        (raw["_HASH"] != "")
        |
        (raw["_LENGTH"].notna())
    )
].copy()


print()
print("SEQUENCE DATA")
print("-" * 110)
print(f"Sequence-bearing rows : {len(sequence_rows)}")
print(
    f"KJUs represented      : "
    f"{sequence_rows['_KJU'].nunique()}"
)


# ============================================================
# MASTER LOOKUP
# ============================================================

master_lookup = (
    master
    .set_index("_KJU")
    .to_dict("index")
)


# ============================================================
# STEP 79A
# PER-KJU SEQUENCE AUDIT
# ============================================================

kju_rows = []

for kju in sorted(master["_KJU"].unique()):

    m = master_lookup.get(kju, {})

    accepted_acc = clean(
        m.get("_ACCESSION", "")
    )

    accepted_length = m.get("_LENGTH")

    accepted_tax = clean(
        m.get("_TAXONOMY", "")
    )

    biological_identity = clean(
        m.get("_IDENTITY", "")
    )

    r = sequence_rows[
        sequence_rows["_KJU"] == kju
    ].copy()

    raw_records = len(r)

    unique_accessions = sorted(
        {
            x
            for x in r["_ACCESSION"]
            if clean(x)
        }
    )

    unique_hashes = sorted(
        {
            x
            for x in r["_HASH"]
            if clean(x)
        }
    )

    lengths = sorted(
        {
            int(x)
            for x in r["_LENGTH"]
            if pd.notna(x)
        }
    )

    accession_match = (
        accepted_acc != ""
        and any(
            x == accepted_acc
            for x in unique_accessions
        )
    )

    length_match = (
        pd.notna(accepted_length)
        and int(accepted_length) in lengths
    )

    if raw_records == 0:
        status = "NO_SEQUENCE"

    elif len(unique_hashes) == 1:
        status = "ONE_UNIQUE_SEQUENCE"

    elif len(unique_hashes) > 1:
        status = "MULTIPLE_UNIQUE_SEQUENCES"

    else:
        status = "ACCESSION_OR_LENGTH_ONLY"

    kju_rows.append(
        {
            "KJU_ID": kju,
            "Accepted_accession": accepted_acc,
            "Accepted_length": (
                int(accepted_length)
                if pd.notna(accepted_length)
                else ""
            ),
            "Accepted_taxonomy": accepted_tax,
            "Biological_identity": biological_identity,
            "Raw_record_count": raw_records,
            "Unique_accession_count": len(
                unique_accessions
            ),
            "Unique_accessions": ";".join(
                unique_accessions
            ),
            "Unique_sequence_count": len(
                unique_hashes
            ),
            "Sequence_hashes": ";".join(
                unique_hashes
            ),
            "Observed_lengths": ";".join(
                map(str, lengths)
            ),
            "Accession_match": (
                "YES"
                if accession_match
                else "NO"
            ),
            "Length_match": (
                "YES"
                if length_match
                else "NO"
            ),
            "Sequence_status": status,
        }
    )


kju_audit = pd.DataFrame(kju_rows)

kju_audit.to_csv(
    KJU_AUDIT,
    index=False
)


# ============================================================
# STEP 79B
# LENGTH MISMATCH ANALYSIS
# ============================================================

length_problem_rows = []

for _, row in kju_audit.iterrows():

    if (
        row["Accepted_length"] != ""
        and row["Length_match"] == "NO"
        and int(row["Raw_record_count"]) > 0
    ):

        length_problem_rows.append(
            row.to_dict()
        )


length_problems = pd.DataFrame(
    length_problem_rows
)

if len(length_problems) > 0:
    length_problems.to_csv(
        LENGTH_MISMATCHES,
        index=False
    )
else:
    pd.DataFrame(
        columns=kju_audit.columns
    ).to_csv(
        LENGTH_MISMATCHES,
        index=False
    )


# ============================================================
# STEP 79C
# ACCESSION-LEVEL AUDIT
# ============================================================

accession_rows = []

for accession in sorted(
    {
        x
        for x in sequence_rows["_ACCESSION"]
        if clean(x)
    }
):

    r = sequence_rows[
        sequence_rows["_ACCESSION"] == accession
    ].copy()

    kjus = sorted(
        {
            x
            for x in r["_KJU"]
            if clean(x)
        }
    )

    hashes = sorted(
        {
            x
            for x in r["_HASH"]
            if clean(x)
        }
    )

    lengths = sorted(
        {
            int(x)
            for x in r["_LENGTH"]
            if pd.notna(x)
        }
    )

    if len(kjus) > 1 and len(hashes) > 1:
        status = "CONFLICTING_ACCESSION_AND_SEQUENCE"

    elif len(kjus) > 1:
        status = "ACCESSION_SHARED_ACROSS_KJUS"

    elif len(hashes) > 1:
        status = "SAME_ACCESSION_MULTIPLE_SEQUENCES"

    else:
        status = "CONSISTENT"

    accession_rows.append(
        {
            "Normalized_accession": accession,
            "KJU_count": len(kjus),
            "KJUs": ";".join(kjus),
            "Unique_sequence_count": len(hashes),
            "Sequence_hashes": ";".join(hashes),
            "Observed_lengths": ";".join(
                map(str, lengths)
            ),
            "Raw_record_count": len(r),
            "Status": status,
        }
    )


accession_audit = pd.DataFrame(
    accession_rows
)

accession_audit.to_csv(
    ACCESSION_AUDIT,
    index=False
)


# ============================================================
# STEP 79D
# TRUE SEQUENCE CONFLICTS
# ============================================================

conflict_rows = []

for _, row in accession_audit.iterrows():

    if row["Status"] in {
        "CONFLICTING_ACCESSION_AND_SEQUENCE",
        "SAME_ACCESSION_MULTIPLE_SEQUENCES",
    }:

        conflict_rows.append(row.to_dict())


conflicts = pd.DataFrame(
    conflict_rows
)

if len(conflicts) > 0:
    conflicts.to_csv(
        CONFLICTS,
        index=False
    )
else:
    pd.DataFrame(
        columns=accession_audit.columns
    ).to_csv(
        CONFLICTS,
        index=False
    )


# ============================================================
# STEP 79E
# DUPLICATE SEQUENCE GROUPS
# ============================================================

duplicate_rows = []

for seq_hash, r in sequence_rows[
    sequence_rows["_HASH"] != ""
].groupby("_HASH"):

    kjus = sorted(
        set(r["_KJU"])
    )

    accessions = sorted(
        {
            x
            for x in r["_ACCESSION"]
            if clean(x)
        }
    )

    lengths = sorted(
        {
            int(x)
            for x in r["_LENGTH"]
            if pd.notna(x)
        }
    )

    if len(kjus) > 1 or len(accessions) > 1:

        if len(kjus) == 1:
            duplicate_type = (
                "SAME_KJU_MULTIPLE_RECORDS"
            )

        elif len(accessions) == 1:
            duplicate_type = (
                "SAME_SEQUENCE_DIFFERENT_KJU"
            )

        else:
            duplicate_type = (
                "DIFFERENT_ACCESSION_IDENTICAL_SEQUENCE"
            )

        duplicate_rows.append(
            {
                "Sequence_SHA256": seq_hash,
                "KJU_count": len(kjus),
                "KJUs": ";".join(kjus),
                "Accession_count": len(accessions),
                "Accessions": ";".join(accessions),
                "Observed_lengths": ";".join(
                    map(str, lengths)
                ),
                "Raw_record_count": len(r),
                "Duplicate_type": duplicate_type,
            }
        )


duplicates = pd.DataFrame(
    duplicate_rows
)

duplicates.to_csv(
    DUPLICATES,
    index=False
)


# ============================================================
# SUMMARY COUNTS
# ============================================================

total_kjus = len(master)

resolved_mask = master["_TAXONOMY"] != ""

resolved_kjus = int(
    resolved_mask.sum()
)

unresolved_kjus = int(
    (~resolved_mask).sum()
)

with_sequence = int(
    (kju_audit["Raw_record_count"] > 0).sum()
)

without_sequence = total_kjus - with_sequence

resolved_with_sequence = int(
    (
        (kju_audit["Raw_record_count"] > 0)
        &
        (kju_audit["Accepted_taxonomy"] != "")
    ).sum()
)

unresolved_with_sequence = int(
    (
        (kju_audit["Raw_record_count"] > 0)
        &
        (kju_audit["Accepted_taxonomy"] == "")
    ).sum()
)

multiple_sequence_kjus = int(
    (
        kju_audit["Unique_sequence_count"]
        .astype(int)
        > 1
    ).sum()
)

single_sequence_kjus = int(
    (
        kju_audit["Unique_sequence_count"]
        .astype(int)
        == 1
    ).sum()
)

length_mismatch_count = len(
    length_problems
)

accession_conflict_count = len(
    conflicts
)

duplicate_group_count = len(
    duplicates
)


# ============================================================
# REPORT
# ============================================================

lines = []

lines.append("=" * 110)
lines.append(
    "STEP 79: SEQUENCE DISCREPANCY INVESTIGATION"
)
lines.append("=" * 110)

lines.append("")
lines.append("MASTER")
lines.append("-" * 110)
lines.append(
    f"Total KJU records             : {total_kjus}"
)
lines.append(
    f"Resolved records              : {resolved_kjus}"
)
lines.append(
    f"Unresolved records            : {unresolved_kjus}"
)

lines.append("")
lines.append("SEQUENCE COVERAGE")
lines.append("-" * 110)
lines.append(
    f"KJUs with sequence            : {with_sequence}"
)
lines.append(
    f"KJUs without sequence         : {without_sequence}"
)
lines.append(
    f"Resolved + sequence           : {resolved_with_sequence}"
)
lines.append(
    f"Unresolved + sequence         : {unresolved_with_sequence}"
)

lines.append("")
lines.append("PER-KJU SEQUENCE STRUCTURE")
lines.append("-" * 110)
lines.append(
    f"KJUs with exactly one sequence: {single_sequence_kjus}"
)
lines.append(
    f"KJUs with multiple sequences  : {multiple_sequence_kjus}"
)

lines.append("")
lines.append("LENGTH AUDIT")
lines.append("-" * 110)
lines.append(
    f"Length mismatch KJUs          : {length_mismatch_count}"
)

lines.append("")
lines.append("ACCESSION AUDIT")
lines.append("-" * 110)
lines.append(
    f"Accession conflict groups     : {accession_conflict_count}"
)

lines.append("")
lines.append("DUPLICATE SEQUENCE AUDIT")
lines.append("-" * 110)
lines.append(
    f"Duplicate sequence groups     : {duplicate_group_count}"
)

lines.append("")
lines.append("IMPORTANT INTERPRETATION")
lines.append("-" * 110)

if unresolved_with_sequence == 0:
    lines.append(
        "PASS — no unresolved KJU has a sequence."
    )
else:
    lines.append(
        "WARNING — unresolved KJUs with sequence detected."
    )

if length_mismatch_count == 0:
    lines.append(
        "PASS — no KJU-level length discrepancy."
    )
else:
    lines.append(
        "REVIEW — length discrepancies require record-level inspection."
    )

if accession_conflict_count == 0:
    lines.append(
        "PASS — no accession-level sequence conflicts."
    )
else:
    lines.append(
        "REVIEW — accession-level sequence conflicts detected."
    )

lines.append("")
lines.append("OUTPUT FILES")
lines.append("-" * 110)
lines.append(f"KJU audit       : {KJU_AUDIT}")
lines.append(f"Accession audit : {ACCESSION_AUDIT}")
lines.append(f"Length problems : {LENGTH_MISMATCHES}")
lines.append(f"Conflicts       : {CONFLICTS}")
lines.append(f"Duplicates      : {DUPLICATES}")
lines.append(f"QC report       : {REPORT}")

lines.append("")
lines.append("=" * 110)
lines.append(
    "FROZEN V1 MASTER WAS NOT MODIFIED"
)
lines.append("=" * 110)


REPORT.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8"
)


# ============================================================
# CONSOLE
# ============================================================

print()
print("\n".join(lines))


# ============================================================
# FINAL SAFETY CHECK
# ============================================================

if len(master) != 204:
    raise ValueError(
        "Frozen master row count changed unexpectedly."
    )

if master["_KJU"].nunique() != 204:
    raise ValueError(
        "Frozen master KJU uniqueness check failed."
    )

print()
print("=" * 110)
print("STEP 79 COMPLETE")
print("=" * 110)
print()
print(
    "NO FILE IN FINAL_KJU_MASTER_FREEZE_V1 WAS MODIFIED."
)
print()
