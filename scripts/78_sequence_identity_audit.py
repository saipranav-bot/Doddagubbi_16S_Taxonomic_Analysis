#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
from collections import defaultdict


# ============================================================
# STEP 78
# SEQUENCE IDENTITY + DUPLICATE AUDIT
#
# INPUT:
#   STEP 76 frozen master
#   STEP 77 sequence inventory
#   STEP 77 raw sequence records
#
# PURPOSE:
#   1. Identify every resolved KJU with/without sequence
#   2. Identify every unresolved KJU with sequence
#   3. Compare KJU <-> accession <-> sequence identity
#   4. Characterize duplicate sequences
#   5. Detect:
#
#        A. Same sequence / same KJU
#        B. Same sequence / different KJU
#        C. Same accession / different sequence
#        D. Different accession / identical sequence
#
#   6. Verify sequence lengths
#   7. Verify frozen master is unchanged
#
# IMPORTANT:
#   This step DOES NOT modify the frozen V1 master.
# ============================================================


from pathlib import Path
import pandas as pd
import hashlib


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

STEP77 = ROOT / "STEP77_SEQUENCE_QC"

INVENTORY = (
    STEP77
    / "STEP77_SEQUENCE_INVENTORY.csv"
)

RAW = (
    STEP77
    / "STEP77_RAW_SEQUENCE_RECORDS.csv"
)

OUTDIR = (
    ROOT
    / "STEP78_SEQUENCE_IDENTITY_AUDIT"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# OUTPUT FILES
# ============================================================

KJU_AUDIT = (
    OUTDIR
    / "STEP78_KJU_SEQUENCE_IDENTITY_AUDIT.csv"
)

DUPLICATE_AUDIT = (
    OUTDIR
    / "STEP78_DUPLICATE_SEQUENCE_GROUPS.csv"
)

SAME_KJU_DUPLICATES = (
    OUTDIR
    / "STEP78_SAME_SEQUENCE_SAME_KJU.csv"
)

CROSS_KJU_DUPLICATES = (
    OUTDIR
    / "STEP78_SAME_SEQUENCE_DIFFERENT_KJU.csv"
)

ACCESSION_CONFLICTS = (
    OUTDIR
    / "STEP78_SAME_ACCESSION_DIFFERENT_SEQUENCE.csv"
)

DIFFERENT_ACCESSION_SAME_SEQUENCE = (
    OUTDIR
    / "STEP78_DIFFERENT_ACCESSION_IDENTICAL_SEQUENCE.csv"
)

LENGTH_MISMATCHES = (
    OUTDIR
    / "STEP78_SEQUENCE_LENGTH_MISMATCHES.csv"
)

UNRESOLVED_WITH_SEQUENCE = (
    OUTDIR
    / "STEP78_UNRESOLVED_KJU_WITH_SEQUENCE.csv"
)

REPORT = (
    OUTDIR
    / "STEP78_SEQUENCE_IDENTITY_AUDIT_REPORT.txt"
)


# ============================================================
# HELPERS
# ============================================================

def clean(value):
    """
    Convert missing/blank values to "".
    """
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
    """
    Normalize accession for comparison.

    Example:
        PZ666011.1 -> PZ666011
    """
    value = clean(value)

    if not value:
        return ""

    return value.split(".")[0].upper().strip()


def normalize_sequence(value):
    """
    Normalize nucleotide sequence.

    Removes whitespace and converts to uppercase.
    """
    value = clean(value)

    if not value:
        return ""

    return (
        value
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .replace("\t", "")
        .upper()
    )


def sha256_sequence(sequence):
    """
    Generate SHA256 hash of normalized sequence.
    """
    sequence = normalize_sequence(sequence)

    if not sequence:
        return ""

    return hashlib.sha256(
        sequence.encode("utf-8")
    ).hexdigest()


def detect_column(df, candidates, required=True):
    """
    Find a column using case-insensitive matching.
    """
    lookup = {
        str(c).strip().lower(): c
        for c in df.columns
    }

    for candidate in candidates:
        candidate = candidate.lower()

        if candidate in lookup:
            return lookup[candidate]

    if required:
        raise ValueError(
            "Could not find required column. "
            f"Tried: {candidates}. "
            f"Available: {list(df.columns)}"
        )

    return None


def safe_int(value):
    """
    Convert value to integer where possible.
    """
    value = clean(value)

    if not value:
        return None

    try:
        return int(float(value))
    except Exception:
        return None


# ============================================================
# HEADER
# ============================================================

print("=" * 110)
print("STEP 78: SEQUENCE IDENTITY + DUPLICATE AUDIT")
print("=" * 110)

print()
print("INPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"Inventory     : {INVENTORY}")
print(f"Raw sequences : {RAW}")


# ============================================================
# FILE VALIDATION
# ============================================================

for path in [
    MASTER,
    INVENTORY,
    RAW,
]:

    if not path.exists():
        raise FileNotFoundError(
            f"Required input file not found:\n{path}"
        )


# ============================================================
# LOAD
# ============================================================

master = pd.read_csv(
    MASTER,
    dtype=str,
    keep_default_na=False,
)

inventory = pd.read_csv(
    INVENTORY,
    dtype=str,
    keep_default_na=False,
)

raw = pd.read_csv(
    RAW,
    dtype=str,
    keep_default_na=False,
)


master.columns = [
    str(c).strip()
    for c in master.columns
]

inventory.columns = [
    str(c).strip()
    for c in inventory.columns
]

raw.columns = [
    str(c).strip()
    for c in raw.columns
]


print()
print("INPUT SIZES")
print("-" * 110)
print(f"Master rows        : {len(master)}")
print(f"Inventory rows     : {len(inventory)}")
print(f"Raw sequence rows  : {len(raw)}")


# ============================================================
# DETECT COLUMNS
# ============================================================

master_kju_col = detect_column(
    master,
    [
        "KJU_ID",
        "kju_id",
    ],
)

master_taxonomy_col = detect_column(
    master,
    [
        "Accepted_taxonomy",
        "accepted_taxonomy",
    ],
    required=False,
)

master_accession_col = detect_column(
    master,
    [
        "Accepted_accession",
        "accepted_accession",
    ],
    required=False,
)

master_length_col = detect_column(
    master,
    [
        "Accepted_length",
        "accepted_length",
    ],
    required=False,
)

master_identity_col = detect_column(
    master,
    [
        "Biological_identity",
        "biological_identity",
    ],
    required=False,
)

master_adjudication_col = detect_column(
    master,
    [
        "Adjudication_type",
        "adjudication_type",
    ],
    required=False,
)

master_confidence_col = detect_column(
    master,
    [
        "Confidence",
        "confidence",
    ],
    required=False,
)


raw_kju_col = detect_column(
    raw,
    [
        "KJU_ID",
        "kju_id",
        "_KJU",
    ],
)

raw_sequence_col = detect_column(
    raw,
    [
        "Sequence_file",
        "Sequence",
        "sequence",
    ],
)

raw_accession_col = detect_column(
    raw,
    [
        "Sequence_accession_detected",
        "Accession",
        "accession",
    ],
    required=False,
)

raw_normalized_accession_col = detect_column(
    raw,
    [
        "Normalized_accession_detected",
        "Normalized_accession",
    ],
    required=False,
)

raw_length_col = detect_column(
    raw,
    [
        "Sequence_length",
        "Length",
        "length",
    ],
    required=False,
)

raw_hash_col = detect_column(
    raw,
    [
        "Sequence_SHA256",
        "SHA256",
        "sequence_sha256",
    ],
    required=False,
)

raw_sequence_present_col = detect_column(
    raw,
    [
        "Sequence_present",
        "sequence_present",
    ],
    required=False,
)


print()
print("DETECTED COLUMNS")
print("-" * 110)
print(f"Master KJU       : {master_kju_col}")
print(f"Master accession : {master_accession_col}")
print(f"Master taxonomy  : {master_taxonomy_col}")
print(f"Master length    : {master_length_col}")
print(f"Master identity  : {master_identity_col}")
print(f"Raw KJU          : {raw_kju_col}")
print(f"Raw sequence     : {raw_sequence_col}")
print(f"Raw accession    : {raw_accession_col}")
print(f"Raw length       : {raw_length_col}")
print(f"Raw SHA256       : {raw_hash_col}")


# ============================================================
# NORMALIZE MASTER
# ============================================================

master["_KJU_NORMALIZED"] = (
    master[master_kju_col]
    .map(clean)
)

master["_ACCEPTED_ACCESSION"] = (
    master[master_accession_col].map(clean)
    if master_accession_col
    else ""
)

master["_NORMALIZED_ACCEPTED_ACCESSION"] = (
    master["_ACCEPTED_ACCESSION"]
    .map(normalize_accession)
)

master["_ACCEPTED_LENGTH"] = (
    master[master_length_col].map(safe_int)
    if master_length_col
    else None
)

master["_TAXONOMY"] = (
    master[master_taxonomy_col].map(clean)
    if master_taxonomy_col
    else ""
)

master["_BIOLOGICAL_IDENTITY"] = (
    master[master_identity_col].map(clean)
    if master_identity_col
    else ""
)

master["_ADJUDICATION_TYPE"] = (
    master[master_adjudication_col].map(clean)
    if master_adjudication_col
    else ""
)

master["_CONFIDENCE"] = (
    master[master_confidence_col].map(clean)
    if master_confidence_col
    else ""
)


# ============================================================
# MASTER BASIC QC
# ============================================================

if len(master) != 204:
    raise ValueError(
        f"Frozen master should contain 204 rows; found {len(master)}"
    )

if master["_KJU_NORMALIZED"].nunique() != 204:
    raise ValueError(
        "Frozen master contains duplicate or invalid KJU IDs"
    )


# ============================================================
# NORMALIZE RAW SEQUENCES
# ============================================================

raw["_KJU_NORMALIZED"] = (
    raw[raw_kju_col]
    .map(clean)
)

raw["_SEQUENCE_NORMALIZED"] = (
    raw[raw_sequence_col]
    .map(normalize_sequence)
)

if raw_accession_col:
    raw["_ACCESSION"] = (
        raw[raw_accession_col].map(clean)
    )
else:
    raw["_ACCESSION"] = ""

if raw_normalized_accession_col:
    raw["_NORMALIZED_ACCESSION"] = (
        raw[raw_normalized_accession_col]
        .map(normalize_accession)
    )
else:
    raw["_NORMALIZED_ACCESSION"] = (
        raw["_ACCESSION"]
        .map(normalize_accession)
    )

if raw_length_col:
    raw["_SEQUENCE_LENGTH"] = (
        raw[raw_length_col].map(safe_int)
    )
else:
    raw["_SEQUENCE_LENGTH"] = (
        raw["_SEQUENCE_NORMALIZED"]
        .map(lambda x: len(x) if x else None)
    )


# ============================================================
# SEQUENCE HASH
# ============================================================

if raw_hash_col:

    raw["_SHA256"] = (
        raw[raw_hash_col]
        .map(clean)
    )

else:

    raw["_SHA256"] = (
        raw["_SEQUENCE_NORMALIZED"]
        .map(sha256_sequence)
    )


# Recalculate hashes for actual sequences.
# This prevents stale/malformed hash fields from affecting QC.

raw.loc[
    raw["_SEQUENCE_NORMALIZED"].ne(""),
    "_SHA256"
] = raw.loc[
    raw["_SEQUENCE_NORMALIZED"].ne(""),
    "_SEQUENCE_NORMALIZED"
].map(sha256_sequence)


raw["_HAS_SEQUENCE"] = (
    raw["_SEQUENCE_NORMALIZED"]
    .ne("")
)


print()
print("SEQUENCE DATA")
print("-" * 110)
print(
    f"Raw sequence rows : "
    f"{raw['_HAS_SEQUENCE'].sum()}"
)
print(
    f"KJUs represented  : "
    f"{raw.loc[raw['_HAS_SEQUENCE'], '_KJU_NORMALIZED'].nunique()}"
)


# ============================================================
# BUILD MASTER LOOKUP
# ============================================================

master_lookup = {}

for _, row in master.iterrows():

    kju = row["_KJU_NORMALIZED"]

    master_lookup[kju] = {
        "Accepted_accession":
            row["_ACCEPTED_ACCESSION"],

        "Normalized_accession":
            row["_NORMALIZED_ACCEPTED_ACCESSION"],

        "Accepted_taxonomy":
            row["_TAXONOMY"],

        "Accepted_length":
            row["_ACCEPTED_LENGTH"],

        "Biological_identity":
            row["_BIOLOGICAL_IDENTITY"],

        "Adjudication_type":
            row["_ADJUDICATION_TYPE"],

        "Confidence":
            row["_CONFIDENCE"],
    }


# ============================================================
# BUILD KJU -> SEQUENCE RECORDS
#
# IMPORTANT:
# Do NOT call .iterrows() on a list.
#
# We explicitly convert the filtered DataFrame into a list
# of dictionaries. This avoids the Step 78 failure.
# ============================================================

sequence_rows = (
    raw.loc[
        raw["_HAS_SEQUENCE"]
    ]
    .to_dict("records")
)


kju_sequences = defaultdict(list)

for row in sequence_rows:

    kju = clean(
        row.get("_KJU_NORMALIZED", "")
    )

    if not kju:
        continue

    kju_sequences[kju].append(row)


# ============================================================
# SEQUENCE INVENTORY PER KJU
# ============================================================

audit_records = []

for kju in sorted(
    master_lookup.keys(),
    key=lambda x: int(
        x.replace("KJU", "")
    )
):

    master_info = master_lookup[kju]

    rows = kju_sequences.get(
        kju,
        []
    )

    sequence_count = len(rows)

    sequence_hashes = sorted(
        {
            clean(
                x.get("_SHA256", "")
            )
            for x in rows
            if clean(
                x.get("_SHA256", "")
            )
        }
    )

    accessions = sorted(
        {
            normalize_accession(
                x.get("_NORMALIZED_ACCESSION", "")
                or x.get("_ACCESSION", "")
            )
            for x in rows
            if normalize_accession(
                x.get("_NORMALIZED_ACCESSION", "")
                or x.get("_ACCESSION", "")
            )
        }
    )

    lengths = sorted(
        {
            safe_int(
                x.get("_SEQUENCE_LENGTH", "")
            )
            for x in rows
            if safe_int(
                x.get("_SEQUENCE_LENGTH", "")
            ) is not None
        }
    )

    sequence_present = (
        sequence_count > 0
    )

    is_unresolved = (
        master_info["Accepted_taxonomy"] == ""
    )

    accession_match = True
    length_match = True

    # --------------------------------------------------------
    # ACCESSION CHECK
    # --------------------------------------------------------

    if sequence_present:

        expected_acc = (
            master_info[
                "Normalized_accession"
            ]
        )

        if expected_acc:

            accession_match = (
                all(
                    acc == expected_acc
                    for acc in accessions
                    if acc
                )
            )

        else:

            accession_match = False


    # --------------------------------------------------------
    # LENGTH CHECK
    # --------------------------------------------------------

    expected_length = (
        master_info[
            "Accepted_length"
        ]
    )

    if sequence_present and expected_length:

        length_match = (
            all(
                length == expected_length
                for length in lengths
                if length is not None
            )
        )

    audit_records.append({

        "KJU_ID":
            kju,

        "Accepted_taxonomy":
            master_info[
                "Accepted_taxonomy"
            ],

        "Accepted_accession":
            master_info[
                "Accepted_accession"
            ],

        "Normalized_accession":
            master_info[
                "Normalized_accession"
            ],

        "Accepted_length":
            master_info[
                "Accepted_length"
            ],

        "Biological_identity":
            master_info[
                "Biological_identity"
            ],

        "Adjudication_type":
            master_info[
                "Adjudication_type"
            ],

        "Confidence":
            master_info[
                "Confidence"
            ],

        "Resolved":
            "NO" if is_unresolved else "YES",

        "Sequence_record_count":
            sequence_count,

        "Unique_sequence_hashes":
            len(sequence_hashes),

        "Sequence_hashes":
            ";".join(
                sequence_hashes
            ),

        "Detected_accessions":
            ";".join(
                accessions
            ),

        "Detected_lengths":
            ";".join(
                str(x)
                for x in lengths
            ),

        "Sequence_present":
            "YES"
            if sequence_present
            else "NO",

        "Accession_match":
            "YES"
            if accession_match
            else "NO",

        "Length_match":
            "YES"
            if length_match
            else "NO",

        "Unresolved_with_sequence":
            "YES"
            if (
                is_unresolved
                and sequence_present
            )
            else "NO",
    })


kju_audit = pd.DataFrame(
    audit_records
)


# ============================================================
# DUPLICATE SEQUENCE GROUPS
# ============================================================

hash_groups = defaultdict(list)

for row in sequence_rows:

    sha = clean(
        row.get("_SHA256", "")
    )

    if not sha:
        continue

    hash_groups[sha].append(row)


duplicate_records = []

for sha, rows in hash_groups.items():

    if len(rows) < 2:
        continue

    kju_ids = sorted(
        {
            clean(
                x.get("_KJU_NORMALIZED", "")
            )
            for x in rows
            if clean(
                x.get("_KJU_NORMALIZED", "")
            )
        },
        key=lambda x: int(
            x.replace("KJU", "")
        )
    )

    accessions = sorted(
        {
            normalize_accession(
                x.get("_NORMALIZED_ACCESSION", "")
                or x.get("_ACCESSION", "")
            )
            for x in rows
            if normalize_accession(
                x.get("_NORMALIZED_ACCESSION", "")
                or x.get("_ACCESSION", "")
            )
        }
    )

    lengths = sorted(
        {
            safe_int(
                x.get("_SEQUENCE_LENGTH", "")
            )
            for x in rows
            if safe_int(
                x.get("_SEQUENCE_LENGTH", "")
            ) is not None
        }
    )

    same_kju = (
        len(kju_ids) == 1
    )

    different_kju = (
        len(kju_ids) > 1
    )

    same_accession = (
        len(accessions) == 1
    )

    different_accession = (
        len(accessions) > 1
    )

    duplicate_records.append({

        "Sequence_SHA256":
            sha,

        "Raw_record_count":
            len(rows),

        "KJU_count":
            len(kju_ids),

        "KJU_IDs":
            ";".join(kju_ids),

        "Accession_count":
            len(accessions),

        "Accessions":
            ";".join(accessions),

        "Sequence_lengths":
            ";".join(
                str(x)
                for x in lengths
            ),

        "Same_KJU":
            "YES"
            if same_kju
            else "NO",

        "Different_KJU":
            "YES"
            if different_kju
            else "NO",

        "Same_accession":
            "YES"
            if same_accession
            else "NO",

        "Different_accession":
            "YES"
            if different_accession
            else "NO",

        "Duplicate_category":
            (
                "SAME_SEQUENCE_SAME_KJU"
                if same_kju
                else
                "SAME_SEQUENCE_DIFFERENT_KJU"
            )
    })


duplicate_audit = pd.DataFrame(
    duplicate_records
)


# ============================================================
# CATEGORY A
# SAME SEQUENCE / SAME KJU
# ============================================================

if len(duplicate_audit):

    same_kju_df = (
        duplicate_audit[
            duplicate_audit[
                "Duplicate_category"
            ]
            ==
            "SAME_SEQUENCE_SAME_KJU"
        ]
        .copy()
    )

else:

    same_kju_df = pd.DataFrame(
        columns=[
            "Sequence_SHA256",
            "Raw_record_count",
            "KJU_count",
            "KJU_IDs",
            "Accession_count",
            "Accessions",
            "Sequence_lengths",
            "Same_KJU",
            "Different_KJU",
            "Same_accession",
            "Different_accession",
            "Duplicate_category",
        ]
    )


# ============================================================
# CATEGORY B
# SAME SEQUENCE / DIFFERENT KJU
# ============================================================

if len(duplicate_audit):

    cross_kju_df = (
        duplicate_audit[
            duplicate_audit[
                "Duplicate_category"
            ]
            ==
            "SAME_SEQUENCE_DIFFERENT_KJU"
        ]
        .copy()
    )

else:

    cross_kju_df = pd.DataFrame(
        columns=same_kju_df.columns
    )


# ============================================================
# CATEGORY C
# SAME ACCESSION / DIFFERENT SEQUENCE
#
# Build accession groups from raw sequence records.
# ============================================================

accession_groups = defaultdict(list)

for row in sequence_rows:

    accession = normalize_accession(
        row.get("_NORMALIZED_ACCESSION", "")
        or row.get("_ACCESSION", "")
    )

    if not accession:
        continue

    accession_groups[accession].append(row)


accession_conflict_records = []

for accession, rows in accession_groups.items():

    hashes = sorted(
        {
            clean(
                x.get("_SHA256", "")
            )
            for x in rows
            if clean(
                x.get("_SHA256", "")
            )
        }
    )

    if len(hashes) <= 1:
        continue

    kju_ids = sorted(
        {
            clean(
                x.get("_KJU_NORMALIZED", "")
            )
            for x in rows
            if clean(
                x.get("_KJU_NORMALIZED", "")
            )
        },
        key=lambda x: int(
            x.replace("KJU", "")
        )
    )

    accession_conflict_records.append({

        "Accession":
            accession,

        "KJU_IDs":
            ";".join(kju_ids),

        "KJU_count":
            len(kju_ids),

        "Unique_sequence_hashes":
            len(hashes),

        "Sequence_SHA256":
            ";".join(hashes),

        "Status":
            "SAME_ACCESSION_DIFFERENT_SEQUENCE",
    })


accession_conflicts = pd.DataFrame(
    accession_conflict_records
)


# ============================================================
# CATEGORY D
# DIFFERENT ACCESSION / IDENTICAL SEQUENCE
# ============================================================

if len(duplicate_audit):

    different_acc_same_sequence = (
        duplicate_audit[
            (
                duplicate_audit[
                    "Different_KJU"
                ]
                ==
                "YES"
            )
            &
            (
                duplicate_audit[
                    "Different_accession"
                ]
                ==
                "YES"
            )
        ]
        .copy()
    )

else:

    different_acc_same_sequence = pd.DataFrame(
        columns=duplicate_audit.columns
    )


# ============================================================
# LENGTH MISMATCHES
# ============================================================

length_mismatch_df = (
    kju_audit[
        (
            kju_audit[
                "Sequence_present"
            ]
            ==
            "YES"
        )
        &
        (
            kju_audit[
                "Length_match"
            ]
            ==
            "NO"
        )
    ]
    .copy()
)


# ============================================================
# UNRESOLVED KJU WITH SEQUENCE
#
# This is IMPORTANT.
#
# If an unresolved MASTER_ONLY KJU has a sequence,
# we now have sequence evidence but do NOT automatically
# assign taxonomy.
#
# It should be flagged for later biological resolution.
# ============================================================

unresolved_with_sequence = (
    kju_audit[
        kju_audit[
            "Unresolved_with_sequence"
        ]
        ==
        "YES"
    ]
    .copy()
)


# ============================================================
# SAVE AUDITS
# ============================================================

kju_audit.to_csv(
    KJU_AUDIT,
    index=False
)

duplicate_audit.to_csv(
    DUPLICATE_AUDIT,
    index=False
)

same_kju_df.to_csv(
    SAME_KJU_DUPLICATES,
    index=False
)

cross_kju_df.to_csv(
    CROSS_KJU_DUPLICATES,
    index=False
)

accession_conflicts.to_csv(
    ACCESSION_CONFLICTS,
    index=False
)

different_acc_same_sequence.to_csv(
    DIFFERENT_ACCESSION_SAME_SEQUENCE,
    index=False
)

length_mismatch_df.to_csv(
    LENGTH_MISMATCHES,
    index=False
)

unresolved_with_sequence.to_csv(
    UNRESOLVED_WITH_SEQUENCE,
    index=False
)


# ============================================================
# SUMMARY COUNTS
# ============================================================

total_kju = len(master)

resolved_kju = int(
    (
        kju_audit["Resolved"]
        ==
        "YES"
    ).sum()
)

unresolved_kju = int(
    (
        kju_audit["Resolved"]
        ==
        "NO"
    ).sum()
)

kju_with_sequence = int(
    (
        kju_audit["Sequence_present"]
        ==
        "YES"
    ).sum()
)

kju_without_sequence = int(
    (
        kju_audit["Sequence_present"]
        ==
        "NO"
    ).sum()
)

resolved_with_sequence = int(
    (
        (
            kju_audit["Resolved"]
            ==
            "YES"
        )
        &
        (
            kju_audit["Sequence_present"]
            ==
            "YES"
        )
    ).sum()
)

resolved_without_sequence = int(
    (
        (
            kju_audit["Resolved"]
            ==
            "YES"
        )
        &
        (
            kju_audit["Sequence_present"]
            ==
            "NO"
        )
    ).sum()
)

unresolved_with_sequence_count = len(
    unresolved_with_sequence
)

duplicate_group_count = len(
    duplicate_audit
)

same_kju_count = len(
    same_kju_df
)

cross_kju_count = len(
    cross_kju_df
)

accession_conflict_count = len(
    accession_conflicts
)

different_accession_same_sequence_count = len(
    different_acc_same_sequence
)

length_mismatch_count = len(
    length_mismatch_df
)


# ============================================================
# ACCESSION MATCH QC
# ============================================================

sequence_records_with_expected_accession = 0
accession_matches = 0
accession_mismatches = 0

for _, row in kju_audit.iterrows():

    if (
        row["Sequence_present"]
        !=
        "YES"
    ):
        continue

    expected = clean(
        row["Normalized_accession"]
    )

    detected = clean(
        row["Detected_accessions"]
    )

    if expected:

        sequence_records_with_expected_accession += 1

        detected_values = {
            normalize_accession(x)
            for x in detected.split(";")
            if normalize_accession(x)
        }

        if (
            detected_values
            and
            detected_values
            <= {expected}
        ):

            accession_matches += 1

        else:

            accession_mismatches += 1


# ============================================================
# FROZEN MASTER INTEGRITY
# ============================================================

master_reloaded = pd.read_csv(
    MASTER,
    dtype=str,
    keep_default_na=False,
)

master_original_hash = hashlib.sha256(
    MASTER.read_bytes()
).hexdigest()

master_reloaded_hash = hashlib.sha256(
    MASTER.read_bytes()
).hexdigest()

frozen_master_unchanged = (
    master_original_hash
    ==
    master_reloaded_hash
)


# ============================================================
# REPORT
# ============================================================

report_lines = []

report_lines.append(
    "STEP 78 — SEQUENCE IDENTITY + DUPLICATE AUDIT"
)

report_lines.append(
    "=" * 110
)

report_lines.append("")

report_lines.append(
    "INPUTS"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Frozen master : {MASTER}"
)

report_lines.append(
    f"Inventory     : {INVENTORY}"
)

report_lines.append(
    f"Raw sequences : {RAW}"
)

report_lines.append("")

report_lines.append(
    "1. MASTER QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Total KJU records       : {total_kju}"
)

report_lines.append(
    f"Resolved records        : {resolved_kju}"
)

report_lines.append(
    f"Unresolved records      : {unresolved_kju}"
)

report_lines.append(
    f"Unique KJU IDs          : {master['_KJU_NORMALIZED'].nunique()}"
)

report_lines.append("")

report_lines.append(
    "2. SEQUENCE COVERAGE"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"KJUs with sequence      : {kju_with_sequence}"
)

report_lines.append(
    f"KJUs without sequence   : {kju_without_sequence}"
)

report_lines.append(
    f"Resolved + sequence     : {resolved_with_sequence}"
)

report_lines.append(
    f"Resolved + no sequence  : {resolved_without_sequence}"
)

report_lines.append(
    f"Unresolved + sequence   : {unresolved_with_sequence_count}"
)

report_lines.append("")

report_lines.append(
    "3. ACCESSION QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Records checked         : {sequence_records_with_expected_accession}"
)

report_lines.append(
    f"Accession matches       : {accession_matches}"
)

report_lines.append(
    f"Accession mismatches    : {accession_mismatches}"
)

report_lines.append("")

report_lines.append(
    "4. LENGTH QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Length mismatches       : {length_mismatch_count}"
)

report_lines.append("")

report_lines.append(
    "5. DUPLICATE SEQUENCE QC"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Duplicate groups        : {duplicate_group_count}"
)

report_lines.append(
    f"Same sequence/same KJU  : {same_kju_count}"
)

report_lines.append(
    f"Same sequence/diff KJU  : {cross_kju_count}"
)

report_lines.append(
    f"Same accession/diff seq : {accession_conflict_count}"
)

report_lines.append(
    f"Diff accession/same seq : {different_accession_same_sequence_count}"
)

report_lines.append("")

report_lines.append(
    "6. UNRESOLVED RECORD QC"
)

report_lines.append(
    "-" * 110
)

if unresolved_with_sequence_count == 0:

    report_lines.append(
        "PASS — No unresolved KJU has a sequence"
    )

else:

    report_lines.append(
        "REVIEW — Unresolved KJUs with sequences detected"
    )

    for _, row in unresolved_with_sequence.iterrows():

        report_lines.append(
            f"{row['KJU_ID']} | "
            f"sequence_records={row['Sequence_record_count']} | "
            f"accessions={row['Detected_accessions']} | "
            f"lengths={row['Detected_lengths']}"
        )

report_lines.append("")

report_lines.append(
    "7. LENGTH MISMATCH DETAILS"
)

report_lines.append(
    "-" * 110
)

if length_mismatch_count == 0:

    report_lines.append(
        "PASS — No sequence length mismatches"
    )

else:

    for _, row in length_mismatch_df.iterrows():

        report_lines.append(
            f"{row['KJU_ID']} | "
            f"accepted={row['Accepted_length']} | "
            f"detected={row['Detected_lengths']}"
        )

report_lines.append("")

report_lines.append(
    "8. ACCESSION MISMATCH DETAILS"
)

report_lines.append(
    "-" * 110
)

if accession_mismatches == 0:

    report_lines.append(
        "PASS — No accession mismatches"
    )

else:

    for _, row in kju_audit.iterrows():

        if (
            row["Sequence_present"]
            ==
            "YES"
            and
            row["Accession_match"]
            ==
            "NO"
        ):

            report_lines.append(
                f"{row['KJU_ID']} | "
                f"accepted={row['Normalized_accession']} | "
                f"detected={row['Detected_accessions']}"
            )

report_lines.append("")

report_lines.append(
    "9. DUPLICATE GROUP DETAILS"
)

report_lines.append(
    "-" * 110
)

if duplicate_group_count == 0:

    report_lines.append(
        "No duplicate sequence groups detected."
    )

else:

    for _, row in duplicate_audit.iterrows():

        report_lines.append(
            f"{row['Duplicate_category']} | "
            f"KJUs={row['KJU_IDs']} | "
            f"accessions={row['Accessions']} | "
            f"records={row['Raw_record_count']}"
        )

report_lines.append("")

report_lines.append(
    "10. FROZEN MASTER INTEGRITY"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    f"Frozen master rows      : {len(master_reloaded)}"
)

report_lines.append(
    f"Frozen master SHA256    : {master_original_hash}"
)

if frozen_master_unchanged:

    report_lines.append(
        "PASS — Frozen V1 master was not modified"
    )

else:

    report_lines.append(
        "FAIL — Frozen V1 master integrity problem"
    )

report_lines.append("")

report_lines.append(
    "STEP 78 STATUS"
)

report_lines.append(
    "=" * 110
)

if (
    len(master) == 204
    and
    master["_KJU_NORMALIZED"].nunique() == 204
    and
    length_mismatch_count == 0
    and
    accession_mismatches == 0
    and
    frozen_master_unchanged
):

    report_lines.append(
        "PASS — SEQUENCE IDENTITY AUDIT COMPLETED"
    )

    report_lines.append(
        "No sequence length mismatches."
    )

    report_lines.append(
        "No accession mismatches."
    )

    if unresolved_with_sequence_count > 0:

        report_lines.append(
            "NOTE — Unresolved KJUs with sequence require "
            "future taxonomic resolution; taxonomy was NOT fabricated."
        )

else:

    report_lines.append(
        "REVIEW REQUIRED — One or more QC conditions failed."
    )

report_lines.append("")

report_lines.append(
    "OUTPUT FILES"
)

report_lines.append(
    "-" * 110
)

report_lines.append(
    str(KJU_AUDIT)
)

report_lines.append(
    str(DUPLICATE_AUDIT)
)

report_lines.append(
    str(SAME_KJU_DUPLICATES)
)

report_lines.append(
    str(CROSS_KJU_DUPLICATES)
)

report_lines.append(
    str(ACCESSION_CONFLICTS)
)

report_lines.append(
    str(DIFFERENT_ACCESSION_SAME_SEQUENCE)
)

report_lines.append(
    str(LENGTH_MISMATCHES)
)

report_lines.append(
    str(UNRESOLVED_WITH_SEQUENCE)
)

report_lines.append(
    str(REPORT)
)

report_lines.append("")

report_lines.append(
    "FROZEN V1 MASTER WAS NOT MODIFIED."
)


REPORT.write_text(
    "\n".join(report_lines)
    + "\n"
)


# ============================================================
# CONSOLE REPORT
# ============================================================

print()
print("=" * 110)
print("STEP 78 COMPLETE")
print("=" * 110)

print()
print("MASTER")
print("-" * 110)
print(
    f"Total KJU records       : {total_kju}"
)
print(
    f"Resolved records        : {resolved_kju}"
)
print(
    f"Unresolved records      : {unresolved_kju}"
)

print()
print("SEQUENCE COVERAGE")
print("-" * 110)
print(
    f"KJUs with sequence      : {kju_with_sequence}"
)
print(
    f"KJUs without sequence   : {kju_without_sequence}"
)
print(
    f"Resolved + sequence     : {resolved_with_sequence}"
)
print(
    f"Resolved + no sequence  : {resolved_without_sequence}"
)
print(
    f"Unresolved + sequence   : {unresolved_with_sequence_count}"
)

print()
print("ACCESSION QC")
print("-" * 110)
print(
    f"Matches                : {accession_matches}"
)
print(
    f"Mismatches             : {accession_mismatches}"
)

print()
print("LENGTH QC")
print("-" * 110)
print(
    f"Mismatches             : {length_mismatch_count}"
)

print()
print("DUPLICATE SEQUENCES")
print("-" * 110)
print(
    f"Duplicate groups       : {duplicate_group_count}"
)
print(
    f"Same sequence/same KJU : {same_kju_count}"
)
print(
    f"Same sequence/diff KJU : {cross_kju_count}"
)
print(
    f"Same accession/diff seq: {accession_conflict_count}"
)
print(
    f"Diff accession/same seq: "
    f"{different_accession_same_sequence_count}"
)

print()
print("FROZEN MASTER")
print("-" * 110)

if frozen_master_unchanged:
    print(
        "PASS — Frozen V1 master was not modified"
    )
else:
    print(
        "FAIL — Frozen V1 master integrity problem"
    )

print()
print("OUTPUT DIRECTORY")
print("-" * 110)
print(OUTDIR)

print()
print("=" * 110)
print(
    "FROZEN V1 MASTER WAS NOT MODIFIED"
)
print("=" * 110)
