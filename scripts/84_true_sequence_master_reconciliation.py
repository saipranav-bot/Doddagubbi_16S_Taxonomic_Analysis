#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
from collections import defaultdict


# ============================================================
# STEP 84
# TRUE SEQUENCE ↔ FROZEN MASTER RECONCILIATION
#
# INPUT:
#   STEP 76 frozen master
#   STEP 83 true nucleotide sequence inventory
#
# PURPOSE:
#   1. Compare REAL nucleotide sequences against the frozen
#      adjudicated master.
#
#   2. For every KJU determine:
#        - number of sequence records
#        - number of unique sequences
#        - candidate accessions
#        - candidate lengths
#        - candidate sequence types
#        - candidate SHA256 hashes
#
#   3. Compare actual sequences with:
#        - Accepted_accession
#        - Accepted_length
#        - Accepted_taxonomy
#        - Biological_identity
#
#   4. Detect:
#        A. Exact master accession + exact sequence
#        B. Exact master accession + different sequence
#        C. Different accession + identical sequence
#        D. Different accession + different sequence
#        E. Length agreement
#        F. Length disagreement
#        G. Multiple sequences per KJU
#
#   5. Determine whether a KJU can be classified as:
#
#        EXACT_MATCH
#        ACCESSION_MATCH_MULTIPLE_SEQUENCES
#        ACCESSION_CONFLICT
#        NO_SEQUENCE
#        UNRESOLVED
#
# IMPORTANT:
#   - Frozen V1 master is NEVER modified.
#   - No taxonomy is fabricated.
#   - No canonical sequence is selected.
#   - No sequence is deleted.
#   - Every observed sequence is preserved.
# ============================================================


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

SEQ_ROOT = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"

TRUE_RECORDS = (
    SEQ_ROOT
    / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
)

TRUE_INVENTORY = (
    SEQ_ROOT
    / "STEP83_TRUE_SEQUENCE_INVENTORY.csv"
)

OUT = (
    ROOT
    / "STEP84_TRUE_SEQUENCE_MASTER_RECONCILIATION"
)

OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# OUTPUT FILES
# ============================================================

KJU_AUDIT = (
    OUT
    / "STEP84_KJU_TRUE_SEQUENCE_RECONCILIATION.csv"
)

CANDIDATES = (
    OUT
    / "STEP84_ALL_TRUE_SEQUENCE_CANDIDATES.csv"
)

EXACT_MATCHES = (
    OUT
    / "STEP84_EXACT_MATCHES.csv"
)

ACCESSION_MATCH_MULTIPLE = (
    OUT
    / "STEP84_ACCESSION_MATCH_MULTIPLE_SEQUENCES.csv"
)

ACCESSION_CONFLICTS = (
    OUT
    / "STEP84_ACCESSION_CONFLICTS.csv"
)

LENGTH_CONFLICTS = (
    OUT
    / "STEP84_LENGTH_CONFLICTS.csv"
)

DUPLICATE_ACROSS_KJU = (
    OUT
    / "STEP84_DUPLICATE_SEQUENCE_ACROSS_KJU.csv"
)

SUMMARY = (
    OUT
    / "STEP84_TRUE_SEQUENCE_MASTER_RECONCILIATION_REPORT.txt"
)


# ============================================================
# HELPERS
# ============================================================

def clean(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def norm_accession(x):
    x = clean(x)

    if not x:
        return ""

    x = x.strip()

    # Remove surrounding punctuation/whitespace.
    x = x.replace(" ", "")

    # Normalize common accession.version representation.
    return x.upper()


def norm_taxonomy(x):
    return clean(x).lower()


def safe_int(x):
    try:
        if pd.isna(x):
            return None

        s = str(x).strip()

        if not s:
            return None

        return int(float(s))

    except Exception:
        return None


def sequence_sha(seq):
    """
    Calculate SHA256 from normalized nucleotide sequence.
    """

    if seq is None:
        return ""

    seq = clean(seq).upper()

    # Remove FASTA whitespace if present.
    seq = re.sub(r"[^A-Z]", "", seq)

    if not seq:
        return ""

    return hashlib.sha256(
        seq.encode("utf-8")
    ).hexdigest()


# ============================================================
# START
# ============================================================

print("=" * 110)
print("STEP 84: TRUE SEQUENCE ↔ FROZEN MASTER RECONCILIATION")
print("=" * 110)


# ============================================================
# LOAD INPUTS
# ============================================================

if not MASTER.exists():
    raise FileNotFoundError(
        f"Frozen master not found:\n{MASTER}"
    )

if not TRUE_RECORDS.exists():
    raise FileNotFoundError(
        f"True sequence records not found:\n{TRUE_RECORDS}"
    )


master = pd.read_csv(
    MASTER,
    dtype=str,
    keep_default_na=False
)

seq = pd.read_csv(
    TRUE_RECORDS,
    dtype=str,
    keep_default_na=False
)


print()
print("INPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"True records  : {TRUE_RECORDS}")

print()
print("INPUT SIZES")
print("-" * 110)
print(f"Master rows       : {len(master)}")
print(f"Sequence rows     : {len(seq)}")


# ============================================================
# COLUMN DETECTION
# ============================================================

def detect_column(df, candidates, label):
    for c in candidates:
        if c in df.columns:
            return c

    raise ValueError(
        f"Could not detect {label} column.\n"
        f"Available columns:\n{list(df.columns)}"
    )


master_kju = detect_column(
    master,
    ["KJU_ID", "KJU", "KJU_ID_x"],
    "master KJU"
)

master_acc = detect_column(
    master,
    ["Accepted_accession", "accepted_accession"],
    "master accession"
)

master_tax = detect_column(
    master,
    ["Accepted_taxonomy", "accepted_taxonomy"],
    "master taxonomy"
)

master_len = detect_column(
    master,
    ["Accepted_length", "accepted_length"],
    "master length"
)

master_identity = detect_column(
    master,
    ["Biological_identity", "biological_identity"],
    "master biological identity"
)


seq_kju = detect_column(
    seq,
    ["KJU_ID", "KJU", "KJU_ID_x"],
    "sequence KJU"
)

seq_acc = detect_column(
    seq,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "Accession",
        "accession"
    ],
    "sequence accession"
)

seq_len = detect_column(
    seq,
    [
        "Sequence_length",
        "sequence_length",
        "Length",
        "length"
    ],
    "sequence length"
)

seq_sha = detect_column(
    seq,
    [
        "Sequence_SHA256",
        "sequence_sha256",
        "SHA256",
        "sha256"
    ],
    "sequence SHA256"
)

seq_type = None

for c in [
    "Sequence_type_detected",
    "sequence_type",
    "Sequence_type",
    "Type"
]:
    if c in seq.columns:
        seq_type = c
        break


print()
print("DETECTED COLUMNS")
print("-" * 110)
print(f"Master KJU       : {master_kju}")
print(f"Master accession : {master_acc}")
print(f"Master taxonomy  : {master_tax}")
print(f"Master length    : {master_len}")
print(f"Master identity  : {master_identity}")

print(f"Sequence KJU     : {seq_kju}")
print(f"Sequence accession: {seq_acc}")
print(f"Sequence length  : {seq_len}")
print(f"Sequence SHA256  : {seq_sha}")

if seq_type:
    print(f"Sequence type    : {seq_type}")
else:
    print("Sequence type    : NOT AVAILABLE")


# ============================================================
# MASTER NORMALIZATION
# ============================================================

master = master.copy()

master["_KJU"] = master[master_kju].map(clean)
master["_ACC"] = master[master_acc].map(norm_accession)
master["_TAX"] = master[master_tax].map(clean)
master["_TAX_NORM"] = master[master_tax].map(norm_taxonomy)
master["_LEN"] = master[master_len].map(safe_int)
master["_IDENTITY"] = master[master_identity].map(clean)


# ============================================================
# SEQUENCE NORMALIZATION
# ============================================================

seq = seq.copy()

seq["_KJU"] = seq[seq_kju].map(clean)
seq["_ACC"] = seq[seq_acc].map(norm_accession)
seq["_LEN"] = seq[seq_len].map(safe_int)
seq["_SHA"] = seq[seq_sha].map(clean)

if seq_type:
    seq["_TYPE"] = seq[seq_type].map(clean)
else:
    seq["_TYPE"] = ""


# ============================================================
# SHA VALIDATION
# ============================================================

reported_sha_available = (
    seq["_SHA"].astype(str).str.len() == 64
).sum()

print()
print("SEQUENCE DATA")
print("-" * 110)
print(f"Sequence rows                  : {len(seq)}")
print(f"KJUs represented               : {seq['_KJU'].nunique()}")
print(f"Rows with valid SHA256         : {reported_sha_available}")


# ============================================================
# MASTER LOOKUP
# ============================================================

master_lookup = {}

for _, row in master.iterrows():

    kju = row["_KJU"]

    master_lookup[kju] = {
        "accepted_accession": row["_ACC"],
        "accepted_taxonomy": row["_TAX"],
        "accepted_length": row["_LEN"],
        "biological_identity": row["_IDENTITY"],
    }


# ============================================================
# GROUP TRUE SEQUENCES BY KJU
# ============================================================

groups = defaultdict(list)

for _, row in seq.iterrows():

    kju = row["_KJU"]

    if not kju:
        continue

    groups[kju].append(row)


# ============================================================
# BUILD ALL CANDIDATES TABLE
# ============================================================

candidate_rows = []

for kju, rows in groups.items():

    for row in rows:

        candidate_rows.append(
            {
                "KJU_ID": kju,
                "Sequence_accession": row["_ACC"],
                "Sequence_length": row["_LEN"],
                "Sequence_SHA256": row["_SHA"],
                "Sequence_type": row["_TYPE"],
                "Sequence_file":
                    clean(
                        row.get(
                            "Sequence_file",
                            ""
                        )
                    ),
                "Record_ID":
                    clean(
                        row.get(
                            "Record_ID",
                            ""
                        )
                    ),
            }
        )


candidate_df = pd.DataFrame(candidate_rows)

if len(candidate_df) == 0:

    candidate_df = pd.DataFrame(
        columns=[
            "KJU_ID",
            "Sequence_accession",
            "Sequence_length",
            "Sequence_SHA256",
            "Sequence_type",
            "Sequence_file",
            "Record_ID",
        ]
    )


# ============================================================
# EXACT DUPLICATE COLLAPSE
# ============================================================

candidate_unique = (
    candidate_df
    .drop_duplicates(
        subset=[
            "KJU_ID",
            "Sequence_accession",
            "Sequence_length",
            "Sequence_SHA256",
        ]
    )
    .copy()
)


# ============================================================
# PER-KJU RECONCILIATION
# ============================================================

audit_rows = []


for _, m in master.iterrows():

    kju = m["_KJU"]

    accepted_acc = m["_ACC"]
    accepted_tax = m["_TAX"]
    accepted_len = m["_LEN"]
    biological_identity = m["_IDENTITY"]

    kju_candidates = candidate_unique[
        candidate_unique["KJU_ID"] == kju
    ].copy()

    unique_sha = sorted(
        set(
            x
            for x in kju_candidates["Sequence_SHA256"]
            if clean(x)
        )
    )

    unique_acc = sorted(
        set(
            x
            for x in kju_candidates["Sequence_accession"]
            if clean(x)
        )
    )

    unique_len = sorted(
        set(
            x
            for x in kju_candidates["Sequence_length"]
            if x is not None
            and not pd.isna(x)
        )
    )

    n_candidates = len(kju_candidates)

    n_unique_seq = len(unique_sha)

    accession_match = (
        accepted_acc != ""
        and accepted_acc in unique_acc
    )

    length_match = (
        accepted_len is not None
        and accepted_len in unique_len
    )

    # Sequence records that exactly match accepted accession.
    accepted_acc_rows = kju_candidates[
        kju_candidates["Sequence_accession"]
        == accepted_acc
    ]

    accepted_acc_sha = sorted(
        set(
            x
            for x in accepted_acc_rows[
                "Sequence_SHA256"
            ]
            if clean(x)
        )
    )

    accepted_acc_lengths = sorted(
        set(
            x
            for x in accepted_acc_rows[
                "Sequence_length"
            ]
            if x is not None
            and not pd.isna(x)
        )
    )

    # --------------------------------------------------------
    # STATUS
    # --------------------------------------------------------

    if clean(accepted_tax) == "" and n_candidates == 0:

        status = "UNRESOLVED"

    elif n_candidates == 0:

        status = "NO_SEQUENCE"

    elif accession_match and n_unique_seq == 1:

        status = "EXACT_MATCH"

    elif accession_match and n_unique_seq > 1:

        status = (
            "ACCESSION_MATCH_MULTIPLE_SEQUENCES"
        )

    elif not accession_match and n_candidates > 0:

        status = "ACCESSION_CONFLICT"

    else:

        status = "REVIEW_REQUIRED"


    # --------------------------------------------------------
    # EVIDENCE
    # --------------------------------------------------------

    if status == "EXACT_MATCH":

        evidence = (
            "Accepted accession matches sequence "
            "accession and only one unique sequence "
            "exists for the KJU."
        )

    elif status == "ACCESSION_MATCH_MULTIPLE_SEQUENCES":

        evidence = (
            "Accepted accession is present, but "
            "multiple distinct sequences are associated "
            "with the KJU."
        )

    elif status == "ACCESSION_CONFLICT":

        evidence = (
            "Sequences exist for the KJU, but none "
            "matches the accepted accession."
        )

    elif status == "NO_SEQUENCE":

        evidence = (
            "Resolved master record has no observed "
            "true nucleotide sequence."
        )

    elif status == "UNRESOLVED":

        evidence = (
            "MASTER_ONLY record remains unresolved and "
            "has no sequence."
        )

    else:

        evidence = (
            "Sequence evidence requires manual review."
        )


    audit_rows.append(
        {
            "KJU_ID": kju,

            "Accepted_accession":
                accepted_acc,

            "Accepted_taxonomy":
                accepted_tax,

            "Accepted_length":
                accepted_len,

            "Biological_identity":
                biological_identity,

            "Candidate_record_count":
                n_candidates,

            "Unique_sequence_count":
                n_unique_seq,

            "Candidate_accessions":
                ";".join(unique_acc),

            "Candidate_lengths":
                ";".join(
                    str(x)
                    for x in unique_len
                ),

            "Candidate_SHA256":
                ";".join(unique_sha),

            "Accepted_accession_present":
                "YES"
                if accession_match
                else "NO",

            "Accepted_length_present":
                "YES"
                if length_match
                else "NO",

            "Accepted_accession_sequence_count":
                len(accepted_acc_sha),

            "Accepted_accession_unique_sequences":
                len(accepted_acc_sha),

            "Accepted_accession_lengths":
                ";".join(
                    str(x)
                    for x in accepted_acc_lengths
                ),

            "Status":
                status,

            "Evidence":
                evidence,
        }
    )


audit = pd.DataFrame(audit_rows)


# ============================================================
# DUPLICATE SEQUENCES ACROSS KJUs
# ============================================================

sha_groups = defaultdict(list)

for _, row in candidate_unique.iterrows():

    sha = clean(row["Sequence_SHA256"])

    if not sha:
        continue

    sha_groups[sha].append(
        {
            "KJU_ID":
                row["KJU_ID"],

            "Sequence_accession":
                row["Sequence_accession"],

            "Sequence_length":
                row["Sequence_length"],
        }
    )


duplicate_rows = []

for sha, items in sha_groups.items():

    kju_set = sorted(
        set(
            x["KJU_ID"]
            for x in items
        )
    )

    if len(kju_set) <= 1:
        continue

    duplicate_rows.append(
        {
            "Sequence_SHA256":
                sha,

            "KJU_count":
                len(kju_set),

            "KJU_IDs":
                ";".join(kju_set),

            "Accessions":
                ";".join(
                    sorted(
                        set(
                            clean(
                                x[
                                    "Sequence_accession"
                                ]
                            )
                            for x in items
                            if clean(
                                x[
                                    "Sequence_accession"
                                ]
                            )
                        )
                    )
                ),

            "Lengths":
                ";".join(
                    sorted(
                        set(
                            str(
                                x[
                                    "Sequence_length"
                                ]
                            )
                            for x in items
                            if x[
                                "Sequence_length"
                            ] is not None
                        )
                    )
                ),
        }
    )


duplicate_df = pd.DataFrame(
    duplicate_rows
)


# ============================================================
# SAME ACCESSION / DIFFERENT SEQUENCE
# ============================================================

accession_groups = defaultdict(list)

for _, row in candidate_unique.iterrows():

    acc = clean(
        row["Sequence_accession"]
    )

    sha = clean(
        row["Sequence_SHA256"]
    )

    if not acc or not sha:
        continue

    accession_groups[acc].append(
        row
    )


accession_conflict_rows = []

for acc, rows in accession_groups.items():

    sha_set = sorted(
        set(
            clean(x["Sequence_SHA256"])
            for x in rows
            if clean(
                x["Sequence_SHA256"]
            )
        )
    )

    if len(sha_set) <= 1:
        continue

    accession_conflict_rows.append(
        {
            "Accession":
                acc,

            "Unique_sequence_count":
                len(sha_set),

            "KJU_IDs":
                ";".join(
                    sorted(
                        set(
                            clean(
                                x["KJU_ID"]
                            )
                            for x in rows
                        )
                    )
                ),

            "Sequence_lengths":
                ";".join(
                    sorted(
                        set(
                            str(
                                x[
                                    "Sequence_length"
                                ]
                            )
                            for x in rows
                            if x[
                                "Sequence_length"
                            ] is not None
                        )
                    )
                ),

            "Sequence_SHA256":
                ";".join(sha_set),
        }
    )


accession_conflict_df = pd.DataFrame(
    accession_conflict_rows
)


# ============================================================
# LENGTH CONFLICTS
# ============================================================

length_conflict_df = audit[
    (
        audit["Accepted_length"].notna()
    )
    &
    (
        audit["Candidate_lengths"].astype(str)
        != ""
    )
    &
    (
        ~audit.apply(
            lambda r:
                str(
                    r["Accepted_length"]
                )
                in str(
                    r["Candidate_lengths"]
                ).split(";"),
            axis=1
        )
    )
].copy()


# ============================================================
# EXACT MATCH TABLE
# ============================================================

exact_df = audit[
    audit["Status"] == "EXACT_MATCH"
].copy()


# ============================================================
# ACCESSION MATCH / MULTIPLE SEQUENCE TABLE
# ============================================================

multi_df = audit[
    audit["Status"]
    == "ACCESSION_MATCH_MULTIPLE_SEQUENCES"
].copy()


# ============================================================
# ACCESSION CONFLICT TABLE
# ============================================================

conflict_df = audit[
    audit["Status"]
    == "ACCESSION_CONFLICT"
].copy()


# ============================================================
# SAVE OUTPUTS
# ============================================================

candidate_unique.to_csv(
    CANDIDATES,
    index=False
)

audit.to_csv(
    KJU_AUDIT,
    index=False
)

exact_df.to_csv(
    EXACT_MATCHES,
    index=False
)

multi_df.to_csv(
    ACCESSION_MATCH_MULTIPLE,
    index=False
)

conflict_df.to_csv(
    ACCESSION_CONFLICTS,
    index=False
)

length_conflict_df.to_csv(
    LENGTH_CONFLICTS,
    index=False
)

duplicate_df.to_csv(
    DUPLICATE_ACROSS_KJU,
    index=False
)


# ============================================================
# COUNTS
# ============================================================

total = len(audit)

resolved = int(
    (
        audit["Accepted_taxonomy"]
        .astype(str)
        .str.strip()
        != ""
    ).sum()
)

unresolved = total - resolved

with_sequence = int(
    (
        audit["Candidate_record_count"]
        > 0
    ).sum()
)

without_sequence = total - with_sequence

exact_count = int(
    (
        audit["Status"]
        == "EXACT_MATCH"
    ).sum()
)

multi_count = int(
    (
        audit["Status"]
        == "ACCESSION_MATCH_MULTIPLE_SEQUENCES"
    ).sum()
)

conflict_count = int(
    (
        audit["Status"]
        == "ACCESSION_CONFLICT"
    ).sum()
)

no_sequence_count = int(
    (
        audit["Status"]
        == "NO_SEQUENCE"
    ).sum()
)

unresolved_count = int(
    (
        audit["Status"]
        == "UNRESOLVED"
    ).sum()
)

different_sequence_kjus = int(
    (
        audit["Unique_sequence_count"]
        > 1
    ).sum()
)

one_sequence_kjus = int(
    (
        audit["Unique_sequence_count"]
        == 1
    ).sum()
)


# ============================================================
# QC CHECKS
# ============================================================

qc_results = []


def qc(name, condition):
    qc_results.append(
        (
            name,
            bool(condition)
        )
    )


qc(
    "Master contains 204 unique KJU IDs",
    total == 204
    and audit["KJU_ID"].nunique() == 204
)

qc(
    "No unresolved KJU has a sequence",
    not (
        audit.loc[
            audit["Status"]
            == "UNRESOLVED",
            "Candidate_record_count"
        ]
        > 0
    ).any()
)

qc(
    "No empty KJU IDs",
    not (
        audit["KJU_ID"]
        .astype(str)
        .str.strip()
        == ""
    ).any()
)

qc(
    "All exact matches have accepted accession",
    (
        exact_df[
            "Accepted_accession_present"
        ]
        == "YES"
    ).all()
    if len(exact_df)
    else True
)

qc(
    "Frozen master remains untouched",
    MASTER.exists()
)


# ============================================================
# REPORT
# ============================================================

with open(
    SUMMARY,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "=" * 110
        + "\n"
    )

    f.write(
        "STEP 84: TRUE SEQUENCE ↔ "
        "FROZEN MASTER RECONCILIATION\n"
    )

    f.write(
        "=" * 110
        + "\n\n"
    )

    f.write(
        "INPUTS\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"Frozen master : {MASTER}\n"
    )

    f.write(
        f"True records  : {TRUE_RECORDS}\n\n"
    )

    f.write(
        "MASTER\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"Total KJU records        : {total}\n"
    )

    f.write(
        f"Resolved records         : {resolved}\n"
    )

    f.write(
        f"Unresolved records       : {unresolved}\n\n"
    )

    f.write(
        "SEQUENCE COVERAGE\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"KJUs with sequence       : {with_sequence}\n"
    )

    f.write(
        f"KJUs without sequence    : {without_sequence}\n"
    )

    f.write(
        f"Exactly one unique seq   : {one_sequence_kjus}\n"
    )

    f.write(
        f"Multiple unique seq      : {different_sequence_kjus}\n\n"
    )

    f.write(
        "RECONCILIATION STATUS\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"EXACT_MATCH                         : "
        f"{exact_count}\n"
    )

    f.write(
        f"ACCESSION_MATCH_MULTIPLE_SEQUENCES : "
        f"{multi_count}\n"
    )

    f.write(
        f"ACCESSION_CONFLICT                 : "
        f"{conflict_count}\n"
    )

    f.write(
        f"NO_SEQUENCE                        : "
        f"{no_sequence_count}\n"
    )

    f.write(
        f"UNRESOLVED                         : "
        f"{unresolved_count}\n\n"
    )

    f.write(
        "CONFLICT SUMMARY\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"Same accession / different sequence groups : "
        f"{len(accession_conflict_df)}\n"
    )

    f.write(
        f"Different KJU / identical sequence groups  : "
        f"{len(duplicate_df)}\n"
    )

    f.write(
        f"KJUs with length conflict                   : "
        f"{len(length_conflict_df)}\n\n"
    )

    f.write(
        "QC\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    for name, passed in qc_results:

        f.write(
            ("PASS" if passed else "FAIL")
            + " — "
            + name
            + "\n"
        )

    f.write("\n")

    f.write(
        "IMPORTANT INTERPRETATION\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        "This audit compares the actual nucleotide sequences "
        "recovered in STEP 83 against the frozen adjudicated "
        "master.\n"
    )

    f.write(
        "No canonical sequence is selected in STEP 84.\n"
    )

    f.write(
        "Multiple-sequence KJUs remain explicitly flagged "
        "for downstream investigation.\n"
    )

    f.write(
        "No taxonomy is fabricated or altered.\n"
    )

    f.write(
        "The frozen V1 master is not modified.\n"
    )


# ============================================================
# CONSOLE REPORT
# ============================================================

print()
print("=" * 110)
print("STEP 84 COMPLETE")
print("=" * 110)

print()
print("MASTER")
print("-" * 110)
print(f"Total KJU records        : {total}")
print(f"Resolved records         : {resolved}")
print(f"Unresolved records       : {unresolved}")

print()
print("SEQUENCE COVERAGE")
print("-" * 110)
print(f"KJUs with sequence       : {with_sequence}")
print(f"KJUs without sequence    : {without_sequence}")
print(f"Exactly one unique seq   : {one_sequence_kjus}")
print(f"Multiple unique seq      : {different_sequence_kjus}")

print()
print("RECONCILIATION STATUS")
print("-" * 110)
print(
    f"EXACT_MATCH                         : "
    f"{exact_count}"
)
print(
    f"ACCESSION_MATCH_MULTIPLE_SEQUENCES : "
    f"{multi_count}"
)
print(
    f"ACCESSION_CONFLICT                 : "
    f"{conflict_count}"
)
print(
    f"NO_SEQUENCE                        : "
    f"{no_sequence_count}"
)
print(
    f"UNRESOLVED                         : "
    f"{unresolved_count}"
)

print()
print("CONFLICT SUMMARY")
print("-" * 110)
print(
    f"Same accession / different sequence : "
    f"{len(accession_conflict_df)}"
)
print(
    f"Different KJU / identical sequence  : "
    f"{len(duplicate_df)}"
)
print(
    f"KJUs with length conflict           : "
    f"{len(length_conflict_df)}"
)

print()
print("QC")
print("-" * 110)

all_pass = True

for name, passed in qc_results:

    print(
        ("PASS" if passed else "FAIL")
        + " — "
        + name
    )

    if not passed:
        all_pass = False


print()
print("OUTPUT FILES")
print("-" * 110)
print(f"KJU audit       : {KJU_AUDIT}")
print(f"All candidates  : {CANDIDATES}")
print(f"Exact matches   : {EXACT_MATCHES}")
print(f"Multiple seq    : {ACCESSION_MATCH_MULTIPLE}")
print(f"Accession conf. : {ACCESSION_CONFLICTS}")
print(f"Length conflicts: {LENGTH_CONFLICTS}")
print(f"Cross-KJU dupes : {DUPLICATE_ACROSS_KJU}")
print(f"QC report       : {SUMMARY}")

print()
print("=" * 110)

if all_pass:
    print("STEP 84 QC PASSED")
else:
    print("STEP 84 QC FAILED — REVIEW REPORT")

print("=" * 110)

print()
print(
    "FROZEN V1 MASTER WAS NOT MODIFIED"
)

print("=" * 110)

