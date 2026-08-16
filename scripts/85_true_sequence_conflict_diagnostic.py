#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
from collections import defaultdict


# ============================================================
# STEP 85
# TRUE SEQUENCE CONFLICT + LENGTH DIAGNOSTIC
#
# INPUT:
#   STEP 76 frozen master
#   STEP 83 true nucleotide sequence records
#   STEP 84 reconciliation outputs
#
# PURPOSE:
#   Diagnose the 109 KJUs reported by STEP 84 as having
#   length conflicts.
#
# IMPORTANT:
#   This is READ-ONLY.
#   The frozen V1 master is NEVER modified.
#
# NO canonical sequence is selected in this step.
# NO taxonomy is changed.
# NO accession is changed.
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

FREEZE_DIR = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"

MASTER = (
    FREEZE_DIR
    / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
)

TRUE_DIR = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"

TRUE_RECORDS = TRUE_DIR / "STEP83_TRUE_SEQUENCE_RECORDS.csv"

STEP84_DIR = ROOT / "STEP84_TRUE_SEQUENCE_MASTER_RECONCILIATION"

STEP84_KJU = (
    STEP84_DIR
    / "STEP84_KJU_TRUE_SEQUENCE_RECONCILIATION.csv"
)

OUT = ROOT / "STEP85_TRUE_SEQUENCE_CONFLICT_DIAGNOSTIC"
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# OUTPUTS
# ============================================================

KJU_AUDIT = OUT / "STEP85_KJU_CONFLICT_AUDIT.csv"

LENGTH_AUDIT = OUT / "STEP85_LENGTH_CONFLICT_AUDIT.csv"

ACCESSION_AUDIT = OUT / "STEP85_ACCESSION_CONFLICT_AUDIT.csv"

CANDIDATE_AUDIT = OUT / "STEP85_SEQUENCE_CANDIDATE_AUDIT.csv"

SAME_ACCESSION_DIFF_LENGTH = (
    OUT / "STEP85_SAME_ACCESSION_DIFFERENT_LENGTH.csv"
)

SAME_ACCESSION_DIFF_SEQUENCE = (
    OUT / "STEP85_SAME_ACCESSION_DIFFERENT_SEQUENCE.csv"
)

CROSS_KJU_IDENTICAL = (
    OUT / "STEP85_CROSS_KJU_IDENTICAL_SEQUENCES.csv"
)

SUMMARY = OUT / "STEP85_TRUE_SEQUENCE_CONFLICT_DIAGNOSTIC_REPORT.txt"


# ============================================================
# HELPERS
# ============================================================

def clean_string(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def clean_int(x):
    if pd.isna(x):
        return None

    s = str(x).strip()

    if not s:
        return None

    try:
        return int(float(s))
    except Exception:
        return None


def normalize_accession(x):
    s = clean_string(x)

    if not s:
        return ""

    s = s.strip()

    # Remove whitespace
    s = re.sub(r"\s+", "", s)

    # Normalize common accession.version form
    return s.upper()


def normalize_sequence(x):
    """
    Defensive sequence normalization.

    Removes:
      whitespace
      FASTA headers
      non-IUPAC characters

    Keeps:
      A C G T U N
      ambiguity codes
    """

    s = clean_string(x)

    if not s:
        return ""

    # Remove FASTA headers if present
    lines = s.splitlines()

    if len(lines) > 1:
        lines = [
            line.strip()
            for line in lines
            if not line.strip().startswith(">")
        ]
        s = "".join(lines)

    # Remove whitespace
    s = re.sub(r"\s+", "", s)

    # Uppercase
    s = s.upper()

    # Keep nucleotide / ambiguity symbols
    s = re.sub(
        r"[^ACGTUNRYSWKMBDHV\-]",
        "",
        s
    )

    # Remove alignment gaps
    s = s.replace("-", "")

    return s


def sha256_sequence(seq):
    if not seq:
        return ""

    return hashlib.sha256(
        seq.encode("utf-8")
    ).hexdigest()


def unique_sorted(values):
    values = {
        clean_string(v)
        for v in values
        if clean_string(v)
    }
    return sorted(values)


def join_values(values):
    vals = unique_sorted(values)
    return ";".join(vals)


def status_for_candidate(
    accepted_accession,
    accepted_length,
    candidate_accession,
    candidate_length,
    accepted_sha,
    candidate_sha,
):
    accession_match = (
        bool(accepted_accession)
        and bool(candidate_accession)
        and
        normalize_accession(accepted_accession)
        ==
        normalize_accession(candidate_accession)
    )

    length_match = (
        accepted_length is not None
        and candidate_length is not None
        and accepted_length == candidate_length
    )

    sha_match = (
        bool(accepted_sha)
        and bool(candidate_sha)
        and accepted_sha == candidate_sha
    )

    if accession_match and length_match and sha_match:
        return "EXACT_ACCESSION_LENGTH_SEQUENCE_MATCH"

    if accession_match and length_match:
        return "ACCESSION_LENGTH_MATCH"

    if accession_match and not length_match:
        return "ACCESSION_MATCH_LENGTH_CONFLICT"

    if not accession_match and length_match:
        return "LENGTH_MATCH_ACCESSION_CONFLICT"

    if accession_match:
        return "ACCESSION_MATCH"

    if length_match:
        return "LENGTH_MATCH"

    return "ACCESSION_AND_LENGTH_CONFLICT"


# ============================================================
# HEADER
# ============================================================

print("=" * 110)
print("STEP 85: TRUE SEQUENCE CONFLICT + LENGTH DIAGNOSTIC")
print("=" * 110)


# ============================================================
# INPUT VALIDATION
# ============================================================

if not MASTER.exists():
    raise FileNotFoundError(
        f"Frozen master not found:\n{MASTER}"
    )

if not TRUE_RECORDS.exists():
    raise FileNotFoundError(
        f"True sequence records not found:\n{TRUE_RECORDS}"
    )

print()
print("INPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"True records  : {TRUE_RECORDS}")

if STEP84_KJU.exists():
    print(f"STEP84 audit  : {STEP84_KJU}")
else:
    print("STEP84 audit  : NOT FOUND — continuing independently")


# ============================================================
# LOAD
# ============================================================

master = pd.read_csv(
    MASTER,
    dtype=str,
    keep_default_na=False
)

raw = pd.read_csv(
    TRUE_RECORDS,
    dtype=str,
    keep_default_na=False
)

print()
print("INPUT SIZES")
print("-" * 110)
print(f"Master rows       : {len(master)}")
print(f"True sequence rows: {len(raw)}")


# ============================================================
# COLUMN DETECTION
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
        raise KeyError(
            f"Could not detect required column. "
            f"Tried: {candidates}\n"
            f"Available: {list(df.columns)}"
        )

    return None


MASTER_KJU = detect_column(
    master,
    ["KJU_ID", "KJU", "KJU_ID_FINAL"]
)

MASTER_ACCESSION = detect_column(
    master,
    [
        "Accepted_accession",
        "accepted_accession",
        "Accession"
    ]
)

MASTER_LENGTH = detect_column(
    master,
    [
        "Accepted_length",
        "accepted_length",
        "Length"
    ]
)

MASTER_TAXONOMY = detect_column(
    master,
    [
        "Accepted_taxonomy",
        "accepted_taxonomy",
        "Taxonomy"
    ]
)

MASTER_IDENTITY = detect_column(
    master,
    [
        "Biological_identity",
        "biological_identity"
    ]
)

RAW_KJU = detect_column(
    raw,
    ["KJU_ID", "KJU"]
)

RAW_ACCESSION = detect_column(
    raw,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "Accession"
    ],
    required=False
)

RAW_LENGTH = detect_column(
    raw,
    [
        "Sequence_length",
        "sequence_length",
        "Length"
    ],
    required=False
)

RAW_SHA = detect_column(
    raw,
    [
        "Sequence_SHA256",
        "sequence_sha256",
        "SHA256"
    ],
    required=False
)

RAW_TYPE = detect_column(
    raw,
    [
        "Sequence_type",
        "Sequence_type_detected"
    ],
    required=False
)

RAW_RECORD = detect_column(
    raw,
    [
        "Record_ID",
        "record_id"
    ],
    required=False
)

RAW_FILE = detect_column(
    raw,
    [
        "Sequence_file",
        "Sequence_file_name"
    ],
    required=False
)

print()
print("DETECTED COLUMNS")
print("-" * 110)
print(f"Master KJU       : {MASTER_KJU}")
print(f"Master accession : {MASTER_ACCESSION}")
print(f"Master length    : {MASTER_LENGTH}")
print(f"Master taxonomy  : {MASTER_TAXONOMY}")
print(f"Master identity  : {MASTER_IDENTITY}")
print(f"Sequence KJU     : {RAW_KJU}")
print(f"Sequence accession: {RAW_ACCESSION}")
print(f"Sequence length  : {RAW_LENGTH}")
print(f"Sequence SHA256  : {RAW_SHA}")
print(f"Sequence type    : {RAW_TYPE}")
print(f"Record ID        : {RAW_RECORD}")


# ============================================================
# NORMALIZE MASTER
# ============================================================

master["_KJU"] = master[MASTER_KJU].map(clean_string)
master["_accepted_accession"] = (
    master[MASTER_ACCESSION].map(normalize_accession)
)

master["_accepted_length"] = (
    master[MASTER_LENGTH].map(clean_int)
)

master["_taxonomy"] = (
    master[MASTER_TAXONOMY].map(clean_string)
)

master["_biological_identity"] = (
    master[MASTER_IDENTITY].map(clean_string)
)


# ============================================================
# NORMALIZE TRUE SEQUENCES
# ============================================================

raw["_KJU"] = raw[RAW_KJU].map(clean_string)

if RAW_ACCESSION:
    raw["_candidate_accession"] = (
        raw[RAW_ACCESSION].map(normalize_accession)
    )
else:
    raw["_candidate_accession"] = ""

if RAW_LENGTH:
    raw["_candidate_length"] = (
        raw[RAW_LENGTH].map(clean_int)
    )
else:
    raw["_candidate_length"] = None

if RAW_SHA:
    raw["_candidate_sha"] = (
        raw[RAW_SHA]
        .map(clean_string)
        .str.lower()
    )
else:
    raw["_candidate_sha"] = ""

if RAW_TYPE:
    raw["_sequence_type"] = (
        raw[RAW_TYPE].map(clean_string)
    )
else:
    raw["_sequence_type"] = ""

if RAW_RECORD:
    raw["_record_id"] = (
        raw[RAW_RECORD].map(clean_string)
    )
else:
    raw["_record_id"] = ""

if RAW_FILE:
    raw["_sequence_file"] = (
        raw[RAW_FILE].map(clean_string)
    )
else:
    raw["_sequence_file"] = ""


# ============================================================
# REMOVE EMPTY KJU
# ============================================================

raw = raw[
    raw["_KJU"].astype(str).str.strip() != ""
].copy()

print()
print("SEQUENCE DATA")
print("-" * 110)
print(f"Valid sequence rows : {len(raw)}")
print(
    f"KJUs represented    : "
    f"{raw['_KJU'].nunique()}"
)


# ============================================================
# MASTER LOOKUP
# ============================================================

master_lookup = {}

for _, row in master.iterrows():

    kju = row["_KJU"]

    master_lookup[kju] = {
        "accepted_accession":
            row["_accepted_accession"],

        "accepted_length":
            row["_accepted_length"],

        "taxonomy":
            row["_taxonomy"],

        "biological_identity":
            row["_biological_identity"]
    }


# ============================================================
# GLOBAL SEQUENCE HASH NORMALIZATION
# ============================================================

raw["_sha_norm"] = raw["_candidate_sha"].map(
    lambda x: clean_string(x).lower()
)

# If SHA is absent, derive from available sequence columns
# where possible. STEP83 should already contain valid SHA256.
missing_sha = raw["_sha_norm"] == ""

if missing_sha.any():

    print()
    print(
        f"WARNING: {missing_sha.sum()} rows have blank SHA256."
    )

    raw.loc[missing_sha, "_sha_norm"] = (
        raw.loc[missing_sha]
        .apply(
            lambda r: hashlib.sha256(
                (
                    str(r.get("Sequence", ""))
                    .strip()
                    .upper()
                    .encode("utf-8")
                )
            ).hexdigest()
            if str(r.get("Sequence", "")).strip()
            else "",
            axis=1
        )
    )


# ============================================================
# BUILD CANDIDATE AUDIT
# ============================================================

candidate_rows = []

for _, row in raw.iterrows():

    kju = row["_KJU"]

    m = master_lookup.get(
        kju,
        {
            "accepted_accession": "",
            "accepted_length": None,
            "taxonomy": "",
            "biological_identity": ""
        }
    )

    accepted_accession = m["accepted_accession"]
    accepted_length = m["accepted_length"]

    candidate_accession = row["_candidate_accession"]
    candidate_length = row["_candidate_length"]

    accession_match = (
        bool(accepted_accession)
        and bool(candidate_accession)
        and
        accepted_accession
        ==
        candidate_accession
    )

    length_match = (
        accepted_length is not None
        and candidate_length is not None
        and
        accepted_length
        ==
        candidate_length
    )

    candidate_rows.append({
        "KJU_ID": kju,
        "Accepted_accession":
            accepted_accession,
        "Accepted_length":
            accepted_length,
        "Accepted_taxonomy":
            m["taxonomy"],
        "Biological_identity":
            m["biological_identity"],
        "Candidate_accession":
            candidate_accession,
        "Candidate_length":
            candidate_length,
        "Candidate_SHA256":
            row["_sha_norm"],
        "Sequence_type":
            row["_sequence_type"],
        "Record_ID":
            row["_record_id"],
        "Sequence_file":
            row["_sequence_file"],
        "Accession_match":
            accession_match,
        "Length_match":
            length_match,
        "Length_difference":
            (
                candidate_length - accepted_length
                if (
                    candidate_length is not None
                    and accepted_length is not None
                )
                else None
            ),
        "Candidate_status":
            status_for_candidate(
                accepted_accession,
                accepted_length,
                candidate_accession,
                candidate_length,
                "",
                row["_sha_norm"]
            )
    })


candidate_df = pd.DataFrame(candidate_rows)


# ============================================================
# KJU-LEVEL AUDIT
# ============================================================

audit_rows = []

all_kjus = sorted(
    set(master["_KJU"])
)

for kju in all_kjus:

    m = master_lookup[kju]

    subset = candidate_df[
        candidate_df["KJU_ID"] == kju
    ]

    if subset.empty:

        audit_rows.append({
            "KJU_ID": kju,
            "Accepted_accession":
                m["accepted_accession"],
            "Accepted_length":
                m["accepted_length"],
            "Accepted_taxonomy":
                m["taxonomy"],
            "Biological_identity":
                m["biological_identity"],
            "Sequence_rows": 0,
            "Unique_sequences": 0,
            "Candidate_accessions": "",
            "Candidate_lengths": "",
            "Candidate_SHA256": "",
            "Accession_matches": 0,
            "Length_matches": 0,
            "Length_conflicts": 0,
            "Status": (
                "UNRESOLVED"
                if not m["accepted_accession"]
                and not m["taxonomy"]
                else "NO_SEQUENCE"
            )
        })

        continue

    unique_sha = sorted(
        set(
            subset["Candidate_SHA256"]
            .astype(str)
            .str.strip()
        )
        - {""}
    )

    candidate_accessions = sorted(
        set(
            subset["Candidate_accession"]
            .astype(str)
            .str.strip()
        )
        - {""}
    )

    candidate_lengths = sorted(
        set(
            pd.to_numeric(
                subset["Candidate_length"],
                errors="coerce"
            )
            .dropna()
            .astype(int)
            .tolist()
        )
    )

    accession_matches = int(
        subset["Accession_match"].sum()
    )

    length_matches = int(
        subset["Length_match"].sum()
    )

    length_conflicts = int(
        (
            subset["Length_match"] == False
        ).sum()
    )

    if (
        len(unique_sha) == 1
        and accession_matches > 0
        and length_matches > 0
    ):
        status = "EXACT_SINGLE_SEQUENCE_MATCH"

    elif (
        accession_matches > 0
        and length_matches > 0
        and len(unique_sha) > 1
    ):
        status = "MULTIPLE_SEQUENCES_WITH_MATCHING_CANDIDATE"

    elif (
        accession_matches > 0
        and length_matches == 0
    ):
        status = "ACCESSION_MATCH_LENGTH_CONFLICT"

    elif (
        accession_matches == 0
        and length_matches > 0
    ):
        status = "LENGTH_MATCH_ACCESSION_CONFLICT"

    elif (
        accession_matches == 0
        and length_matches == 0
    ):
        status = "ACCESSION_AND_LENGTH_CONFLICT"

    else:
        status = "REVIEW"

    audit_rows.append({
        "KJU_ID": kju,
        "Accepted_accession":
            m["accepted_accession"],
        "Accepted_length":
            m["accepted_length"],
        "Accepted_taxonomy":
            m["taxonomy"],
        "Biological_identity":
            m["biological_identity"],
        "Sequence_rows":
            len(subset),
        "Unique_sequences":
            len(unique_sha),
        "Candidate_accessions":
            ";".join(candidate_accessions),
        "Candidate_lengths":
            ";".join(
                str(x)
                for x in candidate_lengths
            ),
        "Candidate_SHA256":
            ";".join(unique_sha),
        "Accession_matches":
            accession_matches,
        "Length_matches":
            length_matches,
        "Length_conflicts":
            length_conflicts,
        "Status":
            status
    })


kju_audit = pd.DataFrame(audit_rows)


# ============================================================
# LENGTH CONFLICT TABLE
# ============================================================

length_conflict = candidate_df[
    candidate_df["Length_match"] == False
].copy()

length_conflict = length_conflict[
    length_conflict["Accepted_length"].notna()
    &
    length_conflict["Candidate_length"].notna()
].copy()

length_conflict["Absolute_length_difference"] = (
    length_conflict["Length_difference"]
    .abs()
)

length_conflict = length_conflict.sort_values(
    [
        "KJU_ID",
        "Absolute_length_difference"
    ],
    ascending=[True, False]
)


# ============================================================
# ACCESSION CONFLICT TABLE
# ============================================================

accession_conflict = candidate_df[
    (
        candidate_df["Accepted_accession"].astype(str).str.strip()
        != ""
    )
    &
    (
        candidate_df["Candidate_accession"].astype(str).str.strip()
        != ""
    )
    &
    (
        candidate_df["Accession_match"] == False
    )
].copy()


# ============================================================
# SAME ACCESSION / DIFFERENT LENGTH
# ============================================================

same_acc_diff_length = candidate_df[
    (
        candidate_df["Accepted_accession"].astype(str).str.strip()
        != ""
    )
    &
    (
        candidate_df["Candidate_accession"]
        ==
        candidate_df["Accepted_accession"]
    )
    &
    (
        candidate_df["Length_match"] == False
    )
].copy()


# ============================================================
# SAME ACCESSION / DIFFERENT SEQUENCE
# ============================================================

same_accession_groups = defaultdict(list)

for _, row in candidate_df.iterrows():

    acc = row["Candidate_accession"]
    sha = row["Candidate_SHA256"]

    if acc and sha:
        same_accession_groups[acc].append(
            (
                row["KJU_ID"],
                sha,
                row["Candidate_length"]
            )
        )


same_acc_diff_seq_rows = []

for acc, items in same_accession_groups.items():

    unique_sha = sorted(
        set(
            item[1]
            for item in items
            if item[1]
        )
    )

    if len(unique_sha) > 1:

        for kju, sha, length in items:

            same_acc_diff_seq_rows.append({
                "Accession": acc,
                "KJU_ID": kju,
                "SHA256": sha,
                "Length": length,
                "Unique_sequences_for_accession":
                    len(unique_sha)
            })


same_acc_diff_seq = pd.DataFrame(
    same_acc_diff_seq_rows
)


# ============================================================
# CROSS-KJU IDENTICAL SEQUENCES
# ============================================================

sha_groups = defaultdict(list)

for _, row in candidate_df.iterrows():

    sha = row["Candidate_SHA256"]

    if sha:
        sha_groups[sha].append(
            (
                row["KJU_ID"],
                row["Candidate_accession"],
                row["Candidate_length"]
            )
        )


cross_kju_rows = []

for sha, items in sha_groups.items():

    kju_set = sorted(
        set(
            item[0]
            for item in items
        )
    )

    if len(kju_set) > 1:

        for kju, acc, length in items:

            cross_kju_rows.append({
                "SHA256": sha,
                "KJU_ID": kju,
                "Candidate_accession": acc,
                "Candidate_length": length,
                "KJUs_sharing_sequence":
                    ";".join(kju_set),
                "Number_of_KJUs":
                    len(kju_set)
            })


cross_kju_df = pd.DataFrame(
    cross_kju_rows
)


# ============================================================
# SAVE OUTPUTS
# ============================================================

kju_audit.to_csv(
    KJU_AUDIT,
    index=False
)

length_conflict.to_csv(
    LENGTH_AUDIT,
    index=False
)

accession_conflict.to_csv(
    ACCESSION_AUDIT,
    index=False
)

candidate_df.to_csv(
    CANDIDATE_AUDIT,
    index=False
)

same_acc_diff_length.to_csv(
    SAME_ACCESSION_DIFF_LENGTH,
    index=False
)

same_acc_diff_seq.to_csv(
    SAME_ACCESSION_DIFF_SEQUENCE,
    index=False
)

cross_kju_df.to_csv(
    CROSS_KJU_IDENTICAL,
    index=False
)


# ============================================================
# SUMMARY COUNTS
# ============================================================

n_master = len(master)

resolved = int(
    (
        master["_accepted_accession"].astype(str).str.strip()
        != ""
    ).sum()
)

unresolved = n_master - resolved

with_sequence = int(
    (
        kju_audit["Sequence_rows"] > 0
    ).sum()
)

without_sequence = n_master - with_sequence

exact_single = int(
    (
        kju_audit["Status"]
        ==
        "EXACT_SINGLE_SEQUENCE_MATCH"
    ).sum()
)

accession_length_conflict = int(
    (
        kju_audit["Status"]
        ==
        "ACCESSION_MATCH_LENGTH_CONFLICT"
    ).sum()
)

length_accession_conflict = int(
    (
        kju_audit["Status"]
        ==
        "LENGTH_MATCH_ACCESSION_CONFLICT"
    ).sum()
)

both_conflict = int(
    (
        kju_audit["Status"]
        ==
        "ACCESSION_AND_LENGTH_CONFLICT"
    ).sum()
)

multiple_with_match = int(
    (
        kju_audit["Status"]
        ==
        "MULTIPLE_SEQUENCES_WITH_MATCHING_CANDIDATE"
    ).sum()
)

unique_sequences_per_kju = (
    kju_audit["Unique_sequences"]
)

exactly_one_unique = int(
    (
        unique_sequences_per_kju == 1
    ).sum()
)

multiple_unique = int(
    (
        unique_sequences_per_kju > 1
    ).sum()
)

same_acc_diff_seq_count = (
    len(same_acc_diff_seq)
)

cross_kju_count = (
    len(cross_kju_df)
)


# ============================================================
# UNRESOLVED QC
# ============================================================

unresolved_kjus = set(
    master[
        master["_accepted_accession"].astype(str).str.strip() == ""
    ]["_KJU"]
)

unresolved_with_sequence = set(
    kju_audit[
        kju_audit["KJU_ID"].isin(unresolved_kjus)
        &
        (kju_audit["Sequence_rows"] > 0)
    ]["KJU_ID"]
)


# ============================================================
# REPORT
# ============================================================

lines = []

lines.append("=" * 110)
lines.append(
    "STEP 85: TRUE SEQUENCE CONFLICT + LENGTH DIAGNOSTIC"
)
lines.append("=" * 110)

lines.append("")
lines.append("INPUTS")
lines.append("-" * 110)
lines.append(f"Frozen master       : {MASTER}")
lines.append(f"True sequence file  : {TRUE_RECORDS}")

lines.append("")
lines.append("INPUT SIZES")
lines.append("-" * 110)
lines.append(f"Master records      : {n_master}")
lines.append(f"True sequence rows  : {len(raw)}")
lines.append(
    f"KJUs with sequences : {with_sequence}"
)

lines.append("")
lines.append("SEQUENCE COVERAGE")
lines.append("-" * 110)
lines.append(
    f"KJUs with sequence       : {with_sequence}"
)
lines.append(
    f"KJUs without sequence    : {without_sequence}"
)
lines.append(
    f"Exactly one unique seq   : {exactly_one_unique}"
)
lines.append(
    f"Multiple unique seq      : {multiple_unique}"
)

lines.append("")
lines.append("RECONCILIATION CLASSIFICATION")
lines.append("-" * 110)
lines.append(
    f"Exact single sequence match              : {exact_single}"
)
lines.append(
    f"Multiple sequences with matching candidate: {multiple_with_match}"
)
lines.append(
    f"Accession match + length conflict        : {accession_length_conflict}"
)
lines.append(
    f"Length match + accession conflict        : {length_accession_conflict}"
)
lines.append(
    f"Accession + length conflict              : {both_conflict}"
)

lines.append("")
lines.append("CONFLICT COUNTS")
lines.append("-" * 110)
lines.append(
    f"Candidate rows with length conflict      : {len(length_conflict)}"
)
lines.append(
    f"Candidate rows with accession conflict   : {len(accession_conflict)}"
)
lines.append(
    f"Same accession / different length rows   : {len(same_acc_diff_length)}"
)
lines.append(
    f"Same accession / different sequence rows : {same_acc_diff_seq_count}"
)
lines.append(
    f"Cross-KJU identical sequence rows        : {cross_kju_count}"
)

lines.append("")
lines.append("UNRESOLVED QC")
lines.append("-" * 110)

if unresolved_with_sequence:
    lines.append(
        "FAIL — unresolved KJUs contain actual sequences:"
    )
    for kju in sorted(unresolved_with_sequence):
        lines.append(f"  {kju}")
else:
    lines.append(
        "PASS — no unresolved KJU has an actual sequence"
    )

lines.append("")
lines.append("INTERPRETATION")
lines.append("-" * 110)

if multiple_unique > 0:
    lines.append(
        "REVIEW REQUIRED — multiple distinct nucleotide sequences "
        "exist for one or more KJUs."
    )

if same_acc_diff_seq_count > 0:
    lines.append(
        "REVIEW REQUIRED — at least one accession is associated "
        "with multiple distinct sequences."
    )

if len(same_acc_diff_length) > 0:
    lines.append(
        "REVIEW REQUIRED — accepted accession matches exist with "
        "different sequence lengths."
    )

if exact_single > 0:
    lines.append(
        f"Potential deterministic matches identified: {exact_single}"
    )

lines.append("")
lines.append("IMPORTANT")
lines.append("-" * 110)
lines.append(
    "STEP 85 DOES NOT SELECT OR MODIFY ANY CANONICAL SEQUENCE."
)
lines.append(
    "STEP 85 DOES NOT MODIFY TAXONOMY."
)
lines.append(
    "STEP 85 DOES NOT MODIFY THE FROZEN V1 MASTER."
)

lines.append("")
lines.append("OUTPUT FILES")
lines.append("-" * 110)
lines.append(f"KJU audit            : {KJU_AUDIT}")
lines.append(f"Length audit         : {LENGTH_AUDIT}")
lines.append(f"Accession audit      : {ACCESSION_AUDIT}")
lines.append(f"Candidate audit      : {CANDIDATE_AUDIT}")
lines.append(
    f"Same accession/length: {SAME_ACCESSION_DIFF_LENGTH}"
)
lines.append(
    f"Same accession/seq   : {SAME_ACCESSION_DIFF_SEQUENCE}"
)
lines.append(
    f"Cross-KJU identical  : {CROSS_KJU_IDENTICAL}"
)

lines.append("")
lines.append("=" * 110)
lines.append("STEP 85 COMPLETE")
lines.append("FROZEN V1 MASTER WAS NOT MODIFIED")
lines.append("=" * 110)

SUMMARY.write_text(
    "\n".join(lines),
    encoding="utf-8"
)


# ============================================================
# CONSOLE OUTPUT
# ============================================================

print()
print("=" * 110)
print("STEP 85 COMPLETE")
print("=" * 110)

print()
print("MASTER")
print("-" * 110)
print(f"Total KJU records        : {n_master}")
print(f"Resolved records         : {resolved}")
print(f"Unresolved records       : {unresolved}")

print()
print("SEQUENCE COVERAGE")
print("-" * 110)
print(f"KJUs with sequence       : {with_sequence}")
print(f"KJUs without sequence    : {without_sequence}")
print(f"Exactly one unique seq   : {exactly_one_unique}")
print(f"Multiple unique seq      : {multiple_unique}")

print()
print("RECONCILIATION")
print("-" * 110)
print(
    f"Exact single match       : {exact_single}"
)
print(
    f"Multiple with match      : {multiple_with_match}"
)
print(
    f"Accession + length conflict : "
    f"{accession_length_conflict}"
)
print(
    f"Length + accession conflict : "
    f"{length_accession_conflict}"
)
print(
    f"Both accession/length conflict: "
    f"{both_conflict}"
)

print()
print("CONFLICTS")
print("-" * 110)
print(
    f"Length-conflict rows     : {len(length_conflict)}"
)
print(
    f"Accession-conflict rows  : {len(accession_conflict)}"
)
print(
    f"Same accession/diff length: "
    f"{len(same_acc_diff_length)}"
)
print(
    f"Same accession/diff seq  : "
    f"{same_acc_diff_seq_count}"
)
print(
    f"Cross-KJU identical seq : "
    f"{cross_kju_count}"
)

print()
print("UNRESOLVED QC")
print("-" * 110)

if unresolved_with_sequence:
    print(
        "FAIL — unresolved KJUs contain actual sequences"
    )
else:
    print(
        "PASS — no unresolved KJU has an actual sequence"
    )

print()
print("OUTPUT DIRECTORY")
print("-" * 110)
print(OUT)

print()
print("=" * 110)
print("FROZEN V1 MASTER WAS NOT MODIFIED")
print("=" * 110)

