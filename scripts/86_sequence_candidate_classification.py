#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
import math


# ============================================================
# STEP 86
# SEQUENCE CANDIDATE CLASSIFICATION + REGION CONSISTENCY AUDIT
#
# INPUT:
#   STEP 76 frozen master
#   STEP 83 true nucleotide sequence records
#   STEP 84 reconciliation outputs
#
# PURPOSE:
#   1. Classify every REAL sequence candidate.
#   2. Handle missing/NaN accepted lengths safely.
#   3. Compare sequence candidates against:
#        - accepted accession
#        - accepted length
#        - accepted taxonomy
#        - biological identity
#   4. Detect whether multiple sequences represent:
#        - exact duplicates
#        - different lengths
#        - same accession / different sequences
#        - different accession / same sequence
#        - likely region variants
#   5. Produce a deterministic audit.
#
# IMPORTANT:
#   - READ-ONLY with respect to frozen master.
#   - NO taxonomy changes.
#   - NO accession changes.
#   - NO canonical sequence is written into the master.
#   - Missing accepted lengths are NOT treated as zero.
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

FREEZE = (
    ROOT
    / "FINAL_KJU_MASTER_FREEZE_V1"
)

STEP83 = (
    ROOT
    / "STEP83_TRUE_SEQUENCE_INVENTORY"
)

STEP84 = (
    ROOT
    / "STEP84_TRUE_SEQUENCE_MASTER_RECONCILIATION"
)

OUT = (
    ROOT
    / "STEP86_SEQUENCE_CANDIDATE_CLASSIFICATION"
)

OUT.mkdir(parents=True, exist_ok=True)


MASTER = (
    FREEZE
    / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
)

TRUE_RECORDS = (
    STEP83
    / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
)

STEP84_AUDIT = (
    STEP84
    / "STEP84_KJU_TRUE_SEQUENCE_RECONCILIATION.csv"
)


KJU_AUDIT = (
    OUT
    / "STEP86_KJU_CANDIDATE_AUDIT.csv"
)

CANDIDATES = (
    OUT
    / "STEP86_SEQUENCE_CANDIDATE_CLASSIFICATION.csv"
)

REVIEW = (
    OUT
    / "STEP86_MANUAL_REVIEW_REQUIRED.csv"
)

REGION_AUDIT = (
    OUT
    / "STEP86_REGION_CONSISTENCY_AUDIT.csv"
)

DUPLICATES = (
    OUT
    / "STEP86_DUPLICATE_SEQUENCE_GROUPS.csv"
)

ACCESSION_CONFLICTS = (
    OUT
    / "STEP86_ACCESSION_SEQUENCE_CONFLICTS.csv"
)

REPORT = (
    OUT
    / "STEP86_SEQUENCE_CANDIDATE_CLASSIFICATION_REPORT.txt"
)


# ============================================================
# HELPERS
# ============================================================

def clean_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def clean_kju(x):
    return clean_text(x).upper()


def clean_accession(x):
    x = clean_text(x)
    if not x:
        return ""
    return x.upper()


def safe_int(x):
    """
    Convert numeric-like value to int.
    Return None for NaN / blank / malformed values.
    """
    if pd.isna(x):
        return None

    s = str(x).strip()

    if not s:
        return None

    try:
        value = float(s)

        if math.isnan(value):
            return None

        return int(value)

    except Exception:
        return None


def normalize_sequence(seq):
    """
    Normalize nucleotide sequence safely.

    Removes:
      whitespace
      FASTA headers
      non-letter formatting

    Keeps:
      A C G T N and IUPAC ambiguity codes
    """

    if pd.isna(seq):
        return ""

    s = str(seq).strip()

    if not s:
        return ""

    # Remove FASTA headers if present
    lines = s.splitlines()

    if lines and lines[0].startswith(">"):
        lines = lines[1:]

    s = "".join(lines)

    # Remove whitespace
    s = re.sub(r"\s+", "", s)

    # Uppercase
    s = s.upper()

    # Keep nucleotide/IUPAC characters
    s = re.sub(r"[^ACGTUNRYSWKMBDHV]", "", s)

    return s


def sequence_sha256(seq):
    if not seq:
        return ""
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()


def is_valid_sequence(seq):
    if not seq:
        return False

    return bool(
        re.fullmatch(
            r"[ACGTUNRYSWKMBDHV]+",
            seq
        )
    )


def accession_match(candidate, accepted):
    candidate = clean_accession(candidate)
    accepted = clean_accession(accepted)

    if not candidate or not accepted:
        return False

    return candidate == accepted


def length_match(candidate_length, accepted_length):
    if candidate_length is None:
        return False

    if accepted_length is None:
        return False

    return candidate_length == accepted_length


def taxonomy_genus(taxonomy):
    taxonomy = clean_text(taxonomy)

    if not taxonomy:
        return ""

    return taxonomy.split()[0]


def likely_same_region(length_a, length_b):
    """
    Conservative region-consistency heuristic.

    This does NOT claim biological equivalence.

    Used only to flag candidates whose lengths are close enough
    that they may represent overlapping/trimmed versions of
    the same marker region.
    """

    if length_a is None or length_b is None:
        return False

    if length_a == 0 or length_b == 0:
        return False

    difference = abs(length_a - length_b)

    shorter = min(length_a, length_b)

    if shorter == 0:
        return False

    ratio = difference / shorter

    return ratio <= 0.15


# ============================================================
# START
# ============================================================

print("=" * 110)
print("STEP 86: SEQUENCE CANDIDATE CLASSIFICATION + REGION CONSISTENCY AUDIT")
print("=" * 110)


# ============================================================
# INPUT VALIDATION
# ============================================================

for path in [MASTER, TRUE_RECORDS, STEP84_AUDIT]:
    if not path.exists():
        raise FileNotFoundError(
            f"Required input file does not exist:\n{path}"
        )


print()
print("INPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"True records  : {TRUE_RECORDS}")
print(f"STEP84 audit  : {STEP84_AUDIT}")


# ============================================================
# LOAD
# ============================================================

master = pd.read_csv(MASTER, dtype=str, keep_default_na=False)

true = pd.read_csv(
    TRUE_RECORDS,
    dtype=str,
    keep_default_na=False
)

step84 = pd.read_csv(
    STEP84_AUDIT,
    dtype=str,
    keep_default_na=False
)


print()
print("INPUT SIZES")
print("-" * 110)
print(f"Master rows        : {len(master)}")
print(f"True sequence rows : {len(true)}")
print(f"STEP84 audit rows  : {len(step84)}")


# ============================================================
# COLUMN DETECTION
# ============================================================

def detect_column(df, candidates, label, required=True):

    lookup = {
        str(c).strip().lower(): c
        for c in df.columns
    }

    for candidate in candidates:

        key = candidate.lower()

        if key in lookup:
            return lookup[key]

    if required:
        raise RuntimeError(
            f"Could not detect {label} column.\n"
            f"Available columns:\n{list(df.columns)}"
        )

    return None


master_kju = detect_column(
    master,
    ["KJU_ID", "KJU", "KJU_ID_FINAL"],
    "master KJU"
)

master_accession = detect_column(
    master,
    [
        "Accepted_accession",
        "accepted_accession",
        "Accession"
    ],
    "master accession"
)

master_taxonomy = detect_column(
    master,
    [
        "Accepted_taxonomy",
        "accepted_taxonomy",
        "Taxonomy"
    ],
    "master taxonomy"
)

master_length = detect_column(
    master,
    [
        "Accepted_length",
        "accepted_length",
        "Length"
    ],
    "master length"
)

master_identity = detect_column(
    master,
    [
        "Biological_identity",
        "biological_identity"
    ],
    "master biological identity"
)


true_kju = detect_column(
    true,
    ["KJU_ID", "KJU", "KJU_ID_FINAL"],
    "sequence KJU"
)

true_accession = detect_column(
    true,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "Accession"
    ],
    "sequence accession"
)

true_length = detect_column(
    true,
    [
        "Sequence_length",
        "sequence_length",
        "Length"
    ],
    "sequence length"
)

true_sha = detect_column(
    true,
    [
        "Sequence_SHA256",
        "sequence_sha256",
        "SHA256"
    ],
    "sequence SHA256"
)

true_type = detect_column(
    true,
    [
        "Sequence_type",
        "Sequence_type_detected",
        "sequence_type"
    ],
    "sequence type",
    required=False
)

true_record = detect_column(
    true,
    [
        "Record_ID",
        "record_id",
        "Record"
    ],
    "record ID",
    required=False
)

true_sequence = detect_column(
    true,
    [
        "Sequence",
        "sequence",
        "Nucleotide_sequence",
        "Actual_sequence"
    ],
    "sequence field",
    required=False
)


print()
print("DETECTED COLUMNS")
print("-" * 110)
print(f"Master KJU       : {master_kju}")
print(f"Master accession : {master_accession}")
print(f"Master length    : {master_length}")
print(f"Master taxonomy  : {master_taxonomy}")
print(f"Master identity  : {master_identity}")
print(f"Sequence KJU     : {true_kju}")
print(f"Sequence accession: {true_accession}")
print(f"Sequence length  : {true_length}")
print(f"Sequence SHA256  : {true_sha}")
print(f"Sequence type    : {true_type}")
print(f"Record ID        : {true_record}")
print(f"Sequence field   : {true_sequence}")


# ============================================================
# NORMALIZE MASTER
# ============================================================

master["_KJU"] = master[master_kju].map(clean_kju)

master["_accepted_accession"] = master[
    master_accession
].map(clean_accession)

master["_accepted_length"] = master[
    master_length
].map(safe_int)

master["_taxonomy"] = master[
    master_taxonomy
].map(clean_text)

master["_identity"] = master[
    master_identity
].map(clean_text)


# ============================================================
# MASTER VALIDATION
# ============================================================

master_kjus = set(master["_KJU"])

if len(master) != 204:
    print(
        f"WARNING — expected 204 master rows, found {len(master)}"
    )

if len(master_kjus) != len(master):
    raise RuntimeError(
        "FATAL — duplicate KJU IDs detected in frozen master."
    )


# ============================================================
# NORMALIZE TRUE SEQUENCES
# ============================================================

true["_KJU"] = true[true_kju].map(clean_kju)

true["_candidate_accession"] = true[
    true_accession
].map(clean_accession)

true["_candidate_length"] = true[
    true_length
].map(safe_int)

if true_type:
    true["_sequence_type"] = true[
        true_type
    ].map(clean_text)
else:
    true["_sequence_type"] = ""


if true_record:
    true["_record"] = true[
        true_record
    ].map(clean_text)
else:
    true["_record"] = ""


# ============================================================
# NORMALIZE ACTUAL SEQUENCE
# ============================================================

if true_sequence:

    true["_sequence"] = true[
        true_sequence
    ].map(normalize_sequence)

else:

    true["_sequence"] = ""


# ============================================================
# RE-CALCULATE SHA256
# ============================================================

true["_calculated_sha256"] = true[
    "_sequence"
].map(sequence_sha256)


# ============================================================
# VALID SEQUENCE FILTER
# ============================================================

valid = true[
    true["_sequence"].map(is_valid_sequence)
].copy()


print()
print("SEQUENCE DATA")
print("-" * 110)
print(f"Total sequence rows       : {len(true)}")
print(f"Valid nucleotide rows     : {len(valid)}")
print(f"KJUs represented          : {valid['_KJU'].nunique()}")


# ============================================================
# MASTER LOOKUP
# ============================================================

master_lookup = (
    master
    .set_index("_KJU")
    .to_dict("index")
)


# ============================================================
# BUILD CANDIDATE CLASSIFICATION
# ============================================================

candidate_rows = []

for idx, row in valid.iterrows():

    kju = row["_KJU"]

    if kju not in master_lookup:
        status = "KJU_NOT_IN_MASTER"

        accepted_acc = ""
        accepted_len = None
        taxonomy = ""
        identity = ""

    else:

        m = master_lookup[kju]

        accepted_acc = m["_accepted_accession"]
        accepted_len = m["_accepted_length"]
        taxonomy = m["_taxonomy"]
        identity = m["_identity"]

        acc_match = accession_match(
            row["_candidate_accession"],
            accepted_acc
        )

        len_match = length_match(
            row["_candidate_length"],
            accepted_len
        )

        has_acc = bool(row["_candidate_accession"])
        has_accepted_acc = bool(accepted_acc)

        has_len = row["_candidate_length"] is not None
        has_accepted_len = accepted_len is not None

        if acc_match and len_match:
            status = "EXACT_ACCESSION_AND_LENGTH_MATCH"

        elif acc_match and not has_accepted_len:
            status = "ACCESSION_MATCH_ACCEPTED_LENGTH_MISSING"

        elif acc_match and has_accepted_len and not has_len:
            status = "ACCESSION_MATCH_CANDIDATE_LENGTH_MISSING"

        elif acc_match and has_len and has_accepted_len:
            status = "ACCESSION_MATCH_LENGTH_CONFLICT"

        elif (
            not has_acc
            and has_len
            and has_accepted_len
            and len_match
        ):
            status = "LENGTH_MATCH_ACCESSION_MISSING"

        elif (
            has_acc
            and has_accepted_acc
            and has_len
            and has_accepted_len
            and len_match
        ):
            status = "LENGTH_MATCH_ACCESSION_CONFLICT"

        elif (
            has_acc
            and has_accepted_acc
            and has_len
            and has_accepted_len
            and not len_match
        ):
            status = "ACCESSION_AND_LENGTH_CONFLICT"

        elif has_acc and has_accepted_acc:
            status = "ACCESSION_CONFLICT"

        elif has_len and has_accepted_len:
            status = "LENGTH_CONFLICT"

        else:
            status = "INSUFFICIENT_METADATA"

    candidate_rows.append(
        {
            "KJU_ID": kju,
            "Record_ID": row["_record"],
            "Candidate_accession": row["_candidate_accession"],
            "Candidate_length": row["_candidate_length"],
            "Candidate_SHA256": (
                row["_calculated_sha256"]
                or clean_text(row[true_sha])
            ),
            "Sequence_type": row["_sequence_type"],
            "Accepted_accession": accepted_acc,
            "Accepted_length": accepted_len,
            "Accepted_taxonomy": taxonomy,
            "Biological_identity": identity,
            "Accession_match": (
                bool(
                    accession_match(
                        row["_candidate_accession"],
                        accepted_acc
                    )
                )
            ),
            "Length_match": (
                bool(
                    length_match(
                        row["_candidate_length"],
                        accepted_len
                    )
                )
            ),
            "Classification": status,
        }
    )


candidate_df = pd.DataFrame(candidate_rows)


# ============================================================
# UNIQUE CANDIDATE COLLAPSE
# ============================================================

if not candidate_df.empty:

    candidate_df["Candidate_SHA256"] = (
        candidate_df["Candidate_SHA256"]
        .fillna("")
        .astype(str)
    )

    unique_candidates = (
        candidate_df
        .drop_duplicates(
            subset=[
                "KJU_ID",
                "Candidate_SHA256"
            ]
        )
        .copy()
    )

else:

    unique_candidates = candidate_df.copy()


# ============================================================
# PER-KJU ANALYSIS
# ============================================================

kju_rows = []

for kju in sorted(master_kjus):

    m = master_lookup[kju]

    sub = unique_candidates[
        unique_candidates["KJU_ID"] == kju
    ].copy()

    n_candidates = len(sub)

    unique_shas = sorted(
        set(
            x
            for x in sub["Candidate_SHA256"].tolist()
            if x
        )
    )

    unique_accessions = sorted(
        set(
            x
            for x in sub["Candidate_accession"].tolist()
            if x
        )
    )

    lengths = sorted(
        set(
            int(x)
            for x in sub["Candidate_length"].tolist()
            if pd.notna(x)
        )
    )

    accepted_acc = m["_accepted_accession"]
    accepted_len = m["_accepted_length"]

    accession_matches = [
        x
        for x in unique_accessions
        if accession_match(x, accepted_acc)
    ]

    length_matches = [
        x
        for x in lengths
        if length_match(x, accepted_len)
    ]

    exact_candidates = sub[
        (
            sub["Accession_match"] == True
        )
        &
        (
            sub["Length_match"] == True
        )
    ]

    exact_count = len(exact_candidates)

    if n_candidates == 0:

        status = (
            "UNRESOLVED"
            if not m["_taxonomy"]
            else "NO_SEQUENCE"
        )

    elif exact_count == 1:

        status = "AUTO_SUPPORTED_EXACT_MATCH"

    elif exact_count > 1:

        status = "MULTIPLE_EXACT_MATCH_CANDIDATES"

    elif accession_matches:

        status = "ACCESSION_MATCH_BUT_NO_EXACT_LENGTH"

    elif length_matches:

        status = "LENGTH_MATCH_BUT_ACCESSION_CONFLICT"

    else:

        status = "REVIEW_REQUIRED"


    if n_candidates == 1:
        candidate_count_class = "SINGLE_CANDIDATE"

    elif n_candidates > 1:
        candidate_count_class = "MULTIPLE_CANDIDATES"

    else:
        candidate_count_class = "NO_CANDIDATE"


    kju_rows.append(
        {
            "KJU_ID": kju,
            "Accepted_accession": accepted_acc,
            "Accepted_length": accepted_len,
            "Accepted_taxonomy": m["_taxonomy"],
            "Biological_identity": m["_identity"],
            "Candidate_count": n_candidates,
            "Unique_sequence_count": len(unique_shas),
            "Unique_accession_count": len(unique_accessions),
            "Candidate_lengths": ";".join(
                str(x) for x in lengths
            ),
            "Candidate_accessions": ";".join(
                unique_accessions
            ),
            "Exact_accession_length_candidates": exact_count,
            "Candidate_count_class": candidate_count_class,
            "Resolution_status": status,
        }
    )


kju_df = pd.DataFrame(kju_rows)


# ============================================================
# REGION CONSISTENCY AUDIT
# ============================================================

region_rows = []

for kju in sorted(master_kjus):

    sub = unique_candidates[
        unique_candidates["KJU_ID"] == kju
    ].copy()

    if sub.empty:
        continue

    lengths = sorted(
        set(
            int(x)
            for x in sub["Candidate_length"].tolist()
            if pd.notna(x)
        )
    )

    if len(lengths) <= 1:

        region_class = "LENGTH_CONSISTENT"

    else:

        compatible_pairs = 0
        total_pairs = 0

        for i in range(len(lengths)):
            for j in range(i + 1, len(lengths)):

                total_pairs += 1

                if likely_same_region(
                    lengths[i],
                    lengths[j]
                ):
                    compatible_pairs += 1

        if total_pairs == 0:

            region_class = "UNKNOWN"

        elif compatible_pairs == total_pairs:

            region_class = "LIKELY_SAME_REGION_LENGTH_VARIANTS"

        elif compatible_pairs > 0:

            region_class = "MIXED_REGION_PATTERN"

        else:

            region_class = "LIKELY_DIFFERENT_REGION_OR_SEQUENCE_SOURCE"


    region_rows.append(
        {
            "KJU_ID": kju,
            "Unique_lengths": ";".join(
                str(x) for x in lengths
            ),
            "Number_unique_lengths": len(lengths),
            "Region_consistency_class": region_class,
        }
    )


region_df = pd.DataFrame(region_rows)


# ============================================================
# EXACT DUPLICATE GROUPS
# ============================================================

duplicate_rows = []

if not unique_candidates.empty:

    for sha, group in (
        unique_candidates
        .groupby("Candidate_SHA256")
    ):

        if not sha:
            continue

        if len(group) <= 1:
            continue

        duplicate_rows.append(
            {
                "Sequence_SHA256": sha,
                "KJU_count": group["KJU_ID"].nunique(),
                "KJU_IDs": ";".join(
                    sorted(
                        group["KJU_ID"]
                        .dropna()
                        .unique()
                    )
                ),
                "Accession_count": group[
                    "Candidate_accession"
                ].replace("", pd.NA).nunique(
                    dropna=True
                ),
                "Accessions": ";".join(
                    sorted(
                        set(
                            x
                            for x in group[
                                "Candidate_accession"
                            ]
                            if x
                        )
                    )
                ),
                "Lengths": ";".join(
                    sorted(
                        set(
                            str(x)
                            for x in group[
                                "Candidate_length"
                            ]
                            if pd.notna(x)
                        )
                    )
                ),
            }
        )


duplicate_df = pd.DataFrame(duplicate_rows)


# ============================================================
# ACCESSION CONFLICT AUDIT
# ============================================================

accession_rows = []

if not unique_candidates.empty:

    temp = unique_candidates[
        unique_candidates["Candidate_accession"] != ""
    ].copy()

    for accession, group in (
        temp.groupby("Candidate_accession")
    ):

        shas = sorted(
            set(
                x
                for x in group[
                    "Candidate_SHA256"
                ]
                if x
            )
        )

        if len(shas) <= 1:
            continue

        accession_rows.append(
            {
                "Accession": accession,
                "Unique_sequence_count": len(shas),
                "KJU_count": group["KJU_ID"].nunique(),
                "KJU_IDs": ";".join(
                    sorted(
                        group["KJU_ID"]
                        .unique()
                    )
                ),
                "Lengths": ";".join(
                    sorted(
                        set(
                            str(x)
                            for x in group[
                                "Candidate_length"
                            ]
                            if pd.notna(x)
                        )
                    )
                ),
                "SHA256s": ";".join(shas),
            }
        )


accession_df = pd.DataFrame(accession_rows)


# ============================================================
# MANUAL REVIEW
# ============================================================

review_df = kju_df[
    kju_df["Resolution_status"].isin(
        [
            "MULTIPLE_EXACT_MATCH_CANDIDATES",
            "ACCESSION_MATCH_BUT_NO_EXACT_LENGTH",
            "LENGTH_MATCH_BUT_ACCESSION_CONFLICT",
            "REVIEW_REQUIRED",
        ]
    )
].copy()


# ============================================================
# WRITE OUTPUTS
# ============================================================

candidate_df.to_csv(
    CANDIDATES,
    index=False
)

kju_df.to_csv(
    KJU_AUDIT,
    index=False
)

review_df.to_csv(
    REVIEW,
    index=False
)

region_df.to_csv(
    REGION_AUDIT,
    index=False
)

duplicate_df.to_csv(
    DUPLICATES,
    index=False
)

accession_df.to_csv(
    ACCESSION_CONFLICTS,
    index=False
)


# ============================================================
# SUMMARY COUNTS
# ============================================================

total_kju = len(master)

resolved_kju = sum(
    bool(master_lookup[k]["_taxonomy"])
    for k in master_kjus
)

unresolved_kju = total_kju - resolved_kju

with_candidates = int(
    (kju_df["Candidate_count"] > 0).sum()
)

without_candidates = int(
    (kju_df["Candidate_count"] == 0).sum()
)

single_candidate = int(
    (
        kju_df["Unique_sequence_count"] == 1
    ).sum()
)

multiple_candidates = int(
    (
        kju_df["Unique_sequence_count"] > 1
    ).sum()
)

auto_supported = int(
    (
        kju_df["Resolution_status"]
        == "AUTO_SUPPORTED_EXACT_MATCH"
    ).sum()
)

multiple_exact = int(
    (
        kju_df["Resolution_status"]
        == "MULTIPLE_EXACT_MATCH_CANDIDATES"
    ).sum()
)

accession_only = int(
    (
        kju_df["Resolution_status"]
        == "ACCESSION_MATCH_BUT_NO_EXACT_LENGTH"
    ).sum()
)

length_only = int(
    (
        kju_df["Resolution_status"]
        == "LENGTH_MATCH_BUT_ACCESSION_CONFLICT"
    ).sum()
)

review_required = int(
    (
        kju_df["Resolution_status"]
        == "REVIEW_REQUIRED"
    ).sum()
)

no_sequence = int(
    (
        kju_df["Resolution_status"]
        == "NO_SEQUENCE"
    ).sum()
)

unresolved_status = int(
    (
        kju_df["Resolution_status"]
        == "UNRESOLVED"
    ).sum()
)


# ============================================================
# SPECIAL KJU CHECK
# ============================================================

special_ids = [
    "KJU0038",
    "KJU0039",
    "KJU0105",
    "KJU0224",
    "KJU0225",
]

special_lines = []

for kju in special_ids:

    row = kju_df[
        kju_df["KJU_ID"] == kju
    ]

    if row.empty:
        special_lines.append(
            f"{kju} | NOT FOUND"
        )
        continue

    r = row.iloc[0]

    special_lines.append(
        f"{kju} | "
        f"accepted_accession={r['Accepted_accession']} | "
        f"accepted_length={r['Accepted_length']} | "
        f"unique_sequences={r['Unique_sequence_count']} | "
        f"candidate_count={r['Candidate_count']} | "
        f"status={r['Resolution_status']}"
    )


# ============================================================
# REPORT
# ============================================================

report_lines = []

report_lines.append("=" * 110)
report_lines.append(
    "STEP 86: SEQUENCE CANDIDATE CLASSIFICATION + REGION CONSISTENCY AUDIT"
)
report_lines.append("=" * 110)

report_lines.append("")
report_lines.append("INPUTS")
report_lines.append("-" * 110)
report_lines.append(f"Frozen master : {MASTER}")
report_lines.append(f"True records  : {TRUE_RECORDS}")
report_lines.append(f"STEP84 audit  : {STEP84_AUDIT}")

report_lines.append("")
report_lines.append("MASTER")
report_lines.append("-" * 110)
report_lines.append(
    f"Total KJU records        : {total_kju}"
)
report_lines.append(
    f"Resolved records         : {resolved_kju}"
)
report_lines.append(
    f"Unresolved records       : {unresolved_kju}"
)

report_lines.append("")
report_lines.append("SEQUENCE DATA")
report_lines.append("-" * 110)
report_lines.append(
    f"Valid nucleotide rows    : {len(valid)}"
)
report_lines.append(
    f"KJUs with candidates     : {with_candidates}"
)
report_lines.append(
    f"KJUs without candidates  : {without_candidates}"
)

report_lines.append("")
report_lines.append("UNIQUE CANDIDATES")
report_lines.append("-" * 110)
report_lines.append(
    f"Exactly one unique seq   : {single_candidate}"
)
report_lines.append(
    f"Multiple unique seq      : {multiple_candidates}"
)

report_lines.append("")
report_lines.append("CANONICAL EVIDENCE CLASSIFICATION")
report_lines.append("-" * 110)
report_lines.append(
    f"Auto-supported exact     : {auto_supported}"
)
report_lines.append(
    f"Multiple exact candidates: {multiple_exact}"
)
report_lines.append(
    f"Accession match only     : {accession_only}"
)
report_lines.append(
    f"Length match only        : {length_only}"
)
report_lines.append(
    f"Review required          : {review_required}"
)
report_lines.append(
    f"No sequence              : {no_sequence}"
)
report_lines.append(
    f"Unresolved               : {unresolved_status}"
)

report_lines.append("")
report_lines.append("DUPLICATE / CONFLICT AUDIT")
report_lines.append("-" * 110)
report_lines.append(
    f"Duplicate SHA groups     : {len(duplicate_df)}"
)
report_lines.append(
    f"Accession conflict groups: {len(accession_df)}"
)

report_lines.append("")
report_lines.append("REGION CONSISTENCY")
report_lines.append("-" * 110)

if not region_df.empty:

    for cls, count in (
        region_df[
            "Region_consistency_class"
        ]
        .value_counts()
        .items()
    ):
        report_lines.append(
            f"{cls:<50}: {count}"
        )

else:

    report_lines.append(
        "No region information available."
    )


report_lines.append("")
report_lines.append("SPECIAL KJU CHECK")
report_lines.append("-" * 110)

report_lines.extend(
    special_lines
)


report_lines.append("")
report_lines.append("QC")
report_lines.append("-" * 110)

# Unresolved KJUs must never have candidates
unresolved_with_candidates = 0

for kju in master_kjus:

    m = master_lookup[kju]

    if not m["_taxonomy"]:

        row = kju_df[
            kju_df["KJU_ID"] == kju
        ]

        if not row.empty:

            if int(
                row.iloc[0]["Candidate_count"]
            ) > 0:

                unresolved_with_candidates += 1


if unresolved_with_candidates == 0:

    report_lines.append(
        "PASS — no unresolved KJU has a sequence candidate"
    )

else:

    report_lines.append(
        "FAIL — unresolved KJUs have sequence candidates: "
        + str(unresolved_with_candidates)
    )


# Check KJU count
if len(kju_df) == total_kju:

    report_lines.append(
        "PASS — every frozen KJU is represented in audit"
    )

else:

    report_lines.append(
        "FAIL — KJU audit does not contain all master records"
    )


# Check no empty KJU
if all(
    bool(x)
    for x in kju_df["KJU_ID"]
):

    report_lines.append(
        "PASS — no empty KJU IDs"
    )

else:

    report_lines.append(
        "FAIL — empty KJU ID detected"
    )


report_lines.append("")
report_lines.append("OUTPUT FILES")
report_lines.append("-" * 110)
report_lines.append(
    f"KJU audit       : {KJU_AUDIT}"
)
report_lines.append(
    f"Candidates      : {CANDIDATES}"
)
report_lines.append(
    f"Manual review   : {REVIEW}"
)
report_lines.append(
    f"Region audit    : {REGION_AUDIT}"
)
report_lines.append(
    f"Duplicates      : {DUPLICATES}"
)
report_lines.append(
    f"Accession conf. : {ACCESSION_CONFLICTS}"
)
report_lines.append(
    f"QC report       : {REPORT}"
)

report_lines.append("")
report_lines.append("=" * 110)
report_lines.append(
    "FROZEN V1 MASTER WAS NOT MODIFIED"
)
report_lines.append("=" * 110)


REPORT.write_text(
    "\n".join(report_lines),
    encoding="utf-8"
)


# ============================================================
# FINAL TERMINAL OUTPUT
# ============================================================

print()
print("==============================================================================================================")
print("STEP 86 COMPLETE")
print("==============================================================================================================")

print()
print("MASTER")
print("-" * 110)
print(f"Total KJU records        : {total_kju}")
print(f"Resolved records         : {resolved_kju}")
print(f"Unresolved records       : {unresolved_kju}")

print()
print("SEQUENCE DATA")
print("-" * 110)
print(f"Valid nucleotide rows    : {len(valid)}")
print(f"KJUs with candidates     : {with_candidates}")
print(f"KJUs without candidates  : {without_candidates}")

print()
print("UNIQUE CANDIDATES")
print("-" * 110)
print(f"Exactly one unique seq   : {single_candidate}")
print(f"Multiple unique seq      : {multiple_candidates}")

print()
print("CANONICAL EVIDENCE")
print("-" * 110)
print(f"Auto-supported exact     : {auto_supported}")
print(f"Multiple exact candidates: {multiple_exact}")
print(f"Accession match only     : {accession_only}")
print(f"Length match only        : {length_only}")
print(f"Review required          : {review_required}")
print(f"No sequence              : {no_sequence}")
print(f"Unresolved               : {unresolved_status}")

print()
print("DUPLICATE / CONFLICT")
print("-" * 110)
print(f"Duplicate SHA groups     : {len(duplicate_df)}")
print(f"Accession conflict groups: {len(accession_df)}")

print()
print("SPECIAL KJU CHECK")
print("-" * 110)

for line in special_lines:
    print(line)

print()
print("OUTPUT DIRECTORY")
print("-" * 110)
print(OUT)

print()
print("==============================================================================================================")
print("FROZEN V1 MASTER WAS NOT MODIFIED")
print("==============================================================================================================")

