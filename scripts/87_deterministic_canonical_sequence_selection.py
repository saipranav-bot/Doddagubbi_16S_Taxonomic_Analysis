#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
import math


# ============================================================
# STEP 87
# DETERMINISTIC CANONICAL SEQUENCE SELECTION
#
# INPUT:
#   STEP 76 frozen master
#   STEP 83 true nucleotide sequence records
#   STEP 86 sequence candidate classification
#
# PURPOSE:
#
#   Select a canonical REAL nucleotide sequence for each KJU
#   ONLY when deterministic evidence supports the selection.
#
# Evidence hierarchy:
#
#   PRIORITY 1
#       Accepted accession + accepted length + exact sequence
#
#   PRIORITY 2
#       Accepted accession + exact sequence
#       when accepted length is unavailable
#
#   PRIORITY 3
#       Exactly one unique sequence candidate
#
#   Otherwise:
#       REVIEW_REQUIRED
#
# IMPORTANT:
#   - Frozen V1 master is NEVER modified.
#   - No taxonomy is changed.
#   - No accession is changed.
#   - No sequence is fabricated.
#   - Blank/malformed KJU IDs are preserved in a separate QC file.
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

STEP86 = (
    ROOT
    / "STEP86_SEQUENCE_CANDIDATE_CLASSIFICATION"
)

OUT = (
    ROOT
    / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
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

STEP86_AUDIT = (
    STEP86
    / "STEP86_KJU_CANDIDATE_AUDIT.csv"
)


# ============================================================
# OUTPUT FILES
# ============================================================

KJU_AUDIT = (
    OUT
    / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"
)

CANONICAL = (
    OUT
    / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
)

REVIEW = (
    OUT
    / "STEP87_MANUAL_REVIEW_REQUIRED.csv"
)

CANDIDATES = (
    OUT
    / "STEP87_ALL_SEQUENCE_CANDIDATES.csv"
)

INVALID_KJU = (
    OUT
    / "STEP87_INVALID_OR_BLANK_KJU_RECORDS.csv"
)

MASTER_CHECK = (
    OUT
    / "STEP87_FROZEN_MASTER_INTEGRITY.txt"
)

REPORT = (
    OUT
    / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE_REPORT.txt"
)


# ============================================================
# HELPERS
# ============================================================

def clean_string(value):
    """
    Convert pandas missing values to empty string and strip whitespace.
    """
    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value.lower() in {
        "",
        "nan",
        "none",
        "null",
        "na",
        "n/a"
    }:
        return ""

    return value


def normalize_kju(value):
    """
    Normalize KJU IDs conservatively.

    Expected form:
        KJU0001
        KJU0038
        KJU0204

    No new KJU IDs are fabricated.
    """

    value = clean_string(value).upper()

    if not value:
        return ""

    # Remove surrounding whitespace only.
    value = re.sub(r"\s+", "", value)

    # Accept canonical KJU format only.
    if re.fullmatch(r"KJU\d{4}", value):
        return value

    # Allow KJU followed by 1-3 digits and zero-pad.
    m = re.fullmatch(r"KJU(\d{1,3})", value)

    if m:
        number = int(m.group(1))
        return f"KJU{number:04d}"

    return value


def normalize_accession(value):
    """
    Normalize accession for comparison.

    Example:
        PZ666011
        PZ666011.1

    are treated as the same accession root for conservative matching.
    """

    value = clean_string(value)

    if not value:
        return ""

    value = value.upper()

    # Remove whitespace.
    value = re.sub(r"\s+", "", value)

    # Remove version suffix only for comparison.
    value = re.sub(r"\.\d+$", "", value)

    return value


def numeric_length(value):
    """
    Safely convert accepted/sequence length to integer.

    NaN, blank, invalid values -> None.
    """

    if pd.isna(value):
        return None

    value = clean_string(value)

    if not value:
        return None

    try:
        number = float(value)

        if math.isnan(number):
            return None

        if number <= 0:
            return None

        return int(number)

    except Exception:
        return None


def normalize_sequence(seq):
    """
    Normalize a real nucleotide sequence.
    """

    seq = clean_string(seq)

    if not seq:
        return ""

    seq = seq.upper()

    # Remove FASTA whitespace/newlines if present.
    seq = re.sub(r"\s+", "", seq)

    return seq


def sequence_sha256(seq):
    """
    Calculate SHA256 from normalized nucleotide sequence.
    """

    seq = normalize_sequence(seq)

    if not seq:
        return ""

    return hashlib.sha256(
        seq.encode("utf-8")
    ).hexdigest()


# ============================================================
# START
# ============================================================

print("=" * 110)
print("STEP 87: DETERMINISTIC CANONICAL SEQUENCE SELECTION")
print("=" * 110)


# ============================================================
# VERIFY INPUTS
# ============================================================

for path in [MASTER, TRUE_RECORDS, STEP86_AUDIT]:

    if not path.exists():
        raise FileNotFoundError(
            f"Required input file does not exist:\n{path}"
        )


print("\nINPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"True records  : {TRUE_RECORDS}")
print(f"STEP86 audit  : {STEP86_AUDIT}")


# ============================================================
# FROZEN MASTER HASH BEFORE
# ============================================================

master_hash_before = hashlib.sha256(
    MASTER.read_bytes()
).hexdigest()


# ============================================================
# LOAD MASTER
# ============================================================

master = pd.read_csv(
    MASTER,
    dtype=str,
    low_memory=False
)

true = pd.read_csv(
    TRUE_RECORDS,
    dtype=str,
    low_memory=False
)

step86 = pd.read_csv(
    STEP86_AUDIT,
    dtype=str,
    low_memory=False
)


print("\nINPUT SIZES")
print("-" * 110)
print(f"Master rows       : {len(master)}")
print(f"True sequence rows: {len(true)}")
print(f"STEP86 audit rows : {len(step86)}")


# ============================================================
# REQUIRED COLUMNS
# ============================================================

required_master = [
    "KJU_ID",
    "Accepted_accession",
    "Accepted_taxonomy",
    "Accepted_length",
    "Biological_identity",
]

required_true = [
    "KJU_ID",
    "Sequence_accession_detected",
    "Sequence_length",
    "Sequence_SHA256",
    "Sequence_type",
    "Record_ID",
    "Sequence",
]


for col in required_master:

    if col not in master.columns:
        raise ValueError(
            f"Frozen master missing required column: {col}"
        )


for col in required_true:

    if col not in true.columns:
        raise ValueError(
            f"STEP83 true sequence records missing required column: {col}"
        )


# ============================================================
# NORMALIZE MASTER KJU IDs
# ============================================================

master["_KJU_NORMALIZED"] = (
    master["KJU_ID"]
    .map(normalize_kju)
)


master_ids = set(
    x
    for x in master["_KJU_NORMALIZED"]
    if x
)


if len(master_ids) != len(master):

    duplicate_master_ids = (
        master[
            master["_KJU_NORMALIZED"].duplicated(keep=False)
        ]
        ["_KJU_NORMALIZED"]
        .tolist()
    )

    raise ValueError(
        "Frozen master contains duplicate/invalid KJU IDs. "
        f"IDs involved: {sorted(set(duplicate_master_ids))}"
    )


print("\nMASTER VALIDATION")
print("-" * 110)
print(f"Master KJU IDs           : {len(master_ids)}")
print(f"Expected frozen records  : 204")


# ============================================================
# NORMALIZE TRUE SEQUENCE KJU IDs
# ============================================================

true["_KJU_NORMALIZED"] = (
    true["KJU_ID"]
    .map(normalize_kju)
)


# ============================================================
# IDENTIFY INVALID / BLANK KJU IDs
# ============================================================

invalid_mask = (
    true["_KJU_NORMALIZED"] == ""
)

invalid_rows = true[
    invalid_mask
].copy()


# ============================================================
# IDENTIFY KJU IDS ABSENT FROM MASTER
# ============================================================

unknown_mask = (
    (~invalid_mask)
    &
    (~true["_KJU_NORMALIZED"].isin(master_ids))
)

unknown_rows = true[
    unknown_mask
].copy()


invalid_or_unknown = pd.concat(
    [
        invalid_rows,
        unknown_rows
    ],
    ignore_index=True
)


# ============================================================
# SAVE INVALID / UNKNOWN RECORDS
# ============================================================

if len(invalid_or_unknown) > 0:

    invalid_or_unknown.to_csv(
        INVALID_KJU,
        index=False
    )

else:

    pd.DataFrame(
        columns=true.columns
    ).to_csv(
        INVALID_KJU,
        index=False
    )


print("\nKJU ID VALIDATION")
print("-" * 110)
print(f"Blank/invalid KJU rows : {len(invalid_rows)}")
print(f"Unknown KJU rows       : {len(unknown_rows)}")


# ============================================================
# CRITICAL SAFETY CHECK
# ============================================================
#
# We DO NOT silently discard unknown KJUs.
#
# However, blank rows are allowed to remain outside the
# canonical selection dataset because they cannot be safely
# assigned to any KJU.
#
# Genuine nonblank unknown KJUs remain a HARD ERROR.
# ============================================================

if len(unknown_rows) > 0:

    unknown_ids = sorted(
        set(
            unknown_rows[
                "_KJU_NORMALIZED"
            ]
            .tolist()
        )
    )

    raise ValueError(
        "TRUE sequence records contain nonblank KJU IDs absent "
        "from frozen master:\n"
        + "\n".join(unknown_ids)
        + "\n\nSee:\n"
        + str(INVALID_KJU)
    )


# ============================================================
# EXCLUDE ONLY BLANK/INVALID KJU ROWS
# ============================================================

true_valid = true[
    ~invalid_mask
].copy()


print(
    f"Sequence rows eligible for canonical selection : "
    f"{len(true_valid)}"
)


# ============================================================
# NORMALIZE SEQUENCE FIELDS
# ============================================================

true_valid["_ACCESSION_NORMALIZED"] = (
    true_valid[
        "Sequence_accession_detected"
    ]
    .map(normalize_accession)
)


true_valid["_SEQUENCE_LENGTH"] = (
    true_valid[
        "Sequence_length"
    ]
    .map(numeric_length)
)


true_valid["_SEQUENCE_NORMALIZED"] = (
    true_valid[
        "Sequence"
    ]
    .map(normalize_sequence)
)


# ============================================================
# VALIDATE ACTUAL SEQUENCE
# ============================================================

true_valid["_CALCULATED_SHA256"] = (
    true_valid[
        "_SEQUENCE_NORMALIZED"
    ]
    .map(sequence_sha256)
)


true_valid["_SEQUENCE_VALID"] = (
    true_valid[
        "_SEQUENCE_NORMALIZED"
    ]
    .str.fullmatch(
        r"[ACGTNRYWSMKBDHVU]+",
        na=False
    )
)


invalid_sequence_rows = true_valid[
    (~true_valid["_SEQUENCE_VALID"])
    |
    (true_valid["_SEQUENCE_NORMALIZED"] == "")
].copy()


if len(invalid_sequence_rows) > 0:

    raise ValueError(
        "STEP83 contains invalid nucleotide sequences. "
        f"Invalid rows: {len(invalid_sequence_rows)}"
    )


# ============================================================
# USE CALCULATED SHA AS TRUE SEQUENCE IDENTITY
# ============================================================

true_valid["_TRUE_SHA256"] = (
    true_valid["_CALCULATED_SHA256"]
)


# ============================================================
# LOAD MASTER LOOKUP
# ============================================================

master_lookup = {}

for _, row in master.iterrows():

    kju = row["_KJU_NORMALIZED"]

    master_lookup[kju] = {
        "Accepted_accession":
            clean_string(
                row["Accepted_accession"]
            ),

        "Accepted_accession_normalized":
            normalize_accession(
                row["Accepted_accession"]
            ),

        "Accepted_taxonomy":
            clean_string(
                row["Accepted_taxonomy"]
            ),

        "Accepted_length":
            numeric_length(
                row["Accepted_length"]
            ),

        "Biological_identity":
            clean_string(
                row["Biological_identity"]
            ),
    }


# ============================================================
# BUILD UNIQUE CANDIDATES
# ============================================================
#
# One candidate = one KJU + one true sequence SHA.
#
# Multiple source records representing exactly the same
# nucleotide sequence collapse into one candidate.
# ============================================================

candidate_groups = (
    true_valid
    .groupby(
        [
            "_KJU_NORMALIZED",
            "_TRUE_SHA256"
        ],
        dropna=False
    )
)


candidate_rows = []


for (kju, sha), group in candidate_groups:

    info = master_lookup[kju]

    accessions = sorted(
        set(
            x
            for x in group[
                "_ACCESSION_NORMALIZED"
            ]
            if x
        )
    )

    raw_accessions = sorted(
        set(
            clean_string(x)
            for x in group[
                "Sequence_accession_detected"
            ]
            if clean_string(x)
        )
    )

    lengths = sorted(
        set(
            x
            for x in group[
                "_SEQUENCE_LENGTH"
            ]
            if x is not None
        )
    )

    sequence_types = sorted(
        set(
            clean_string(x)
            for x in group[
                "Sequence_type"
            ]
            if clean_string(x)
        )
    )

    record_ids = sorted(
        set(
            clean_string(x)
            for x in group[
                "Record_ID"
            ]
            if clean_string(x)
        )
    )

    seq = group.iloc[0][
        "_SEQUENCE_NORMALIZED"
    ]

    candidate_rows.append(
        {
            "KJU_ID":
                kju,

            "Sequence_SHA256":
                sha,

            "Sequence_length":
                len(seq),

            "Candidate_accessions":
                ";".join(accessions),

            "Raw_candidate_accessions":
                ";".join(raw_accessions),

            "Candidate_lengths":
                ";".join(
                    str(x)
                    for x in lengths
                ),

            "Sequence_types":
                ";".join(sequence_types),

            "Source_record_count":
                len(group),

            "Source_record_IDs":
                ";".join(record_ids),

            "Sequence":
                seq,

            "Accepted_accession":
                info[
                    "Accepted_accession"
                ],

            "Accepted_accession_normalized":
                info[
                    "Accepted_accession_normalized"
                ],

            "Accepted_length":
                info[
                    "Accepted_length"
                ],

            "Accepted_taxonomy":
                info[
                    "Accepted_taxonomy"
                ],

            "Biological_identity":
                info[
                    "Biological_identity"
                ],
        }
    )


candidates = pd.DataFrame(
    candidate_rows
)


print("\nCANDIDATE COLLAPSE")
print("-" * 110)
print(
    f"Raw valid sequence rows : "
    f"{len(true_valid)}"
)
print(
    f"Unique KJU/SHA candidates: "
    f"{len(candidates)}"
)


# ============================================================
# SAVE ALL CANDIDATES
# ============================================================

candidates.to_csv(
    CANDIDATES,
    index=False
)


# ============================================================
# DETERMINE UNIQUE CANDIDATE COUNT PER KJU
# ============================================================

candidate_count = (
    candidates
    .groupby("KJU_ID")
    .size()
    .to_dict()
)


# ============================================================
# DETERMINISTIC SELECTION
# ============================================================

selection_rows = []


for _, master_row in master.iterrows():

    kju = master_row[
        "_KJU_NORMALIZED"
    ]

    info = master_lookup[kju]

    accepted_acc = (
        info[
            "Accepted_accession_normalized"
        ]
    )

    accepted_len = (
        info[
            "Accepted_length"
        ]
    )

    kju_candidates = candidates[
        candidates[
            "KJU_ID"
        ] == kju
    ].copy()

    n_candidates = len(
        kju_candidates
    )

    selected = None
    status = ""
    evidence = ""

    # --------------------------------------------------------
    # NO SEQUENCE
    # --------------------------------------------------------

    if n_candidates == 0:

        if kju in set(
            true_valid[
                "_KJU_NORMALIZED"
            ]
        ):

            status = "NO_CANONICAL_CANDIDATE"

        else:

            if clean_string(
                master_row[
                    "Accepted_taxonomy"
                ]
            ):

                # Resolved master record without sequence.
                status = "NO_SEQUENCE"

            else:

                status = "UNRESOLVED"

        evidence = "NONE"


    # --------------------------------------------------------
    # PRIORITY 1
    #
    # Accepted accession + accepted length + exact sequence
    # --------------------------------------------------------

    else:

        priority1 = kju_candidates[
            (
                kju_candidates[
                    "Accepted_accession_normalized"
                ].eq(
                    kju_candidates[
                        "Candidate_accessions"
                    ]
                )
            )
            &
            (
                kju_candidates[
                    "Candidate_lengths"
                ].apply(
                    lambda x:
                    (
                        accepted_len is not None
                        and
                        str(accepted_len)
                        in str(x).split(";")
                    )
                )
            )
        ]

        # The candidate accession field may contain multiple
        # accessions. Therefore also perform explicit token test.

        if accepted_acc:

            priority1 = kju_candidates[
                kju_candidates[
                    "Candidate_accessions"
                ].apply(
                    lambda x:
                    accepted_acc in
                    str(x).split(";")
                )
                &
                kju_candidates[
                    "Candidate_lengths"
                ].apply(
                    lambda x:
                    (
                        accepted_len is not None
                        and
                        str(accepted_len)
                        in str(x).split(";")
                    )
                )
            ]

        else:

            priority1 = kju_candidates.iloc[0:0]


        if len(priority1) == 1:

            selected = priority1.iloc[0]

            status = (
                "AUTO_SELECTED_PRIORITY_1"
            )

            evidence = (
                "ACCEPTED_ACCESSION_PLUS_ACCEPTED_LENGTH"
            )


        # ----------------------------------------------------
        # PRIORITY 2
        #
        # Accepted accession + exact sequence
        # when accepted length unavailable
        # ----------------------------------------------------

        elif (
            accepted_acc
            and
            accepted_len is None
        ):

            priority2 = kju_candidates[
                kju_candidates[
                    "Candidate_accessions"
                ].apply(
                    lambda x:
                    accepted_acc in
                    str(x).split(";")
                )
            ]

            if len(priority2) == 1:

                selected = priority2.iloc[0]

                status = (
                    "AUTO_SELECTED_PRIORITY_2"
                )

                evidence = (
                    "ACCEPTED_ACCESSION_ONLY_LENGTH_UNAVAILABLE"
                )


        # ----------------------------------------------------
        # PRIORITY 3
        #
        # Exactly one unique sequence candidate
        # ----------------------------------------------------

        if (
            selected is None
            and
            n_candidates == 1
        ):

            selected = (
                kju_candidates.iloc[0]
            )

            status = (
                "AUTO_SELECTED_PRIORITY_3"
            )

            evidence = (
                "EXACTLY_ONE_UNIQUE_SEQUENCE"
            )


        # ----------------------------------------------------
        # OTHERWISE REVIEW
        # ----------------------------------------------------

        if selected is None:

            status = (
                "REVIEW_REQUIRED"
            )

            evidence = (
                "MULTIPLE_SEQUENCE_CANDIDATES_OR_CONFLICT"
            )


    # ========================================================
    # BUILD AUDIT ROW
    # ========================================================

    row = {
        "KJU_ID":
            kju,

        "Accepted_accession":
            info[
                "Accepted_accession"
            ],

        "Accepted_length":
            info[
                "Accepted_length"
            ],

        "Accepted_taxonomy":
            info[
                "Accepted_taxonomy"
            ],

        "Biological_identity":
            info[
                "Biological_identity"
            ],

        "Unique_candidate_count":
            n_candidates,

        "Selection_status":
            status,

        "Selection_evidence":
            evidence,

        "Canonical_sequence_SHA256":
            "",

        "Canonical_sequence_length":
            "",

        "Canonical_accession":
            "",

        "Canonical_sequence":
            "",
    }


    if selected is not None:

        row[
            "Canonical_sequence_SHA256"
        ] = selected[
            "Sequence_SHA256"
        ]

        row[
            "Canonical_sequence_length"
        ] = selected[
            "Sequence_length"
        ]

        row[
            "Canonical_accession"
        ] = selected[
            "Candidate_accessions"
        ]

        row[
            "Canonical_sequence"
        ] = selected[
            "Sequence"
        ]


    selection_rows.append(row)


selection = pd.DataFrame(
    selection_rows
)


# ============================================================
# SAVE AUDIT
# ============================================================

selection.to_csv(
    KJU_AUDIT,
    index=False
)


# ============================================================
# CANONICAL SEQUENCES
# ============================================================

canonical = selection[
    selection[
        "Selection_status"
    ].str.startswith(
        "AUTO_SELECTED_",
        na=False
    )
].copy()


canonical.to_csv(
    CANONICAL,
    index=False
)


# ============================================================
# MANUAL REVIEW
# ============================================================

review = selection[
    selection[
        "Selection_status"
    ] == "REVIEW_REQUIRED"
].copy()


review.to_csv(
    REVIEW,
    index=False
)


# ============================================================
# SUMMARY COUNTS
# ============================================================

priority1_count = (
    selection[
        "Selection_status"
    ]
    ==
    "AUTO_SELECTED_PRIORITY_1"
).sum()


priority2_count = (
    selection[
        "Selection_status"
    ]
    ==
    "AUTO_SELECTED_PRIORITY_2"
).sum()


priority3_count = (
    selection[
        "Selection_status"
    ]
    ==
    "AUTO_SELECTED_PRIORITY_3"
).sum()


review_count = (
    selection[
        "Selection_status"
    ]
    ==
    "REVIEW_REQUIRED"
).sum()


no_sequence_count = (
    selection[
        "Selection_status"
    ].isin(
        [
            "NO_SEQUENCE",
            "NO_CANONICAL_CANDIDATE"
        ]
    )
).sum()


unresolved_count = (
    selection[
        "Selection_status"
    ]
    ==
    "UNRESOLVED"
).sum()


auto_count = (
    selection[
        "Selection_status"
    ].str.startswith(
        "AUTO_SELECTED_",
        na=False
    )
).sum()


# ============================================================
# SPECIAL KJU CHECKS
# ============================================================

print("\n" + "=" * 110)
print("SPECIAL KJU CHECK")
print("=" * 110)


for kju in [
    "KJU0038",
    "KJU0039",
    "KJU0105",
    "KJU0224",
    "KJU0225",
]:

    x = selection[
        selection[
            "KJU_ID"
        ] == kju
    ]

    if len(x) == 1:

        x = x.iloc[0]

        print(
            f"{kju} | "
            f"accepted_accession="
            f"{x['Accepted_accession']} | "
            f"accepted_length="
            f"{x['Accepted_length']} | "
            f"unique_sequences="
            f"{x['Unique_candidate_count']} | "
            f"status="
            f"{x['Selection_status']} | "
            f"evidence="
            f"{x['Selection_evidence']}"
        )


# ============================================================
# QC CHECKS
# ============================================================

print("\n" + "=" * 110)
print("STEP 87 QC")
print("=" * 110)


# Must have exactly 204 audit rows.
if len(selection) != 204:

    raise ValueError(
        f"Expected 204 selection rows, got {len(selection)}"
    )


# Must have exactly 204 unique KJUs.
if selection["KJU_ID"].nunique() != 204:

    raise ValueError(
        "Canonical selection audit does not contain "
        "204 unique KJU IDs."
    )


# No unresolved KJU should have a sequence.
master_unresolved = set(
    master[
        master[
            "Accepted_taxonomy"
        ]
        .map(clean_string)
        == ""
    ][
        "_KJU_NORMALIZED"
    ]
)

unresolved_with_sequence = selection[
    (
        selection[
            "KJU_ID"
        ].isin(master_unresolved)
    )
    &
    (
        selection[
            "Unique_candidate_count"
        ]
        > 0
    )
]


if len(unresolved_with_sequence) > 0:

    raise ValueError(
        "Unresolved KJUs have sequence candidates:\n"
        + "\n".join(
            unresolved_with_sequence[
                "KJU_ID"
            ].tolist()
        )
    )


# Every automatically selected sequence must have SHA.
auto_without_sha = canonical[
    canonical[
        "Canonical_sequence_SHA256"
    ].map(clean_string)
    == ""
]


if len(auto_without_sha) > 0:

    raise ValueError(
        "Automatically selected records have blank SHA256."
    )


# Every automatically selected sequence must have sequence.
auto_without_sequence = canonical[
    canonical[
        "Canonical_sequence"
    ].map(clean_string)
    == ""
]


if len(auto_without_sequence) > 0:

    raise ValueError(
        "Automatically selected records have blank sequence."
    )


# ============================================================
# FROZEN MASTER HASH AFTER
# ============================================================

master_hash_after = hashlib.sha256(
    MASTER.read_bytes()
).hexdigest()


if master_hash_before != master_hash_after:

    raise RuntimeError(
        "CRITICAL: Frozen V1 master changed during STEP 87."
    )


# ============================================================
# MASTER INTEGRITY REPORT
# ============================================================

with open(
    MASTER_CHECK,
    "w"
) as f:

    f.write(
        "STEP 87 FROZEN MASTER INTEGRITY\n"
    )

    f.write(
        "=" * 90
        + "\n"
    )

    f.write(
        f"SHA256 BEFORE : "
        f"{master_hash_before}\n"
    )

    f.write(
        f"SHA256 AFTER  : "
        f"{master_hash_after}\n"
    )

    f.write(
        f"STATUS        : "
        f"{'PASS' if master_hash_before == master_hash_after else 'FAIL'}\n"
    )


# ============================================================
# FINAL REPORT
# ============================================================

with open(
    REPORT,
    "w"
) as f:

    f.write(
        "STEP 87: DETERMINISTIC CANONICAL SEQUENCE SELECTION\n"
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
        f"True records  : {TRUE_RECORDS}\n"
    )

    f.write(
        f"STEP86 audit  : {STEP86_AUDIT}\n\n"
    )

    f.write(
        "MASTER\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"Total KJU records : {len(master)}\n"
    )

    f.write(
        f"Unique KJU IDs    : {len(master_ids)}\n"
    )

    f.write(
        f"Resolved records  : "
        f"{len(master) - len(master_unresolved)}\n"
    )

    f.write(
        f"Unresolved records: "
        f"{len(master_unresolved)}\n\n"
    )

    f.write(
        "KJU INPUT VALIDATION\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"Blank/invalid KJU rows : "
        f"{len(invalid_rows)}\n"
    )

    f.write(
        f"Unknown KJU rows       : "
        f"{len(unknown_rows)}\n"
    )

    f.write(
        f"Eligible sequence rows: "
        f"{len(true_valid)}\n\n"
    )

    f.write(
        "CANONICAL SELECTION\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"Priority 1 selected : "
        f"{priority1_count}\n"
    )

    f.write(
        f"Priority 2 selected : "
        f"{priority2_count}\n"
    )

    f.write(
        f"Priority 3 selected : "
        f"{priority3_count}\n"
    )

    f.write(
        f"Total automatically selected : "
        f"{auto_count}\n"
    )

    f.write(
        f"Manual review required : "
        f"{review_count}\n"
    )

    f.write(
        f"No sequence available : "
        f"{no_sequence_count}\n"
    )

    f.write(
        f"Unresolved : "
        f"{unresolved_count}\n\n"
    )

    f.write(
        "SAFETY\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        "No taxonomy modified.\n"
    )

    f.write(
        "No accession modified.\n"
    )

    f.write(
        "No sequence fabricated.\n"
    )

    f.write(
        "Frozen V1 master was not modified.\n"
    )

    f.write(
        f"Frozen V1 SHA256: "
        f"{master_hash_after}\n\n"
    )

    f.write(
        "OUTPUTS\n"
    )

    f.write(
        "-" * 110
        + "\n"
    )

    f.write(
        f"KJU audit       : {KJU_AUDIT}\n"
    )

    f.write(
        f"Canonical        : {CANONICAL}\n"
    )

    f.write(
        f"Manual review    : {REVIEW}\n"
    )

    f.write(
        f"All candidates   : {CANDIDATES}\n"
    )

    f.write(
        f"Invalid KJUs     : {INVALID_KJU}\n"
    )

    f.write(
        f"Master integrity : {MASTER_CHECK}\n"
    )

    f.write(
        f"QC report        : {REPORT}\n"
    )


# ============================================================
# FINAL CONSOLE
# ============================================================

print("\n" + "=" * 110)
print("STEP 87 COMPLETE")
print("=" * 110)

print("\nMASTER")
print("-" * 110)
print(f"Total KJU records        : {len(master)}")
print(f"Unique KJU IDs           : {len(master_ids)}")

print("\nINPUT KJU QC")
print("-" * 110)
print(f"Blank/invalid KJU rows   : {len(invalid_rows)}")
print(f"Unknown KJU rows         : {len(unknown_rows)}")
print(f"Eligible sequence rows   : {len(true_valid)}")

print("\nCANONICAL SELECTION")
print("-" * 110)
print(f"Priority 1 selected      : {priority1_count}")
print(f"Priority 2 selected      : {priority2_count}")
print(f"Priority 3 selected      : {priority3_count}")
print(f"Automatically selected   : {auto_count}")
print(f"Manual review required   : {review_count}")
print(f"No sequence available    : {no_sequence_count}")
print(f"Unresolved               : {unresolved_count}")

print("\nQC")
print("-" * 110)
print("PASS — 204 KJU audit rows generated")
print("PASS — no unknown nonblank KJU IDs")
print("PASS — no unresolved KJU has sequence")
print("PASS — automatically selected sequences have SHA256")
print("PASS — frozen V1 master was not modified")

print("\nOUTPUT DIRECTORY")
print("-" * 110)
print(OUT)

print("\nCANONICAL FILE")
print("-" * 110)
print(CANONICAL)

print("\nMANUAL REVIEW FILE")
print("-" * 110)
print(REVIEW)

print("\nQC REPORT")
print("-" * 110)
print(REPORT)

print("\nFROZEN V1 SHA256")
print("-" * 110)
print(master_hash_after)

print("\n" + "=" * 110)
print("STEP 87 COMPLETE — FROZEN V1 MASTER UNTOUCHED")
print("=" * 110)
