#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 80
# CANONICAL SEQUENCE SELECTION + PROVENANCE AUDIT
#
# INPUT:
#   STEP 76 frozen master
#   STEP 77 raw sequence records
#   STEP 77 sequence inventory
#   STEP 79 discrepancy outputs
#
# PURPOSE:
#   1. Collapse exact duplicate sequence records
#   2. Determine unique sequences per KJU
#   3. Determine unique sequences per accession
#   4. Identify KJUs having multiple distinct sequences
#   5. Select a canonical sequence ONLY when selection is
#      supported by deterministic evidence
#   6. Preserve all alternative sequences in provenance
#   7. Never modify the frozen V1 master
#
# IMPORTANT:
#   This script DOES NOT fabricate taxonomy.
#   This script DOES NOT overwrite the frozen master.
#   Multiple-sequence KJUs are flagged for review rather than
#   arbitrarily choosing a sequence.
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

SEQ_ROOT = ROOT / "STEP77_SEQUENCE_QC"

RAW = (
    SEQ_ROOT
    / "STEP77_RAW_SEQUENCE_RECORDS.csv"
)

INVENTORY = (
    SEQ_ROOT
    / "STEP77_SEQUENCE_INVENTORY.csv"
)

OUT = (
    ROOT
    / "STEP80_CANONICAL_SEQUENCE_AUDIT"
)

OUT.mkdir(parents=True, exist_ok=True)


CANONICAL = (
    OUT
    / "STEP80_CANONICAL_SEQUENCE_SELECTION.csv"
)

PROVENANCE = (
    OUT
    / "STEP80_SEQUENCE_PROVENANCE.csv"
)

KJU_AUDIT = (
    OUT
    / "STEP80_KJU_SEQUENCE_AUDIT.csv"
)

ACCESSION_AUDIT = (
    OUT
    / "STEP80_ACCESSION_SEQUENCE_AUDIT.csv"
)

CONFLICTS = (
    OUT
    / "STEP80_SEQUENCE_CONFLICTS.csv"
)

DUPLICATES = (
    OUT
    / "STEP80_EXACT_DUPLICATES.csv"
)

REPORT = (
    OUT
    / "STEP80_CANONICAL_SEQUENCE_AUDIT_REPORT.txt"
)


# ============================================================
# HELPERS
# ============================================================

def clean(value):
    """
    Convert NaN / None / blank-like values to "".
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


def clean_upper(value):
    return clean(value).upper()


def normalize_accession(value):
    """
    Normalize GenBank accession.

    Examples:
        PZ666011.1 -> PZ666011
        PZ666011   -> PZ666011
    """
    value = clean(value)

    if not value:
        return ""

    value = value.strip()

    # Remove whitespace
    value = re.sub(r"\s+", "", value)

    # Remove version suffix
    value = re.sub(r"\.\d+$", "", value)

    return value.upper()


def normalize_sequence(value):
    """
    Normalize nucleotide sequence for hashing/comparison.
    """
    value = clean(value)

    if not value:
        return ""

    # Keep only nucleotide/IUPAC symbols
    value = re.sub(
        r"[^ACGTRYSWKMBDHVN]",
        "",
        value.upper(),
    )

    return value


def sequence_sha256(sequence):
    sequence = normalize_sequence(sequence)

    if not sequence:
        return ""

    return hashlib.sha256(
        sequence.encode("utf-8")
    ).hexdigest()


def to_int(value):
    """
    Safe numeric conversion.

    This is the key fix for the previous:
        TypeError: '>' not supported between instances of 'str' and 'int'
    """
    value = clean(value)

    if not value:
        return None

    try:
        return int(float(value))
    except Exception:
        return None


def first_existing(df, candidates):
    """
    Return first matching column from candidate list.
    """
    for col in candidates:
        if col in df.columns:
            return col
    return None


def unique_preserve(values):
    seen = set()
    result = []

    for value in values:
        value = clean(value)

        if value and value not in seen:
            seen.add(value)
            result.append(value)

    return result


# ============================================================
# HEADER
# ============================================================

print("=" * 110)
print("STEP 80: CANONICAL SEQUENCE SELECTION + PROVENANCE AUDIT")
print("=" * 110)


# ============================================================
# INPUT VALIDATION
# ============================================================

for path in [MASTER, RAW, INVENTORY]:

    if not path.exists():
        raise FileNotFoundError(
            f"Required input file not found:\n{path}"
        )


print()
print("INPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"Raw sequences : {RAW}")
print(f"Inventory     : {INVENTORY}")


# ============================================================
# LOAD INPUTS
# ============================================================

master = pd.read_csv(
    MASTER,
    dtype=str,
    keep_default_na=False,
)

raw = pd.read_csv(
    RAW,
    dtype=str,
    keep_default_na=False,
)

inventory = pd.read_csv(
    INVENTORY,
    dtype=str,
    keep_default_na=False,
)


print()
print("INPUT SIZES")
print("-" * 110)
print(f"Master rows        : {len(master)}")
print(f"Raw sequence rows  : {len(raw)}")
print(f"Inventory rows     : {len(inventory)}")


# ============================================================
# DETECT MASTER COLUMNS
# ============================================================

master_kju_col = first_existing(
    master,
    [
        "KJU_ID",
        "kju_id",
        "KJU",
    ],
)

master_acc_col = first_existing(
    master,
    [
        "Accepted_accession",
        "accepted_accession",
        "canonical_accession",
    ],
)

master_tax_col = first_existing(
    master,
    [
        "Accepted_taxonomy",
        "accepted_taxonomy",
        "canonical_organism",
    ],
)

master_len_col = first_existing(
    master,
    [
        "Accepted_length",
        "accepted_length",
        "canonical_length",
    ],
)

master_identity_col = first_existing(
    master,
    [
        "Biological_identity",
        "biological_identity",
    ],
)


if master_kju_col is None:
    raise ValueError(
        "Could not detect KJU column in frozen master."
    )

if master_acc_col is None:
    raise ValueError(
        "Could not detect accepted accession column."
    )

if master_tax_col is None:
    raise ValueError(
        "Could not detect accepted taxonomy column."
    )


# ============================================================
# DETECT RAW COLUMNS
# ============================================================

raw_kju_col = first_existing(
    raw,
    [
        "KJU_ID",
        "kju_id",
        "KJU",
    ],
)

raw_acc_col = first_existing(
    raw,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "accession",
        "Accession",
    ],
)

raw_len_col = first_existing(
    raw,
    [
        "Sequence_length",
        "sequence_length",
        "length",
    ],
)

raw_seq_col = first_existing(
    raw,
    [
        "Sequence_file",
        "Sequence",
        "sequence",
        "sequence_data",
    ],
)

raw_sha_col = first_existing(
    raw,
    [
        "Sequence_SHA256",
        "sequence_sha256",
    ],
)

raw_record_col = first_existing(
    raw,
    [
        "Record_ID",
        "record_id",
        "record",
    ],
)

raw_file_col = first_existing(
    raw,
    [
        "Sequence_file_name",
        "sequence_file_name",
        "File_name",
    ],
)


if raw_kju_col is None:
    raise ValueError(
        "Could not detect KJU column in raw sequence file."
    )

if raw_seq_col is None:
    raise ValueError(
        "Could not detect sequence column in raw sequence file."
    )


print()
print("DETECTED COLUMNS")
print("-" * 110)
print(f"Master KJU       : {master_kju_col}")
print(f"Master accession : {master_acc_col}")
print(f"Master taxonomy  : {master_tax_col}")
print(f"Master length    : {master_len_col}")
print(f"Master identity  : {master_identity_col}")
print(f"Raw KJU          : {raw_kju_col}")
print(f"Raw accession    : {raw_acc_col}")
print(f"Raw length       : {raw_len_col}")
print(f"Raw sequence     : {raw_seq_col}")
print(f"Raw SHA256       : {raw_sha_col}")


# ============================================================
# NORMALIZE MASTER
# ============================================================

master["_KJU"] = master[master_kju_col].map(clean)

master["_accepted_accession"] = (
    master[master_acc_col]
    .map(normalize_accession)
)

master["_accepted_taxonomy"] = (
    master[master_tax_col]
    .map(clean)
)

if master_len_col:
    master["_accepted_length"] = (
        master[master_len_col]
        .map(to_int)
    )
else:
    master["_accepted_length"] = None

if master_identity_col:
    master["_biological_identity"] = (
        master[master_identity_col]
        .map(clean)
    )
else:
    master["_biological_identity"] = ""


master_lookup = {}

for _, row in master.iterrows():

    kju = clean(row["_KJU"])

    if not kju:
        continue

    master_lookup[kju] = {
        "accepted_accession":
            clean(row["_accepted_accession"]),

        "accepted_taxonomy":
            clean(row["_accepted_taxonomy"]),

        "accepted_length":
            row["_accepted_length"],

        "biological_identity":
            clean(row["_biological_identity"]),
    }


# ============================================================
# PREPARE RAW SEQUENCES
# ============================================================

raw_rows = []

for idx, row in raw.iterrows():

    kju = clean(row[raw_kju_col])

    if not kju:
        continue

    sequence = normalize_sequence(
        row[raw_seq_col]
    )

    if not sequence:
        continue

    accession = ""

    if raw_acc_col:
        accession = normalize_accession(
            row[raw_acc_col]
        )

    length = len(sequence)

    supplied_length = None

    if raw_len_col:
        supplied_length = to_int(
            row[raw_len_col]
        )

    sha = ""

    if raw_sha_col:
        sha = clean(
            row[raw_sha_col]
        )

    if not sha:
        sha = sequence_sha256(sequence)

    record_id = ""

    if raw_record_col:
        record_id = clean(
            row[raw_record_col]
        )

    file_name = ""

    if raw_file_col:
        file_name = clean(
            row[raw_file_col]
        )

    raw_rows.append(
        {
            "raw_row_index": idx,
            "KJU_ID": kju,
            "Accession": accession,
            "Sequence": sequence,
            "Sequence_length": length,
            "Supplied_sequence_length": supplied_length,
            "Sequence_SHA256": sha,
            "Record_ID": record_id,
            "Sequence_file_name": file_name,
        }
    )


seq_df = pd.DataFrame(raw_rows)


if seq_df.empty:

    raise ValueError(
        "No sequence-bearing records were detected."
    )


print()
print("SEQUENCE DATA")
print("-" * 110)
print(
    f"Sequence-bearing rows : {len(seq_df)}"
)
print(
    f"KJUs represented      : "
    f"{seq_df['KJU_ID'].nunique()}"
)


# ============================================================
# EXACT DUPLICATE COLLAPSE
# ============================================================

dedup_cols = [
    "KJU_ID",
    "Accession",
    "Sequence_SHA256",
]

dedup = (
    seq_df
    .drop_duplicates(
        subset=dedup_cols,
        keep="first",
    )
    .copy()
)

seq_df["Exact_duplicate_group"] = (
    seq_df
    .groupby(dedup_cols, dropna=False)
    .ngroup()
)


print()
print("DUPLICATE COLLAPSE")
print("-" * 110)
print(
    f"Raw sequence rows       : {len(seq_df)}"
)
print(
    f"Unique KJU/accession/SHA: {len(dedup)}"
)
print(
    f"Exact duplicate rows    : "
    f"{len(seq_df) - len(dedup)}"
)


# ============================================================
# BUILD KJU-LEVEL AUDIT
# ============================================================

kju_audit_rows = []

for kju in sorted(master_lookup):

    m = master_lookup[kju]

    rows = dedup[
        dedup["KJU_ID"] == kju
    ].copy()

    unique_sha = sorted(
        set(
            x
            for x in rows["Sequence_SHA256"].map(clean)
            if x
        )
    )

    unique_accessions = sorted(
        set(
            x
            for x in rows["Accession"].map(clean)
            if x
        )
    )

    unique_lengths = sorted(
        set(
            int(x)
            for x in rows["Sequence_length"].tolist()
            if x is not None
        )
    )

    n_unique_sequences = int(
        len(unique_sha)
    )

    n_unique_accessions = int(
        len(unique_accessions)
    )

    n_unique_lengths = int(
        len(unique_lengths)
    )

    accepted_accession = (
        m["accepted_accession"]
    )

    accepted_length = (
        m["accepted_length"]
    )

    accepted_taxonomy = (
        m["accepted_taxonomy"]
    )

    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------

    has_sequence = (
        n_unique_sequences > 0
    )

    unresolved = (
        accepted_taxonomy == ""
    )

    # --------------------------------------------------------
    # Accepted accession match
    # --------------------------------------------------------

    accession_match = False

    if accepted_accession and unique_accessions:

        accession_match = (
            accepted_accession
            in unique_accessions
        )

    # --------------------------------------------------------
    # Length consistency
    # --------------------------------------------------------

    length_match = False

    if (
        accepted_length is not None
        and n_unique_lengths > 0
    ):

        length_match = (
            accepted_length
            in unique_lengths
        )

    # --------------------------------------------------------
    # Sequence multiplicity
    # --------------------------------------------------------

    if n_unique_sequences == 0:

        sequence_status = (
            "NO_SEQUENCE"
        )

    elif n_unique_sequences == 1:

        sequence_status = (
            "SINGLE_UNIQUE_SEQUENCE"
        )

    else:

        sequence_status = (
            "MULTIPLE_UNIQUE_SEQUENCES_REVIEW"
        )

    # --------------------------------------------------------
    # Canonical candidate selection
    #
    # IMPORTANT:
    # Never arbitrarily choose among multiple sequences.
    # --------------------------------------------------------

    canonical_sha = ""
    canonical_accession = ""
    canonical_length = None
    canonical_sequence = ""
    selection_status = ""

    if unresolved:

        selection_status = (
            "UNRESOLVED_NO_CANONICAL_SEQUENCE"
        )

    elif n_unique_sequences == 0:

        selection_status = (
            "NO_SEQUENCE_AVAILABLE"
        )

    elif n_unique_sequences == 1:

        candidate = rows.iloc[0]

        canonical_sha = (
            clean(candidate["Sequence_SHA256"])
        )

        canonical_accession = (
            clean(candidate["Accession"])
        )

        canonical_length = int(
            candidate["Sequence_length"]
        )

        canonical_sequence = (
            clean(candidate["Sequence"])
        )

        selection_status = (
            "CANONICAL_UNIQUE_SEQUENCE"
        )

    else:

        # Multiple distinct sequences exist.
        # Only select automatically if exactly one sequence
        # carries the accepted accession.

        accepted_candidates = rows[
            rows["Accession"]
            == accepted_accession
        ].copy()

        accepted_shas = sorted(
            set(
                accepted_candidates[
                    "Sequence_SHA256"
                ]
                .map(clean)
            )
        )

        if (
            accepted_accession
            and len(accepted_shas) == 1
        ):

            candidate_sha = (
                accepted_shas[0]
            )

            candidates = rows[
                rows["Sequence_SHA256"]
                == candidate_sha
            ]

            candidate = (
                candidates
                .sort_values(
                    [
                        "Sequence_length",
                        "raw_row_index",
                    ]
                )
                .iloc[0]
            )

            canonical_sha = (
                clean(candidate["Sequence_SHA256"])
            )

            canonical_accession = (
                clean(candidate["Accession"])
            )

            canonical_length = int(
                candidate["Sequence_length"]
            )

            canonical_sequence = (
                clean(candidate["Sequence"])
            )

            selection_status = (
                "CANONICAL_BY_ACCEPTED_ACCESSION"
            )

        else:

            selection_status = (
                "MULTIPLE_SEQUENCES_MANUAL_REVIEW"
            )

    # --------------------------------------------------------
    # Append
    # --------------------------------------------------------

    kju_audit_rows.append(
        {
            "KJU_ID":
                kju,

            "Accepted_taxonomy":
                accepted_taxonomy,

            "Accepted_accession":
                accepted_accession,

            "Accepted_length":
                accepted_length,

            "Biological_identity":
                m["biological_identity"],

            "Resolved":
                "NO"
                if unresolved
                else "YES",

            "Sequence_present":
                "YES"
                if has_sequence
                else "NO",

            "Unique_sequences":
                n_unique_sequences,

            "Unique_accessions":
                n_unique_accessions,

            "Unique_lengths":
                n_unique_lengths,

            "Accession_match":
                "YES"
                if accession_match
                else "NO",

            "Length_match":
                "YES"
                if length_match
                else "NO",

            "Sequence_status":
                sequence_status,

            "Canonical_sequence_SHA256":
                canonical_sha,

            "Canonical_accession":
                canonical_accession,

            "Canonical_length":
                canonical_length,

            "Canonical_sequence":
                canonical_sequence,

            "Selection_status":
                selection_status,
        }
    )


kju_audit = pd.DataFrame(
    kju_audit_rows
)


# ============================================================
# IMPORTANT TYPE FIX
#
# Force Unique_sequences to numeric BEFORE comparisons.
#
# This directly prevents:
#
# TypeError:
# '>' not supported between instances of 'str' and 'int'
# ============================================================

kju_audit["Unique_sequences"] = pd.to_numeric(
    kju_audit["Unique_sequences"],
    errors="coerce",
).fillna(0).astype(int)

kju_audit["Unique_accessions"] = pd.to_numeric(
    kju_audit["Unique_accessions"],
    errors="coerce",
).fillna(0).astype(int)

kju_audit["Unique_lengths"] = pd.to_numeric(
    kju_audit["Unique_lengths"],
    errors="coerce",
).fillna(0).astype(int)


# ============================================================
# CANONICAL TABLE
# ============================================================

canonical = kju_audit[
    [
        "KJU_ID",
        "Accepted_taxonomy",
        "Accepted_accession",
        "Accepted_length",
        "Biological_identity",
        "Resolved",
        "Sequence_present",
        "Unique_sequences",
        "Canonical_sequence_SHA256",
        "Canonical_accession",
        "Canonical_length",
        "Canonical_sequence",
        "Selection_status",
    ]
].copy()


# ============================================================
# PROVENANCE TABLE
# ============================================================

provenance_rows = []

for kju in sorted(master_lookup):

    rows = seq_df[
        seq_df["KJU_ID"] == kju
    ].copy()

    for _, row in rows.iterrows():

        provenance_rows.append(
            {
                "KJU_ID":
                    kju,

                "Accepted_accession":
                    master_lookup[kju][
                        "accepted_accession"
                    ],

                "Accepted_taxonomy":
                    master_lookup[kju][
                        "accepted_taxonomy"
                    ],

                "Raw_accession":
                    clean(row["Accession"]),

                "Sequence_length":
                    int(row["Sequence_length"]),

                "Sequence_SHA256":
                    clean(row["Sequence_SHA256"]),

                "Record_ID":
                    clean(row["Record_ID"]),

                "Sequence_file_name":
                    clean(row["Sequence_file_name"]),

                "Is_canonical":
                    (
                        "YES"
                        if (
                            clean(row["Sequence_SHA256"])
                            == clean(
                                kju_audit.loc[
                                    kju_audit["KJU_ID"]
                                    == kju,
                                    "Canonical_sequence_SHA256",
                                ].iloc[0]
                            )
                            and clean(
                                kju_audit.loc[
                                    kju_audit["KJU_ID"]
                                    == kju,
                                    "Canonical_sequence_SHA256",
                                ].iloc[0]
                            )
                        )
                        else "NO"
                    ),
            }
        )


provenance = pd.DataFrame(
    provenance_rows
)


# ============================================================
# ACCESSION-LEVEL AUDIT
# ============================================================

accession_rows = []

for accession in sorted(
    set(
        x
        for x in dedup["Accession"].map(clean)
        if x
    )
):

    rows = dedup[
        dedup["Accession"]
        == accession
    ].copy()

    kju_ids = sorted(
        set(
            rows["KJU_ID"].map(clean)
        )
    )

    shas = sorted(
        set(
            rows["Sequence_SHA256"].map(clean)
        )
    )

    lengths = sorted(
        set(
            int(x)
            for x in rows["Sequence_length"]
            if x is not None
        )
    )

    accession_rows.append(
        {
            "Accession":
                accession,

            "KJU_count":
                int(len(kju_ids)),

            "KJU_IDs":
                ";".join(kju_ids),

            "Unique_sequences":
                int(len(shas)),

            "Sequence_SHA256s":
                ";".join(shas),

            "Unique_lengths":
                int(len(lengths)),

            "Lengths":
                ";".join(
                    str(x)
                    for x in lengths
                ),

            "Conflict":
                (
                    "YES"
                    if len(shas) > 1
                    else "NO"
                ),
        }
    )


accession_audit = pd.DataFrame(
    accession_rows
)


# ============================================================
# CONFLICT TABLE
# ============================================================

conflict_rows = []

# ------------------------------------------------------------
# Same accession / different sequences
# ------------------------------------------------------------

for _, row in accession_audit.iterrows():

    if int(row["Unique_sequences"]) > 1:

        conflict_rows.append(
            {
                "Conflict_type":
                    "SAME_ACCESSION_DIFFERENT_SEQUENCE",

                "Accession":
                    row["Accession"],

                "KJU_IDs":
                    row["KJU_IDs"],

                "Unique_sequences":
                    int(row["Unique_sequences"]),

                "Sequence_SHA256s":
                    row["Sequence_SHA256s"],

                "Lengths":
                    row["Lengths"],
            }
        )


# ------------------------------------------------------------
# Different accession / same sequence
# ------------------------------------------------------------

sha_groups = (
    dedup
    .groupby(
        "Sequence_SHA256",
        dropna=False,
    )
)

for sha, rows in sha_groups:

    sha = clean(sha)

    if not sha:
        continue

    accessions = sorted(
        set(
            rows["Accession"].map(clean)
        )
    )

    accessions = [
        x for x in accessions if x
    ]

    if len(accessions) > 1:

        conflict_rows.append(
            {
                "Conflict_type":
                    "DIFFERENT_ACCESSION_SAME_SEQUENCE",

                "Accession":
                    ";".join(accessions),

                "KJU_IDs":
                    ";".join(
                        sorted(
                            set(
                                rows[
                                    "KJU_ID"
                                ].map(clean)
                            )
                        )
                    ),

                "Unique_sequences":
                    1,

                "Sequence_SHA256s":
                    sha,

                "Lengths":
                    ";".join(
                        sorted(
                            set(
                                str(int(x))
                                for x
                                in rows[
                                    "Sequence_length"
                                ]
                                if x is not None
                            )
                        )
                    ),
            }
        )


conflicts = pd.DataFrame(
    conflict_rows
)


# ============================================================
# EXACT DUPLICATE TABLE
# ============================================================

duplicate_groups = (
    seq_df
    .groupby(
        [
            "KJU_ID",
            "Accession",
            "Sequence_SHA256",
        ],
        dropna=False,
    )
    .size()
    .reset_index(
        name="Raw_record_count"
    )
)

duplicate_groups = duplicate_groups[
    duplicate_groups["Raw_record_count"] > 1
].copy()


# ============================================================
# SPECIAL KJU0038 / KJU0039 CHECK
# ============================================================

print()
print("=" * 110)
print("SPECIAL KJU0038 / KJU0039 CHECK")
print("-" * 110)

for kju in ["KJU0038", "KJU0039"]:

    if kju not in master_lookup:
        print(f"{kju} | NOT PRESENT")
        continue

    m = master_lookup[kju]

    rows = dedup[
        dedup["KJU_ID"] == kju
    ]

    shas = sorted(
        set(
            rows["Sequence_SHA256"].map(clean)
        )
    )

    accessions = sorted(
        set(
            rows["Accession"].map(clean)
        )
    )

    print(
        f"{kju} | "
        f"accepted_accession="
        f"{m['accepted_accession']} | "
        f"candidate_accessions="
        f"{';'.join(accessions)} | "
        f"unique_sequences="
        f"{len(shas)}"
    )


# ============================================================
# UNRESOLVED QC
# ============================================================

unresolved_kjus = {
    kju
    for kju, value in master_lookup.items()
    if not value["accepted_taxonomy"]
}

unresolved_with_sequence = sorted(
    set(
        seq_df[
            seq_df["KJU_ID"]
            .isin(unresolved_kjus)
        ]["KJU_ID"]
    )
)


# ============================================================
# SUMMARY COUNTS
# ============================================================

total_kju = int(len(master))

resolved_count = int(
    sum(
        1
        for x in master_lookup.values()
        if x["accepted_taxonomy"]
    )
)

unresolved_count = int(
    total_kju - resolved_count
)

kju_with_sequence = int(
    sum(
        1
        for x in kju_audit["Sequence_present"]
        if x == "YES"
    )
)

kju_without_sequence = int(
    total_kju - kju_with_sequence
)

exactly_one = int(
    (
        kju_audit["Unique_sequences"]
        == 1
    ).sum()
)

multiple_sequences = int(
    (
        kju_audit["Unique_sequences"]
        > 1
    ).sum()
)

same_accession_diff_seq = int(
    sum(
        conflicts["Conflict_type"]
        == "SAME_ACCESSION_DIFFERENT_SEQUENCE"
    )
) if not conflicts.empty else 0

diff_accession_same_seq = int(
    sum(
        conflicts["Conflict_type"]
        == "DIFFERENT_ACCESSION_SAME_SEQUENCE"
    )
) if not conflicts.empty else 0


# ============================================================
# SAVE OUTPUTS
# ============================================================

canonical.to_csv(
    CANONICAL,
    index=False,
)

provenance.to_csv(
    PROVENANCE,
    index=False,
)

kju_audit.to_csv(
    KJU_AUDIT,
    index=False,
)

accession_audit.to_csv(
    ACCESSION_AUDIT,
    index=False,
)

conflicts.to_csv(
    CONFLICTS,
    index=False,
)

duplicate_groups.to_csv(
    DUPLICATES,
    index=False,
)


# ============================================================
# REPORT
# ============================================================

report = []

report.append("=" * 110)
report.append(
    "STEP 80: CANONICAL SEQUENCE SELECTION + PROVENANCE AUDIT"
)
report.append("=" * 110)

report.append("")
report.append("INPUTS")
report.append("-" * 110)
report.append(f"Frozen master : {MASTER}")
report.append(f"Raw sequences : {RAW}")
report.append(f"Inventory     : {INVENTORY}")

report.append("")
report.append("MASTER")
report.append("-" * 110)
report.append(
    f"Total KJU records          : {total_kju}"
)
report.append(
    f"Resolved records           : {resolved_count}"
)
report.append(
    f"Unresolved records         : {unresolved_count}"
)

report.append("")
report.append("SEQUENCE COVERAGE")
report.append("-" * 110)
report.append(
    f"KJUs with sequence         : {kju_with_sequence}"
)
report.append(
    f"KJUs without sequence      : {kju_without_sequence}"
)
report.append(
    f"Exactly one unique seq     : {exactly_one}"
)
report.append(
    f"Multiple unique sequences  : {multiple_sequences}"
)

report.append("")
report.append("ACCESSION / SEQUENCE CONFLICTS")
report.append("-" * 110)
report.append(
    f"Same accession/diff seq    : "
    f"{same_accession_diff_seq}"
)
report.append(
    f"Different accession/same seq : "
    f"{diff_accession_same_seq}"
)

report.append("")
report.append("UNRESOLVED RECORD QC")
report.append("-" * 110)

if unresolved_with_sequence:

    report.append(
        "FAIL — unresolved KJUs with sequences:"
    )

    for kju in unresolved_with_sequence:
        report.append(
            f"  {kju}"
        )

else:

    report.append(
        "PASS — no unresolved KJU has a sequence"
    )


report.append("")
report.append("CANONICAL SELECTION POLICY")
report.append("-" * 110)
report.append(
    "Unique sequence per KJU:"
)
report.append(
    "  -> canonical sequence accepted automatically."
)
report.append(
    "Multiple sequences per KJU:"
)
report.append(
    "  -> accepted-accession sequence selected only "
    "when exactly one unique sequence matches the "
    "accepted accession."
)
report.append(
    "Multiple sequences without deterministic accession "
    "support:"
)
report.append(
    "  -> MANUAL_REVIEW; no arbitrary sequence selected."
)

report.append("")
report.append("SPECIAL KJU0038 / KJU0039")
report.append("-" * 110)

for kju in ["KJU0038", "KJU0039"]:

    if kju in master_lookup:

        m = master_lookup[kju]

        rows = dedup[
            dedup["KJU_ID"] == kju
        ]

        shas = sorted(
            set(
                rows["Sequence_SHA256"].map(clean)
            )
        )

        report.append(
            f"{kju} | "
            f"accepted_accession="
            f"{m['accepted_accession']} | "
            f"unique_sequences="
            f"{len(shas)}"
        )

report.append("")
report.append("OUTPUT FILES")
report.append("-" * 110)
report.append(
    f"Canonical       : {CANONICAL}"
)
report.append(
    f"Provenance      : {PROVENANCE}"
)
report.append(
    f"KJU audit       : {KJU_AUDIT}"
)
report.append(
    f"Accession audit : {ACCESSION_AUDIT}"
)
report.append(
    f"Conflicts       : {CONFLICTS}"
)
report.append(
    f"Duplicates      : {DUPLICATES}"
)
report.append(
    f"QC report       : {REPORT}"
)

report.append("")
report.append("=" * 110)
report.append(
    "FROZEN V1 MASTER WAS NOT MODIFIED"
)
report.append("=" * 110)


REPORT.write_text(
    "\n".join(report),
    encoding="utf-8",
)


# ============================================================
# FINAL TERMINAL OUTPUT
# ============================================================

print()
print("=" * 110)
print("STEP 80 COMPLETE")
print("=" * 110)

print()
print("MASTER")
print("-" * 110)
print(
    f"Total KJU records          : {total_kju}"
)
print(
    f"Resolved records           : {resolved_count}"
)
print(
    f"Unresolved records         : {unresolved_count}"
)

print()
print("SEQUENCE COVERAGE")
print("-" * 110)
print(
    f"KJUs with sequence         : {kju_with_sequence}"
)
print(
    f"KJUs without sequence      : {kju_without_sequence}"
)
print(
    f"Exactly one unique seq     : {exactly_one}"
)
print(
    f"Multiple unique sequences  : {multiple_sequences}"
)

print()
print("ACCESSION / SEQUENCE CONFLICTS")
print("-" * 110)
print(
    f"Same accession/diff seq    : "
    f"{same_accession_diff_seq}"
)
print(
    f"Different accession/same seq : "
    f"{diff_accession_same_seq}"
)

print()
print("UNRESOLVED RECORD QC")
print("-" * 110)

if unresolved_with_sequence:

    print(
        "FAIL — unresolved KJUs with sequences:"
    )

    for kju in unresolved_with_sequence:
        print(f"  {kju}")

else:

    print(
        "PASS — no unresolved KJU has a sequence"
    )


print()
print("CANONICAL SELECTION")
print("-" * 110)

automatic = int(
    (
        kju_audit["Selection_status"]
        .isin(
            [
                "CANONICAL_UNIQUE_SEQUENCE",
                "CANONICAL_BY_ACCEPTED_ACCESSION",
            ]
        )
    ).sum()
)

manual_review = int(
    (
        kju_audit["Selection_status"]
        == "MULTIPLE_SEQUENCES_MANUAL_REVIEW"
    ).sum()
)

no_sequence = int(
    (
        kju_audit["Selection_status"]
        == "NO_SEQUENCE_AVAILABLE"
    ).sum()
)

unresolved_selection = int(
    (
        kju_audit["Selection_status"]
        == "UNRESOLVED_NO_CANONICAL_SEQUENCE"
    ).sum()
)

print(
    f"Automatically selected    : {automatic}"
)

print(
    f"Manual review required    : {manual_review}"
)

print(
    f"No sequence available     : {no_sequence}"
)

print(
    f"Unresolved/no canonical   : "
    f"{unresolved_selection}"
)


print()
print("OUTPUT DIRECTORY")
print("-" * 110)
print(OUT)

print()
print("=" * 110)
print(
    "FROZEN V1 MASTER WAS NOT MODIFIED"
)
print("=" * 110)
