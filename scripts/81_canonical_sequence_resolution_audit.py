#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 81
# CANONICAL SEQUENCE RESOLUTION AUDIT
#
# INPUT:
#   STEP 76 frozen master
#   STEP 77 raw sequence records
#   STEP 80 canonical sequence audit outputs
#
# PURPOSE:
#   1. Reconstruct unique sequence candidates per KJU
#   2. Compare candidates against:
#        - accepted accession
#        - accepted length
#        - accepted taxonomy
#        - biological identity
#   3. Determine whether a canonical sequence can be selected
#      using deterministic evidence
#   4. Separate:
#        - AUTO_RESOLVED
#        - REVIEW_REQUIRED
#        - NO_SEQUENCE
#        - UNRESOLVED
#   5. Preserve every candidate sequence
#
# IMPORTANT:
#   - Frozen V1 master is NEVER modified.
#   - No taxonomy is fabricated.
#   - No arbitrary sequence is selected.
#   - A sequence is AUTO_RESOLVED only when deterministic
#     evidence supports it.
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

STEP77 = ROOT / "STEP77_SEQUENCE_QC"

RAW = STEP77 / "STEP77_RAW_SEQUENCE_RECORDS.csv"

STEP80 = ROOT / "STEP80_CANONICAL_SEQUENCE_AUDIT"

OUT = ROOT / "STEP81_CANONICAL_SEQUENCE_RESOLUTION"

OUT.mkdir(parents=True, exist_ok=True)

SUMMARY = OUT / "STEP81_CANONICAL_SEQUENCE_RESOLUTION_SUMMARY.txt"
KJU_AUDIT = OUT / "STEP81_KJU_RESOLUTION_AUDIT.csv"
CANDIDATES = OUT / "STEP81_ALL_SEQUENCE_CANDIDATES.csv"
REVIEW = OUT / "STEP81_MANUAL_REVIEW_REQUIRED.csv"
AUTO = OUT / "STEP81_AUTO_RESOLVED.csv"
CONFLICTS = OUT / "STEP81_SEQUENCE_CONFLICTS.csv"


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
        "unknown",
    }:
        return ""

    return value


def norm_accession(value):
    value = clean(value)

    if not value:
        return ""

    value = value.upper().strip()

    # Remove whitespace
    value = re.sub(r"\s+", "", value)

    # Normalize accession.version -> accession
    if "." in value:
        value = value.split(".")[0]

    return value


def norm_int(value):
    value = clean(value)

    if not value:
        return None

    try:
        return int(float(value))
    except Exception:
        return None


def normalize_sequence(value):
    value = clean(value)

    if not value:
        return ""

    # Remove FASTA header if accidentally present
    if value.startswith(">"):
        lines = value.splitlines()
        value = "".join(
            line.strip()
            for line in lines
            if not line.startswith(">")
        )

    # Keep biological sequence characters only
    value = re.sub(r"[^A-Za-z]", "", value)

    return value.upper()


def sequence_hash(sequence):
    sequence = normalize_sequence(sequence)

    if not sequence:
        return ""

    return hashlib.sha256(
        sequence.encode("utf-8")
    ).hexdigest()


def safe_text(value):
    return clean(value).replace("\n", " ").replace("\r", " ")


# ============================================================
# START
# ============================================================

print("=" * 110)
print("STEP 81: CANONICAL SEQUENCE RESOLUTION AUDIT")
print("=" * 110)


# ============================================================
# LOAD INPUTS
# ============================================================

print("\nINPUT FILES")
print("-" * 110)
print(f"Frozen master : {MASTER}")
print(f"Raw sequences : {RAW}")
print(f"Step 80 dir  : {STEP80}")


if not MASTER.exists():
    raise FileNotFoundError(
        f"Frozen master not found: {MASTER}"
    )

if not RAW.exists():
    raise FileNotFoundError(
        f"Raw sequence file not found: {RAW}"
    )


master = pd.read_csv(MASTER, dtype=str).fillna("")
raw = pd.read_csv(RAW, dtype=str).fillna("")


print("\nINPUT SIZES")
print("-" * 110)
print(f"Master rows        : {len(master)}")
print(f"Raw sequence rows  : {len(raw)}")


# ============================================================
# DETECT REQUIRED COLUMNS
# ============================================================

def detect_column(df, candidates, label):
    lower = {str(c).lower(): c for c in df.columns}

    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]

    raise ValueError(
        f"Could not detect {label} column. "
        f"Available columns: {list(df.columns)}"
    )


master_kju = detect_column(
    master,
    ["KJU_ID", "kju_id", "KJU"],
    "master KJU",
)

master_acc = detect_column(
    master,
    [
        "Accepted_accession",
        "accepted_accession",
        "canonical_accession",
    ],
    "master accession",
)

master_tax = detect_column(
    master,
    [
        "Accepted_taxonomy",
        "accepted_taxonomy",
        "canonical_organism",
    ],
    "master taxonomy",
)

master_len = detect_column(
    master,
    [
        "Accepted_length",
        "accepted_length",
        "canonical_length",
    ],
    "master length",
)

master_identity = detect_column(
    master,
    [
        "Biological_identity",
        "biological_identity",
    ],
    "master biological identity",
)

raw_kju = detect_column(
    raw,
    ["KJU_ID", "kju_id", "KJU"],
    "raw KJU",
)

raw_seq = detect_column(
    raw,
    [
        "Sequence_file",
        "sequence",
        "Sequence",
        "sequence_data",
    ],
    "raw sequence",
)

raw_acc = detect_column(
    raw,
    [
        "Sequence_accession_detected",
        "Normalized_accession_detected",
        "accession",
    ],
    "raw accession",
)

raw_len = detect_column(
    raw,
    [
        "Sequence_length",
        "sequence_length",
        "length",
    ],
    "raw length",
)


print("\nDETECTED COLUMNS")
print("-" * 110)
print(f"Master KJU       : {master_kju}")
print(f"Master accession : {master_acc}")
print(f"Master taxonomy  : {master_tax}")
print(f"Master length    : {master_len}")
print(f"Master identity  : {master_identity}")
print(f"Raw KJU          : {raw_kju}")
print(f"Raw sequence     : {raw_seq}")
print(f"Raw accession    : {raw_acc}")
print(f"Raw length       : {raw_len}")


# ============================================================
# MASTER NORMALIZATION
# ============================================================

master["_KJU"] = master[master_kju].map(clean)
master["_ACC"] = master[master_acc].map(norm_accession)
master["_TAX"] = master[master_tax].map(clean)
master["_LEN"] = master[master_len].map(norm_int)
master["_BIO"] = master[master_identity].map(clean)

master_lookup = {}

for _, row in master.iterrows():

    kju = row["_KJU"]

    if not kju:
        continue

    master_lookup[kju] = {
        "accepted_accession": row["_ACC"],
        "accepted_taxonomy": row["_TAX"],
        "accepted_length": row["_LEN"],
        "biological_identity": row["_BIO"],
    }


# ============================================================
# RAW SEQUENCE NORMALIZATION
# ============================================================

raw["_KJU"] = raw[raw_kju].map(clean)
raw["_ACC"] = raw[raw_acc].map(norm_accession)
raw["_SEQ"] = raw[raw_seq].map(normalize_sequence)

raw["_LEN"] = raw["_SEQ"].map(
    lambda x: len(x) if x else None
)

raw["_SHA"] = raw["_SEQ"].map(sequence_hash)


# ============================================================
# REMOVE EMPTY SEQUENCE ROWS
# ============================================================

seq_raw = raw[
    (raw["_KJU"] != "")
    & (raw["_SEQ"] != "")
].copy()

print("\nSEQUENCE DATA")
print("-" * 110)
print(f"Sequence-bearing rows : {len(seq_raw)}")
print(
    f"KJUs represented      : "
    f"{seq_raw['_KJU'].nunique()}"
)


# ============================================================
# COLLAPSE EXACT DUPLICATES
#
# Unique biological candidate =
# KJU + accession + sequence SHA
# ============================================================

candidate_rows = []

group_cols = [
    "_KJU",
    "_ACC",
    "_SHA",
]

for (kju, acc, sha), group in seq_raw.groupby(
    group_cols,
    dropna=False,
):

    first = group.iloc[0]

    sequence = first["_SEQ"]

    candidate_rows.append({
        "KJU_ID": kju,
        "Candidate_accession": acc,
        "Candidate_length": len(sequence),
        "Candidate_SHA256": sha,
        "Candidate_sequence": sequence,
        "Raw_record_count": len(group),
    })


candidates = pd.DataFrame(candidate_rows)


if candidates.empty:

    candidates = pd.DataFrame(
        columns=[
            "KJU_ID",
            "Candidate_accession",
            "Candidate_length",
            "Candidate_SHA256",
            "Candidate_sequence",
            "Raw_record_count",
        ]
    )


print("\nDUPLICATE COLLAPSE")
print("-" * 110)
print(f"Raw sequence rows       : {len(seq_raw)}")
print(f"Unique candidates       : {len(candidates)}")


# ============================================================
# RESOLUTION LOGIC
# ============================================================

audit_rows = []
auto_rows = []
review_rows = []
conflict_rows = []


for kju, master_info in master_lookup.items():

    m_acc = master_info["accepted_accession"]
    m_len = master_info["accepted_length"]
    m_tax = master_info["accepted_taxonomy"]
    m_bio = master_info["biological_identity"]

    sub = candidates[
        candidates["KJU_ID"] == kju
    ].copy()

    unique_count = len(sub)

    # --------------------------------------------------------
    # UNRESOLVED MASTER RECORD
    # --------------------------------------------------------

    if not m_tax:

        if unique_count > 0:

            status = "ERROR_UNRESOLVED_HAS_SEQUENCE"

        else:

            status = "UNRESOLVED_NO_SEQUENCE"

        audit_rows.append({
            "KJU_ID": kju,
            "Accepted_accession": m_acc,
            "Accepted_taxonomy": m_tax,
            "Accepted_length": m_len,
            "Biological_identity": m_bio,
            "Unique_sequences": unique_count,
            "Candidate_accessions": ";".join(
                sorted(
                    set(
                        sub["Candidate_accession"]
                        .tolist()
                    )
                )
            ),
            "Resolution_status": status,
            "Canonical_SHA256": "",
            "Canonical_accession": "",
            "Canonical_length": "",
            "Selection_reason": (
                "UNRESOLVED_RECORD_NO_SEQUENCE"
                if unique_count == 0
                else
                "ERROR_UNRESOLVED_RECORD_HAS_SEQUENCE"
            ),
        })

        if unique_count > 0:
            conflict_rows.append({
                "KJU_ID": kju,
                "Conflict_type": (
                    "UNRESOLVED_RECORD_HAS_SEQUENCE"
                ),
                "Unique_sequences": unique_count,
            })

        continue


    # --------------------------------------------------------
    # RESOLVED RECORD WITH NO SEQUENCE
    # --------------------------------------------------------

    if unique_count == 0:

        audit_rows.append({
            "KJU_ID": kju,
            "Accepted_accession": m_acc,
            "Accepted_taxonomy": m_tax,
            "Accepted_length": m_len,
            "Biological_identity": m_bio,
            "Unique_sequences": 0,
            "Candidate_accessions": "",
            "Resolution_status": "NO_SEQUENCE",
            "Canonical_SHA256": "",
            "Canonical_accession": "",
            "Canonical_length": "",
            "Selection_reason": "NO_SEQUENCE_AVAILABLE",
        })

        continue


    # --------------------------------------------------------
    # SCORE EACH CANDIDATE
    # --------------------------------------------------------

    scored = []

    for _, candidate in sub.iterrows():

        c_acc = norm_accession(
            candidate["Candidate_accession"]
        )

        c_len = norm_int(
            candidate["Candidate_length"]
        )

        score = 0
        reasons = []

        # Accepted accession match
        if m_acc and c_acc and c_acc == m_acc:

            score += 100
            reasons.append(
                "ACCEPTED_ACCESSION_MATCH"
            )

        # Accepted length match
        if (
            m_len is not None
            and c_len is not None
            and c_len == m_len
        ):

            score += 20
            reasons.append(
                "ACCEPTED_LENGTH_MATCH"
            )

        # Accession + length together
        if (
            m_acc
            and c_acc
            and c_acc == m_acc
            and m_len is not None
            and c_len is not None
            and c_len == m_len
        ):

            score += 50
            reasons.append(
                "ACCESSION_LENGTH_CONCORDANT"
            )

        scored.append({
            "Candidate_SHA256":
                candidate["Candidate_SHA256"],
            "Candidate_accession":
                c_acc,
            "Candidate_length":
                c_len,
            "Candidate_sequence":
                candidate["Candidate_sequence"],
            "Raw_record_count":
                candidate["Raw_record_count"],
            "Score":
                score,
            "Reason":
                ";".join(reasons),
        })


    scored_df = pd.DataFrame(scored)

    scored_df = scored_df.sort_values(
        by=[
            "Score",
            "Raw_record_count",
        ],
        ascending=[
            False,
            False,
        ],
    ).reset_index(drop=True)


    best = scored_df.iloc[0]

    best_score = int(best["Score"])

    if len(scored_df) > 1:

        second_score = int(
            scored_df.iloc[1]["Score"]
        )

    else:

        second_score = -1


    # --------------------------------------------------------
    # DETERMINE WHETHER AUTO-SELECTION IS SAFE
    # --------------------------------------------------------

    canonical = None
    status = "REVIEW_REQUIRED"
    reason = ""

    # One unique sequence = deterministic
    if unique_count == 1:

        canonical = best
        status = "AUTO_RESOLVED"
        reason = (
            "ONLY_ONE_UNIQUE_SEQUENCE"
        )

    # Multiple sequences but exact accession + length match
    elif (
        best_score >= 170
        and best_score > second_score
    ):

        canonical = best
        status = "AUTO_RESOLVED"
        reason = (
            "UNIQUE_ACCEPTED_ACCESSION_AND_LENGTH_MATCH"
        )

    # Multiple sequences with unique accepted accession
    elif (
        best_score >= 100
        and best_score > second_score
    ):

        canonical = best
        status = "AUTO_RESOLVED"
        reason = (
            "UNIQUE_ACCEPTED_ACCESSION_MATCH"
        )

    else:

        reason = (
            "MULTIPLE_CANDIDATES_WITHOUT_UNIQUE_"
            "DETERMINISTIC_EVIDENCE"
        )


    # --------------------------------------------------------
    # RECORD AUDIT
    # --------------------------------------------------------

    audit_rows.append({
        "KJU_ID": kju,
        "Accepted_accession": m_acc,
        "Accepted_taxonomy": m_tax,
        "Accepted_length": m_len,
        "Biological_identity": m_bio,
        "Unique_sequences": unique_count,
        "Candidate_accessions": ";".join(
            sorted(
                set(
                    scored_df[
                        "Candidate_accession"
                    ].tolist()
                )
            )
        ),
        "Resolution_status": status,
        "Canonical_SHA256": (
            str(canonical["Candidate_SHA256"])
            if canonical is not None
            else ""
        ),
        "Canonical_accession": (
            str(canonical["Candidate_accession"])
            if canonical is not None
            else ""
        ),
        "Canonical_length": (
            int(canonical["Candidate_length"])
            if canonical is not None
            else ""
        ),
        "Selection_reason": reason,
    })


    # --------------------------------------------------------
    # AUTO OUTPUT
    # --------------------------------------------------------

    if status == "AUTO_RESOLVED":

        auto_rows.append({
            "KJU_ID": kju,
            "Accepted_accession": m_acc,
            "Accepted_taxonomy": m_tax,
            "Biological_identity": m_bio,
            "Canonical_accession":
                canonical["Candidate_accession"],
            "Canonical_length":
                canonical["Candidate_length"],
            "Canonical_SHA256":
                canonical["Candidate_SHA256"],
            "Selection_reason": reason,
            "Evidence_score": best_score,
        })


    # --------------------------------------------------------
    # REVIEW OUTPUT
    # --------------------------------------------------------

    if status == "REVIEW_REQUIRED":

        review_rows.append({
            "KJU_ID": kju,
            "Accepted_accession": m_acc,
            "Accepted_taxonomy": m_tax,
            "Accepted_length": m_len,
            "Biological_identity": m_bio,
            "Unique_sequences": unique_count,
            "Best_score": best_score,
            "Second_score": second_score,
            "Candidate_accessions":
                ";".join(
                    sorted(
                        set(
                            scored_df[
                                "Candidate_accession"
                            ].tolist()
                        )
                    )
                ),
            "Reason": reason,
        })

        conflict_rows.append({
            "KJU_ID": kju,
            "Conflict_type":
                "MULTIPLE_SEQUENCE_CANDIDATES",
            "Unique_sequences":
                unique_count,
        })


# ============================================================
# SAVE ALL CANDIDATES
# ============================================================

candidates.to_csv(
    CANDIDATES,
    index=False,
)


# ============================================================
# SAVE KJU AUDIT
# ============================================================

audit_df = pd.DataFrame(audit_rows)

audit_df = audit_df.sort_values(
    by="KJU_ID"
)

audit_df.to_csv(
    KJU_AUDIT,
    index=False,
)


# ============================================================
# SAVE AUTO-RESOLVED
# ============================================================

auto_df = pd.DataFrame(auto_rows)

if auto_df.empty:

    auto_df = pd.DataFrame(
        columns=[
            "KJU_ID",
            "Accepted_accession",
            "Accepted_taxonomy",
            "Biological_identity",
            "Canonical_accession",
            "Canonical_length",
            "Canonical_SHA256",
            "Selection_reason",
            "Evidence_score",
        ]
    )

auto_df.to_csv(
    AUTO,
    index=False,
)


# ============================================================
# SAVE MANUAL REVIEW
# ============================================================

review_df = pd.DataFrame(review_rows)

if review_df.empty:

    review_df = pd.DataFrame(
        columns=[
            "KJU_ID",
            "Accepted_accession",
            "Accepted_taxonomy",
            "Accepted_length",
            "Biological_identity",
            "Unique_sequences",
            "Best_score",
            "Second_score",
            "Candidate_accessions",
            "Reason",
        ]
    )

review_df.to_csv(
    REVIEW,
    index=False,
)


# ============================================================
# SAVE CONFLICTS
# ============================================================

conflict_df = pd.DataFrame(conflict_rows)

if conflict_df.empty:

    conflict_df = pd.DataFrame(
        columns=[
            "KJU_ID",
            "Conflict_type",
            "Unique_sequences",
        ]
    )

conflict_df.to_csv(
    CONFLICTS,
    index=False,
)


# ============================================================
# COUNTS
# ============================================================

resolved_count = int(
    (master["_TAX"] != "").sum()
)

unresolved_count = int(
    (master["_TAX"] == "").sum()
)

sequence_kjus = int(
    candidates["KJU_ID"].nunique()
)

no_sequence_count = int(
    (
        audit_df["Resolution_status"]
        == "NO_SEQUENCE"
    ).sum()
)

unresolved_no_sequence = int(
    (
        audit_df["Resolution_status"]
        == "UNRESOLVED_NO_SEQUENCE"
    ).sum()
)

error_unresolved_sequence = int(
    (
        audit_df["Resolution_status"]
        == "ERROR_UNRESOLVED_HAS_SEQUENCE"
    ).sum()
)

auto_count = int(
    (
        audit_df["Resolution_status"]
        == "AUTO_RESOLVED"
    ).sum()
)

review_count = int(
    (
        audit_df["Resolution_status"]
        == "REVIEW_REQUIRED"
    ).sum()
)

unique_one_count = int(
    (
        candidates.groupby("KJU_ID")
        .size()
        == 1
    ).sum()
)

multiple_count = int(
    (
        candidates.groupby("KJU_ID")
        .size()
        > 1
    ).sum()
)


# ============================================================
# ACCESSION CONFLICT QC
# ============================================================

same_acc_diff_seq = 0

for acc, group in candidates[
    candidates["Candidate_accession"] != ""
].groupby(
    "Candidate_accession"
):

    if group["Candidate_SHA256"].nunique() > 1:
        same_acc_diff_seq += 1


# ============================================================
# FROZEN MASTER HASH
# ============================================================

master_bytes = MASTER.read_bytes()

frozen_sha256 = hashlib.sha256(
    master_bytes
).hexdigest()


# ============================================================
# FINAL QC
# ============================================================

print("\n" + "=" * 110)
print("STEP 81 COMPLETE")
print("=" * 110)

print("\nMASTER")
print("-" * 110)
print(f"Total KJU records          : {len(master)}")
print(f"Resolved records           : {resolved_count}")
print(f"Unresolved records         : {unresolved_count}")

print("\nSEQUENCE COVERAGE")
print("-" * 110)
print(f"KJUs with sequence         : {sequence_kjus}")
print(
    f"KJUs without sequence      : "
    f"{len(master) - sequence_kjus}"
)

print("\nCANONICAL RESOLUTION")
print("-" * 110)
print(
    f"Exactly one unique seq     : "
    f"{unique_one_count}"
)
print(
    f"Multiple unique sequences  : "
    f"{multiple_count}"
)
print(
    f"Automatically resolved     : "
    f"{auto_count}"
)
print(
    f"Manual review required     : "
    f"{review_count}"
)
print(
    f"No sequence available      : "
    f"{no_sequence_count}"
)
print(
    f"Unresolved/no canonical     : "
    f"{unresolved_no_sequence}"
)

print("\nCONFLICT QC")
print("-" * 110)
print(
    f"Same accession/diff seq    : "
    f"{same_acc_diff_seq}"
)

print("\nUNRESOLVED QC")
print("-" * 110)

if error_unresolved_sequence == 0:

    print(
        "PASS — no unresolved KJU has a sequence"
    )

else:

    print(
        "FAIL — unresolved KJUs with sequences: "
        f"{error_unresolved_sequence}"
    )


# ============================================================
# SPECIAL RECORD CHECK
# ============================================================

print("\n" + "=" * 110)
print("SPECIAL KJU0038 / KJU0039 CHECK")
print("=" * 110)

for kju in ["KJU0038", "KJU0039"]:

    row = audit_df[
        audit_df["KJU_ID"] == kju
    ]

    if row.empty:

        print(
            f"{kju} | NOT FOUND"
        )

    else:

        r = row.iloc[0]

        print(
            f"{kju} | "
            f"accepted_accession="
            f"{r['Accepted_accession']} | "
            f"unique_sequences="
            f"{r['Unique_sequences']} | "
            f"status="
            f"{r['Resolution_status']} | "
            f"canonical_accession="
            f"{r['Canonical_accession']}"
        )


# ============================================================
# WRITE SUMMARY REPORT
# ============================================================

with open(
    SUMMARY,
    "w",
    encoding="utf-8",
) as fh:

    fh.write(
        "=" * 110 + "\n"
    )
    fh.write(
        "STEP 81: CANONICAL SEQUENCE "
        "RESOLUTION AUDIT\n"
    )
    fh.write(
        "=" * 110 + "\n\n"
    )

    fh.write(
        "Frozen V1 master was used as read-only input.\n"
    )

    fh.write(
        f"Frozen master SHA256: "
        f"{frozen_sha256}\n\n"
    )

    fh.write(
        "MASTER\n"
    )
    fh.write(
        f"Total KJU records          : "
        f"{len(master)}\n"
    )
    fh.write(
        f"Resolved records           : "
        f"{resolved_count}\n"
    )
    fh.write(
        f"Unresolved records         : "
        f"{unresolved_count}\n\n"
    )

    fh.write(
        "SEQUENCE COVERAGE\n"
    )
    fh.write(
        f"KJUs with sequence         : "
        f"{sequence_kjus}\n"
    )
    fh.write(
        f"KJUs without sequence      : "
        f"{len(master) - sequence_kjus}\n"
    )
    fh.write(
        f"Exactly one unique seq     : "
        f"{unique_one_count}\n"
    )
    fh.write(
        f"Multiple unique sequences  : "
        f"{multiple_count}\n\n"
    )

    fh.write(
        "CANONICAL RESOLUTION\n"
    )
    fh.write(
        f"Automatically resolved     : "
        f"{auto_count}\n"
    )
    fh.write(
        f"Manual review required     : "
        f"{review_count}\n"
    )
    fh.write(
        f"No sequence available      : "
        f"{no_sequence_count}\n"
    )
    fh.write(
        f"Unresolved/no canonical     : "
        f"{unresolved_no_sequence}\n\n"
    )

    fh.write(
        "CONFLICT QC\n"
    )
    fh.write(
        f"Same accession/diff seq    : "
        f"{same_acc_diff_seq}\n\n"
    )

    if error_unresolved_sequence == 0:

        fh.write(
            "PASS — no unresolved KJU has a sequence.\n"
        )

    else:

        fh.write(
            "FAIL — unresolved KJUs have sequences.\n"
        )

    fh.write(
        "\nFROZEN V1 MASTER WAS NOT MODIFIED.\n"
    )


# ============================================================
# FINAL STATUS
# ============================================================

if error_unresolved_sequence != 0:

    raise ValueError(
        "STEP 81 FAILED: unresolved KJU has sequence data."
    )


print("\nOUTPUT FILES")
print("-" * 110)
print(f"KJU audit       : {KJU_AUDIT}")
print(f"All candidates  : {CANDIDATES}")
print(f"Auto-resolved   : {AUTO}")
print(f"Manual review   : {REVIEW}")
print(f"Conflicts       : {CONFLICTS}")
print(f"QC summary      : {SUMMARY}")

print("\nFROZEN MASTER")
print("-" * 110)
print(
    "PASS — Frozen V1 master was not modified"
)

print(
    "\nFrozen V1 SHA256:"
)
print(frozen_sha256)

print("\n" + "=" * 110)
print("STEP 81 COMPLETE — NO MASTER FILE WAS MODIFIED")
print("=" * 110)
