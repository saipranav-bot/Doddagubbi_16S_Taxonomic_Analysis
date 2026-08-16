#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
from collections import defaultdict
import hashlib
import re


# ============================================================
# STEP 77
# SEQUENCE INVENTORY + METADATA CONSISTENCY QC
#
# INPUT:
#   FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv
#
# PURPOSE:
#   1. Inventory all sequence files in the project
#   2. Identify KJU IDs associated with sequences
#   3. Check FASTA/GenBank sequence availability
#   4. Check sequence length against accepted metadata
#   5. Check accession consistency
#   6. Detect duplicate sequences
#   7. Classify 16S / ITS / unknown
#   8. Verify unresolved records remain unresolved
#   9. Produce a standalone QC report
#
# IMPORTANT:
#   The frozen V1 master is NEVER modified.
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

FREEZE = (
    ROOT
    / "FINAL_KJU_MASTER_FREEZE_V1"
)

MASTER_FILE = (
    FREEZE
    / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
)

OUTDIR = (
    ROOT
    / "STEP77_SEQUENCE_QC"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

INVENTORY_FILE = (
    OUTDIR
    / "STEP77_SEQUENCE_INVENTORY.csv"
)

QC_FILE = (
    OUTDIR
    / "STEP77_SEQUENCE_QC_REPORT.txt"
)

DUP_FILE = (
    OUTDIR
    / "STEP77_DUPLICATE_SEQUENCES.csv"
)

MISSING_FILE = (
    OUTDIR
    / "STEP77_MISSING_SEQUENCE_RECORDS.csv"
)

MISMATCH_FILE = (
    OUTDIR
    / "STEP77_LENGTH_MISMATCHES.csv"
)


# ============================================================
# HELPERS
# ============================================================

def clean(x):
    if x is None:
        return ""

    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass

    x = str(x).strip()

    if x.lower() in {
        "",
        "nan",
        "none",
        "null",
        "na",
        "n/a"
    }:
        return ""

    return x


def normalize_accession(x):

    x = clean(x)

    if not x:
        return ""

    x = x.split(".")[0]

    return x.upper()


def normalize_sequence(seq):

    if not seq:
        return ""

    seq = re.sub(
        r"[^A-Za-z]",
        "",
        seq
    )

    return seq.upper()


def sequence_type_from_sequence(seq):

    seq = normalize_sequence(seq)

    if not seq:
        return "NO_SEQUENCE"

    dna = set(seq)

    if dna <= set("ACGTN"):
        return "DNA"

    return "NON_STANDARD"


def classify_expected_sequence_type(
    accepted_type,
    accepted_taxonomy,
    description=""
):

    text = " ".join(
        [
            clean(accepted_type),
            clean(accepted_taxonomy),
            clean(description),
        ]
    ).lower()

    if "its" in text:
        return "ITS"

    if "16s" in text:
        return "16S"

    return "UNKNOWN"


def parse_fasta(path):

    records = []

    current_id = None
    current_desc = ""
    current_seq = []

    try:

        with open(
            path,
            "r",
            errors="ignore"
        ) as fh:

            for line in fh:

                line = line.rstrip()

                if not line:
                    continue

                if line.startswith(">"):

                    if current_id is not None:

                        seq = normalize_sequence(
                            "".join(current_seq)
                        )

                        records.append(
                            (
                                current_id,
                                current_desc,
                                seq
                            )
                        )

                    header = line[1:].strip()

                    parts = header.split(
                        None,
                        1
                    )

                    current_id = parts[0]

                    current_desc = (
                        parts[1]
                        if len(parts) > 1
                        else ""
                    )

                    current_seq = []

                else:

                    current_seq.append(line)

        if current_id is not None:

            seq = normalize_sequence(
                "".join(current_seq)
            )

            records.append(
                (
                    current_id,
                    current_desc,
                    seq
                )
            )

    except Exception as e:

        print(
            f"WARNING: Could not read FASTA {path}: {e}"
        )

    return records


def parse_genbank(path):

    records = []

    try:

        from Bio import SeqIO

        for rec in SeqIO.parse(
            str(path),
            "genbank"
        ):

            records.append(
                (
                    clean(rec.id),
                    clean(rec.description),
                    normalize_sequence(
                        str(rec.seq)
                    )
                )
            )

    except Exception as e:

        print(
            f"WARNING: Could not read GenBank {path}: {e}"
        )

    return records


# ============================================================
# LOAD MASTER
# ============================================================

print("=" * 110)
print("STEP 77: SEQUENCE INVENTORY + METADATA CONSISTENCY QC")
print("=" * 110)

print()
print("INPUT")
print("-" * 100)
print(MASTER_FILE)

if not MASTER_FILE.exists():

    raise FileNotFoundError(
        f"Frozen master not found: {MASTER_FILE}"
    )


df = pd.read_csv(
    MASTER_FILE,
    dtype=str
).fillna("")


df.columns = [
    c.strip()
    for c in df.columns
]


# ============================================================
# MASTER QC
# ============================================================

assert len(df) == 204

assert (
    df["KJU_ID"].nunique()
    == 204
)


print()
print("MASTER VALIDATION")
print("-" * 100)

print(
    f"Rows                 : {len(df)}"
)

print(
    f"Unique KJU IDs       : {df['KJU_ID'].nunique()}"
)

print(
    f"Resolved             : "
    f"{(
        df['Adjudication_type']
        != 'MASTER_ONLY_UNRESOLVED'
    ).sum()}"
)

print(
    f"Unresolved           : "
    f"{(
        df['Adjudication_type']
        == 'MASTER_ONLY_UNRESOLVED'
    ).sum()}"
)


# ============================================================
# FIND ALL SEQUENCE FILES
# ============================================================

print()
print("SEARCHING FOR SEQUENCE FILES...")
print("-" * 100)


sequence_extensions = {
    ".fa",
    ".fasta",
    ".fas",
    ".fna",
    ".ffn",
    ".gb",
    ".gbk",
    ".genbank"
}


files = []

for path in BASE.rglob("*"):

    if not path.is_file():
        continue

    if (
        path.suffix.lower()
        in sequence_extensions
    ):

        files.append(path)


# Remove files generated inside STEP77
files = [
    p for p in files
    if OUTDIR not in p.parents
]


files = sorted(
    set(files)
)


print(
    f"Sequence files found : {len(files)}"
)


# ============================================================
# INVENTORY
# ============================================================

inventory = []

for path in files:

    suffix = path.suffix.lower()

    if suffix in {
        ".gb",
        ".gbk",
        ".genbank"
    }:

        records = parse_genbank(path)

        file_format = "GENBANK"

    else:

        records = parse_fasta(path)

        file_format = "FASTA"


    for record_id, description, sequence in records:

        if not sequence:
            continue


        header_text = (
            f"{record_id} "
            f"{description} "
            f"{path.name}"
        )


        # ----------------------------------------------------
        # KJU ID extraction
        # ----------------------------------------------------

        kju_match = re.search(
            r"\bKJU\d{4}\b",
            header_text,
            flags=re.IGNORECASE
        )

        kju_id = (
            kju_match.group(0).upper()
            if kju_match
            else ""
        )


        # ----------------------------------------------------
        # ACCESSION extraction
        # ----------------------------------------------------

        accession = ""

        accession_patterns = [

            r"\b[A-Z]{1,4}\d{5,8}(?:\.\d+)?\b",

            r"\b[A-Z]{2,5}_\d{5,9}(?:\.\d+)?\b",

        ]


        for pattern in accession_patterns:

            match = re.search(
                pattern,
                header_text
            )

            if match:

                candidate = match.group(0)

                if not candidate.upper().startswith(
                    "KJU"
                ):

                    accession = candidate
                    break


        normalized_acc = (
            normalize_accession(accession)
        )


        sequence_hash = hashlib.sha256(
            sequence.encode()
        ).hexdigest()


        inventory.append({

            "KJU_ID": kju_id,

            "Sequence_file": str(path),

            "Sequence_file_name": path.name,

            "File_format": file_format,

            "Record_ID": record_id,

            "Description": description,

            "Sequence_accession_detected": accession,

            "Normalized_accession_detected":
                normalized_acc,

            "Sequence_length":
                len(sequence),

            "Sequence_type_detected":
                sequence_type_from_sequence(
                    sequence
                ),

            "Sequence_SHA256":
                sequence_hash,

            "Sequence_present":
                "YES",

        })


# ============================================================
# INVENTORY DATAFRAME
# ============================================================

inv = pd.DataFrame(
    inventory
)


if len(inv) == 0:

    print()
    print(
        "WARNING: No sequence records were detected."
    )

    inv = pd.DataFrame(
        columns=[
            "KJU_ID",
            "Sequence_file",
            "Sequence_file_name",
            "File_format",
            "Record_ID",
            "Description",
            "Sequence_accession_detected",
            "Normalized_accession_detected",
            "Sequence_length",
            "Sequence_type_detected",
            "Sequence_SHA256",
            "Sequence_present",
        ]
    )


# ============================================================
# KJU SEQUENCE MAPPING
# ============================================================

kju_to_sequences = {}

for _, row in inv.iterrows():

    kju = clean(
        row["KJU_ID"]
    )

    if not kju:
        continue

    kju_to_sequences.setdefault(
        kju,
        []
    ).append(row)


# ============================================================
# BUILD MASTER-LEVEL SEQUENCE QC
# ============================================================

rows = []

for _, m in df.iterrows():

    kju = clean(
        m["KJU_ID"]
    )

    accepted_acc = normalize_accession(
        m["Accepted_accession"]
    )

    accepted_taxonomy = clean(
        m["Accepted_taxonomy"]
    )

    accepted_length = clean(
        m["Accepted_length"]
    )

    accepted_sequence_type = clean(
        m["Accepted_sequence_type"]
    )

    unresolved = (
        clean(
            m["Adjudication_type"]
        )
        == "MASTER_ONLY_UNRESOLVED"
    )


    matches = kju_to_sequences.get(
        kju,
        []
    )


    sequence_count = len(matches)


    if sequence_count == 0:

        sequence_status = "NO_SEQUENCE"

        detected_accessions = ""

        detected_lengths = ""

        detected_types = ""

        sequence_hashes = ""

    else:

        sequence_status = "SEQUENCE_FOUND"

        detected_accessions = ";".join(
            sorted(
                set(
                    clean(
                        x[
                            "Normalized_accession_detected"
                        ]
                    )
                    for x in matches
                    if clean(
                        x[
                            "Normalized_accession_detected"
                        ]
                    )
                )
            )
        )

        detected_lengths = ";".join(
            sorted(
                set(
                    str(
                        x["Sequence_length"]
                    )
                    for x in matches
                )
            )
        )

        detected_types = ";".join(
            sorted(
                set(
                    clean(
                        x[
                            "Sequence_type_detected"
                        ]
                    )
                    for x in matches
                )
            )
        )

        sequence_hashes = ";".join(
            sorted(
                set(
                    x["Sequence_SHA256"]
                    for x in matches
                )
            )
        )


    # --------------------------------------------------------
    # ACCESSION QC
    # --------------------------------------------------------

    accession_match = "NOT_CHECKED"

    if accepted_acc:

        if detected_accessions:

            if accepted_acc in (
                detected_accessions
                .split(";")
            ):

                accession_match = "MATCH"

            else:

                accession_match = "MISMATCH"

        else:

            accession_match = "NO_SEQUENCE"

    elif unresolved:

        accession_match = "UNRESOLVED"

    else:

        accession_match = "MISSING_METADATA"


    # --------------------------------------------------------
    # LENGTH QC
    # --------------------------------------------------------

    length_status = "NOT_CHECKED"

    if accepted_length:

        try:

            metadata_length = int(
                float(
                    accepted_length
                )
            )

        except Exception:

            metadata_length = None


        if metadata_length is None:

            length_status = (
                "INVALID_METADATA_LENGTH"
            )

        elif sequence_count == 0:

            length_status = "NO_SEQUENCE"

        else:

            observed_lengths = {
                int(
                    x["Sequence_length"]
                )
                for x in matches
            }

            if metadata_length in observed_lengths:

                length_status = "MATCH"

            else:

                length_status = "MISMATCH"

    elif unresolved:

        length_status = "UNRESOLVED"

    else:

        length_status = "MISSING_METADATA"


    # --------------------------------------------------------
    # SEQUENCE TYPE QC
    # --------------------------------------------------------

    if unresolved:

        sequence_type_status = "UNRESOLVED"

    elif sequence_count == 0:

        sequence_type_status = "NO_SEQUENCE"

    elif accepted_sequence_type:

        sequence_type_status = "METADATA_AVAILABLE"

    else:

        sequence_type_status = "MISSING_METADATA"


    # --------------------------------------------------------
    # OVERALL QC
    # --------------------------------------------------------

    if unresolved:

        overall_status = (
            "UNRESOLVED_MASTER_ONLY"
        )

    elif sequence_count == 0:

        overall_status = "NO_SEQUENCE"

    elif accession_match == "MISMATCH":

        overall_status = "ACCESSION_MISMATCH"

    elif length_status == "MISMATCH":

        overall_status = "LENGTH_MISMATCH"

    else:

        overall_status = "PASS"


    rows.append({

        "KJU_ID": kju,

        "Master_status":
            clean(m["Master_status"]),

        "Adjudication_type":
            clean(m["Adjudication_type"]),

        "Accepted_taxonomy":
            accepted_taxonomy,

        "Accepted_accession":
            clean(m["Accepted_accession"]),

        "Normalized_accepted_accession":
            accepted_acc,

        "Accepted_sequence_type":
            accepted_sequence_type,

        "Accepted_length":
            accepted_length,

        "Sequence_count":
            sequence_count,

        "Sequence_status":
            sequence_status,

        "Detected_accessions":
            detected_accessions,

        "Detected_lengths":
            detected_lengths,

        "Detected_sequence_types":
            detected_types,

        "Sequence_hashes":
            sequence_hashes,

        "Accession_match":
            accession_match,

        "Length_status":
            length_status,

        "Sequence_type_status":
            sequence_type_status,

        "Overall_QC":
            overall_status,

    })


qc_df = pd.DataFrame(rows)


# ============================================================
# SAVE INVENTORY
# ============================================================

qc_df.to_csv(
    INVENTORY_FILE,
    index=False
)


inv.to_csv(
    OUTDIR / "STEP77_RAW_SEQUENCE_RECORDS.csv",
    index=False
)


# ============================================================
# MISSING SEQUENCES
# ============================================================

missing = qc_df[
    qc_df["Sequence_count"] == 0
].copy()


missing.to_csv(
    MISSING_FILE,
    index=False
)


# ============================================================
# LENGTH MISMATCHES
# ============================================================

length_mismatch = qc_df[
    qc_df["Length_status"] == "MISMATCH"
].copy()


length_mismatch.to_csv(
    MISMATCH_FILE,
    index=False
)


# ============================================================
# DUPLICATE SEQUENCE QC
# ============================================================

if len(inv):

    duplicate_groups = (
        inv
        .groupby(
            "Sequence_SHA256"
        )
        .agg(
            Sequence_count=(
                "Sequence_SHA256",
                "size"
            ),
            KJU_IDs=(
                "KJU_ID",
                lambda x:
                ";".join(
                    sorted(
                        set(
                            clean(v)
                            for v in x
                            if clean(v)
                        )
                    )
                )
            ),
            Accessions=(
                "Normalized_accession_detected",
                lambda x:
                ";".join(
                    sorted(
                        set(
                            clean(v)
                            for v in x
                            if clean(v)
                        )
                    )
                )
            ),
            Lengths=(
                "Sequence_length",
                lambda x:
                ";".join(
                    sorted(
                        set(
                            str(v)
                            for v in x
                        )
                    )
                )
            )
        )
        .reset_index()
    )


    duplicate_groups = duplicate_groups[
        duplicate_groups["Sequence_count"] > 1
    ].copy()

else:

    duplicate_groups = pd.DataFrame(
        columns=[
            "Sequence_SHA256",
            "Sequence_count",
            "KJU_IDs",
            "Accessions",
            "Lengths",
        ]
    )


duplicate_groups.to_csv(
    DUP_FILE,
    index=False
)


# ============================================================
# SUMMARY COUNTS
# ============================================================

total_kju = len(qc_df)

sequence_found = (
    qc_df["Sequence_count"] > 0
).sum()

no_sequence = (
    qc_df["Sequence_count"] == 0
).sum()

resolved = (
    qc_df["Adjudication_type"]
    != "MASTER_ONLY_UNRESOLVED"
).sum()

unresolved = (
    qc_df["Adjudication_type"]
    == "MASTER_ONLY_UNRESOLVED"
).sum()

resolved_with_sequence = (
    (
        qc_df["Adjudication_type"]
        != "MASTER_ONLY_UNRESOLVED"
    )
    &
    (
        qc_df["Sequence_count"] > 0
    )
).sum()

resolved_without_sequence = (
    (
        qc_df["Adjudication_type"]
        != "MASTER_ONLY_UNRESOLVED"
    )
    &
    (
        qc_df["Sequence_count"] == 0
    )
).sum()

accession_matches = (
    qc_df["Accession_match"]
    == "MATCH"
).sum()

accession_mismatches = (
    qc_df["Accession_match"]
    == "MISMATCH"
).sum()

length_matches = (
    qc_df["Length_status"]
    == "MATCH"
).sum()

length_mismatches = (
    qc_df["Length_status"]
    == "MISMATCH"
).sum()

unresolved_with_sequence = (
    (
        qc_df["Adjudication_type"]
        == "MASTER_ONLY_UNRESOLVED"
    )
    &
    (
        qc_df["Sequence_count"] > 0
    )
).sum()


# ============================================================
# REPORT
# ============================================================

lines = []

lines.append(
    "STEP 77 — SEQUENCE INVENTORY + METADATA CONSISTENCY QC"
)

lines.append("=" * 110)
lines.append("")

lines.append("1. INPUT")
lines.append("-" * 100)
lines.append(str(MASTER_FILE))
lines.append("")

lines.append("2. MASTER COUNTS")
lines.append("-" * 100)
lines.append(
    f"Total KJU records              : {total_kju}"
)
lines.append(
    f"Resolved records               : {resolved}"
)
lines.append(
    f"Unresolved MASTER_ONLY         : {unresolved}"
)
lines.append("")

lines.append("3. SEQUENCE FILE INVENTORY")
lines.append("-" * 100)
lines.append(
    f"Sequence files discovered      : {len(files)}"
)
lines.append(
    f"Raw sequence records           : {len(inv)}"
)
lines.append(
    f"KJUs with sequence             : {sequence_found}"
)
lines.append(
    f"KJUs without sequence          : {no_sequence}"
)
lines.append("")

lines.append("4. RESOLVED RECORD SEQUENCE QC")
lines.append("-" * 100)
lines.append(
    f"Resolved with sequence         : "
    f"{resolved_with_sequence}"
)
lines.append(
    f"Resolved without sequence      : "
    f"{resolved_without_sequence}"
)
lines.append("")

lines.append("5. ACCESSION QC")
lines.append("-" * 100)
lines.append(
    f"Accession matches              : "
    f"{accession_matches}"
)
lines.append(
    f"Accession mismatches           : "
    f"{accession_mismatches}"
)
lines.append("")

lines.append("6. LENGTH QC")
lines.append("-" * 100)
lines.append(
    f"Length matches                 : "
    f"{length_matches}"
)
lines.append(
    f"Length mismatches              : "
    f"{length_mismatches}"
)
lines.append("")

lines.append("7. UNRESOLVED RECORD QC")
lines.append("-" * 100)
lines.append(
    f"Unresolved records             : "
    f"{unresolved}"
)
lines.append(
    f"Unresolved with sequence       : "
    f"{unresolved_with_sequence}"
)
lines.append(
    "NOTE: sequence presence does not automatically "
    "resolve taxonomy."
)
lines.append("")

lines.append("8. DUPLICATE SEQUENCE QC")
lines.append("-" * 100)
lines.append(
    f"Duplicate sequence groups      : "
    f"{len(duplicate_groups)}"
)
lines.append("")

lines.append("9. LENGTH MISMATCH FILE")
lines.append("-" * 100)
lines.append(str(MISMATCH_FILE))
lines.append("")

lines.append("10. MISSING SEQUENCE FILE")
lines.append("-" * 100)
lines.append(str(MISSING_FILE))
lines.append("")

lines.append("11. DUPLICATE SEQUENCE FILE")
lines.append("-" * 100)
lines.append(str(DUP_FILE))
lines.append("")

lines.append("12. FINAL STEP 77 STATUS")
lines.append("-" * 100)

if accession_mismatches == 0:
    lines.append(
        "PASS — no accession mismatches detected"
    )
else:
    lines.append(
        f"REVIEW — {accession_mismatches} accession mismatches detected"
    )

if length_mismatches == 0:
    lines.append(
        "PASS — no sequence length mismatches detected"
    )
else:
    lines.append(
        f"REVIEW — {length_mismatches} sequence length mismatches detected"
    )

if resolved_without_sequence == 0:
    lines.append(
        "PASS — every resolved record has a sequence"
    )
else:
    lines.append(
        f"REVIEW — {resolved_without_sequence} resolved records lack sequences"
    )

lines.append("")
lines.append(
    "Frozen V1 master was NOT modified."
)
lines.append("")


QC_FILE.write_text(
    "\n".join(lines) + "\n"
)


# ============================================================
# CONSOLE REPORT
# ============================================================

print()
print("=" * 110)
print("STEP 77 COMPLETE")
print("=" * 110)

print()
print("MASTER")
print("-" * 100)
print(
    f"Total KJU records       : {total_kju}"
)
print(
    f"Resolved               : {resolved}"
)
print(
    f"Unresolved             : {unresolved}"
)

print()
print("SEQUENCE INVENTORY")
print("-" * 100)
print(
    f"Sequence files         : {len(files)}"
)
print(
    f"Raw sequence records   : {len(inv)}"
)
print(
    f"KJUs with sequence     : {sequence_found}"
)
print(
    f"KJUs without sequence  : {no_sequence}"
)

print()
print("ACCESSION QC")
print("-" * 100)
print(
    f"Matches                : {accession_matches}"
)
print(
    f"Mismatches             : {accession_mismatches}"
)

print()
print("LENGTH QC")
print("-" * 100)
print(
    f"Matches                : {length_matches}"
)
print(
    f"Mismatches             : {length_mismatches}"
)

print()
print("RESOLVED RECORDS")
print("-" * 100)
print(
    f"With sequence          : {resolved_with_sequence}"
)
print(
    f"Without sequence       : {resolved_without_sequence}"
)

print()
print("DUPLICATE SEQUENCES")
print("-" * 100)
print(
    f"Duplicate groups       : {len(duplicate_groups)}"
)

print()
print("OUTPUT FILES")
print("-" * 100)
print(
    f"Inventory              : {INVENTORY_FILE}"
)
print(
    f"Raw sequences          : "
    f"{OUTDIR / 'STEP77_RAW_SEQUENCE_RECORDS.csv'}"
)
print(
    f"Missing sequences      : {MISSING_FILE}"
)
print(
    f"Length mismatches      : {MISMATCH_FILE}"
)
print(
    f"Duplicate sequences    : {DUP_FILE}"
)
print(
    f"QC report              : {QC_FILE}"
)

print()
print("=" * 110)
print("FROZEN V1 MASTER WAS NOT MODIFIED")
print("=" * 110)
