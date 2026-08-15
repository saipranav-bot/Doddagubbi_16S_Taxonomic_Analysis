from pathlib import Path
import pandas as pd
import re
import hashlib
import sys

from Bio import SeqIO

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

FINAL = META / "FINAL_KJU_16S_DATASET_VALIDATED_V1.csv"
VALIDATION = META / "FINAL_109_GENBANK_VALIDATION_V1.csv"

OUTDIR = META / "FINAL_109_SEQUENCE_QC_V1"
OUTDIR.mkdir(exist_ok=True)

FASTA_OUT = OUTDIR / "FINAL_109_KJU_16S.fasta"
QC_OUT = OUTDIR / "FINAL_109_SEQUENCE_QC.csv"
SUMMARY_OUT = OUTDIR / "FINAL_109_SEQUENCE_QC_SUMMARY.txt"
DUP_OUT = OUTDIR / "DUPLICATE_SEQUENCE_GROUPS.csv"

print("=" * 110)
print(" DODDAGGUBBI - 109 SEQUENCE EXTRACTION + QC")
print("=" * 110)

# ------------------------------------------------------------
# Load validated metadata
# ------------------------------------------------------------

if not FINAL.exists():
    raise SystemExit(f"ERROR: Missing:\n{FINAL}")

if not VALIDATION.exists():
    raise SystemExit(f"ERROR: Missing:\n{VALIDATION}")

df = pd.read_csv(FINAL, dtype=str).fillna("")
validation = pd.read_csv(VALIDATION, dtype=str).fillna("")

if len(df) != 109:
    raise SystemExit(
        f"ERROR: Expected 109 validated records, found {len(df)}"
    )

# Only use records that passed GenBank validation
passed = validation[
    validation["validation_status"] == "PASS"
].copy()

if len(passed) != 109:
    raise SystemExit(
        f"ERROR: Expected 109 validation PASS records, found {len(passed)}"
    )

# ------------------------------------------------------------
# Merge validation information
# ------------------------------------------------------------

df = df.merge(
    passed[
        [
            "kju_id",
            "accession",
            "genbank_file",
            "genbank_accession",
            "genbank_length",
            "genbank_organism",
        ]
    ],
    on=["kju_id", "accession"],
    how="left",
    suffixes=("", "_validated"),
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

VALID_DNA = set("ACGTN")
STRICT_DNA = set("ACGT")


def normalize_sequence(seq):
    return re.sub(r"\s+", "", str(seq).upper())


def sha256_sequence(seq):
    return hashlib.sha256(
        seq.encode("utf-8")
    ).hexdigest()


# ------------------------------------------------------------
# Extract sequences
# ------------------------------------------------------------

qc_rows = []
fasta_records = []

for _, row in df.iterrows():

    kju = row["kju_id"].strip()
    accession = row["accession"].strip()
    organism = row["organism"].strip()
    gb_file = row["genbank_file"].strip()

    if not gb_file:
        raise SystemExit(
            f"ERROR: No GenBank file recorded for {kju}"
        )

    gb_path = Path(gb_file)

    if not gb_path.is_absolute():
        gb_path = MASTER / gb_path

    if not gb_path.exists():
        raise SystemExit(
            f"ERROR: GenBank file missing for {kju}:\n{gb_path}"
        )

    records = list(
        SeqIO.parse(
            gb_path,
            "genbank"
        )
    )

    if len(records) != 1:
        raise SystemExit(
            f"ERROR: Expected exactly 1 GenBank record for {kju}, "
            f"found {len(records)}:\n{gb_path}"
        )

    record = records[0]

    seq = normalize_sequence(
        str(record.seq)
    )

    if not seq:
        raise SystemExit(
            f"ERROR: Empty sequence for {kju}"
        )

    length = len(seq)

    invalid_chars = sorted(
        set(seq) - VALID_DNA
    )

    ambiguous_count = seq.count("N")

    gc_count = (
        seq.count("G") +
        seq.count("C")
    )

    gc_content = (
        gc_count / length * 100
        if length
        else 0
    )

    sequence_hash = sha256_sequence(seq)

    kju_in_description = (
        kju.upper()
        in record.description.upper()
    )

    accession_match = (
        accession.split(".")[0].upper()
        ==
        record.id.split(".")[0].upper()
    )

    strict_dna = not bool(
        set(seq) - STRICT_DNA
    )

    # QC status
    problems = []

    if invalid_chars:
        problems.append(
            "INVALID_DNA_CHARACTER"
        )

    if not accession_match:
        problems.append(
            "ACCESSION_MISMATCH"
        )

    if not kju_in_description:
        problems.append(
            "KJU_NOT_IN_DESCRIPTION"
        )

    status = (
        "PASS"
        if not problems
        else "FAIL"
    )

    fasta_header = (
        f">{kju}|{accession}|{organism.replace(' ', '_')}"
    )

    fasta_records.append(
        (
            fasta_header,
            seq
        )
    )

    qc_rows.append({

        "kju_id": kju,

        "accession": accession,

        "organism": organism,

        "genbank_file": str(gb_path),

        "record_id": record.id,

        "description": record.description,

        "sequence_length": length,

        "ambiguous_N_count": ambiguous_count,

        "ambiguous_N_percent": round(
            ambiguous_count / length * 100,
            4
        ),

        "gc_percent": round(
            gc_content,
            4
        ),

        "invalid_characters": "".join(
            invalid_chars
        ),

        "strict_ACGT_only": strict_dna,

        "accession_match": accession_match,

        "kju_in_description": kju_in_description,

        "sequence_sha256": sequence_hash,

        "qc_status": status,

        "qc_problems": ";".join(
            problems
        ),
    })


# ------------------------------------------------------------
# Create QC table
# ------------------------------------------------------------

qc = pd.DataFrame(qc_rows)

qc = qc.sort_values(
    "kju_id",
    key=lambda s:
    s.str.extract(
        r"(\d+)",
        expand=False
    ).astype(int)
)

qc.to_csv(
    QC_OUT,
    index=False
)

# ------------------------------------------------------------
# Write FASTA
# ------------------------------------------------------------

with open(
    FASTA_OUT,
    "w",
    encoding="utf-8"
) as fh:

    for header, seq in fasta_records:

        fh.write(header + "\n")

        for i in range(0, len(seq), 80):
            fh.write(
                seq[i:i+80] + "\n"
            )

# ------------------------------------------------------------
# Duplicate sequence detection
# ------------------------------------------------------------

duplicate_groups = []

for seq_hash, group in qc.groupby(
    "sequence_sha256"
):

    if len(group) > 1:

        duplicate_groups.append({

            "sequence_sha256": seq_hash,

            "count": len(group),

            "kju_ids": ";".join(
                group["kju_id"]
            ),

            "accessions": ";".join(
                group["accession"]
            ),
        })

duplicates = pd.DataFrame(
    duplicate_groups
)

if duplicates.empty:

    duplicates = pd.DataFrame(
        columns=[
            "sequence_sha256",
            "count",
            "kju_ids",
            "accessions",
        ]
    )

duplicates.to_csv(
    DUP_OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

pass_count = (
    qc["qc_status"] == "PASS"
).sum()

fail_count = (
    qc["qc_status"] == "FAIL"
).sum()

duplicate_sequence_count = (
    qc["sequence_sha256"]
    .duplicated(keep=False)
    .sum()
)

summary = []

summary.append(
    "DODDAGGUBBI - FINAL 109 SEQUENCE QC"
)

summary.append("=" * 110)
summary.append("")

summary.append(
    f"Validated input records: {len(qc)}"
)

summary.append(
    f"Sequence QC PASS: {pass_count}"
)

summary.append(
    f"Sequence QC FAIL: {fail_count}"
)

summary.append("")

summary.append("LENGTH")
summary.append("-" * 110)

summary.append(
    f"Minimum length: {qc['sequence_length'].min()}"
)

summary.append(
    f"Maximum length: {qc['sequence_length'].max()}"
)

summary.append(
    f"Mean length: {qc['sequence_length'].mean():.2f}"
)

summary.append(
    f"Median length: {qc['sequence_length'].median():.2f}"
)

summary.append("")

summary.append("AMBIGUOUS BASES")
summary.append("-" * 110)

summary.append(
    f"Sequences containing N: "
    f"{(qc['ambiguous_N_count'] > 0).sum()}"
)

summary.append(
    f"Total ambiguous bases: "
    f"{qc['ambiguous_N_count'].sum()}"
)

summary.append(
    f"Maximum N percentage: "
    f"{qc['ambiguous_N_percent'].max():.4f}%"
)

summary.append("")

summary.append("GC CONTENT")
summary.append("-" * 110)

summary.append(
    f"Minimum GC%: {qc['gc_percent'].min():.2f}"
)

summary.append(
    f"Maximum GC%: {qc['gc_percent'].max():.2f}"
)

summary.append(
    f"Mean GC%: {qc['gc_percent'].mean():.2f}"
)

summary.append("")

summary.append("DUPLICATES")
summary.append("-" * 110)

summary.append(
    f"Duplicate sequence groups: "
    f"{len(duplicates)}"
)

summary.append(
    f"Records belonging to duplicate groups: "
    f"{duplicate_sequence_count}"
)

summary.append("")

summary.append("VALIDATION")
summary.append("-" * 110)

summary.append(
    f"Accession matches: "
    f"{qc['accession_match'].sum()} / {len(qc)}"
)

summary.append(
    f"KJU found in description: "
    f"{qc['kju_in_description'].sum()} / {len(qc)}"
)

summary.append(
    f"Invalid-character sequences: "
    f"{(qc['invalid_characters'] != '').sum()}"
)

summary.append("")

summary.append("FAILED SEQUENCES")
summary.append("-" * 110)

failed = qc[
    qc["qc_status"] == "FAIL"
]

if failed.empty:
    summary.append("NONE")
else:
    for _, r in failed.iterrows():
        summary.append(
            f"{r['kju_id']} | "
            f"{r['accession']} | "
            f"{r['qc_problems']}"
        )

summary.append("")
summary.append("IMPORTANT")
summary.append("-" * 110)

summary.append(
    "Source GenBank files were NOT modified."
)

summary.append(
    "Original master was NOT modified."
)

summary.append(
    "KJU0039 is not part of this 109-record validated sequence set."
)

summary.append("")

summary.append(
    f"FASTA: {FASTA_OUT}"
)

summary.append(
    f"QC table: {QC_OUT}"
)

summary.append(
    f"Duplicate groups: {DUP_OUT}"
)

SUMMARY_OUT.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

print("")
print("=" * 110)
print(" SEQUENCE QC COMPLETE")
print("=" * 110)

print("Records:", len(qc))
print("PASS:", pass_count)
print("FAIL:", fail_count)

print("")
print("Length:")
print("  Min:", qc["sequence_length"].min())
print("  Max:", qc["sequence_length"].max())
print("  Mean:", round(qc["sequence_length"].mean(), 2))
print("  Median:", qc["sequence_length"].median())

print("")
print("Sequences containing N:",
      (qc["ambiguous_N_count"] > 0).sum())

print("Duplicate sequence groups:",
      len(duplicates))

print("")
print("FASTA:", FASTA_OUT)
print("QC:", QC_OUT)
print("SUMMARY:", SUMMARY_OUT)
print("DUPLICATES:", DUP_OUT)

if fail_count:
    print("")
    print("QC FAILURES DETECTED.")
    sys.exit(1)

print("")
print("ALL 109 SEQUENCES PASSED BASIC SEQUENCE QC.")
