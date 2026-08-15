from pathlib import Path
import pandas as pd
import re
from collections import Counter

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

FINAL = META / "FINAL_KJU_16S_DATASET_VALIDATED_V1.csv"
QC_DIR = META / "FINAL_109_SEQUENCE_QC_V2"

QC_DIR.mkdir(exist_ok=True)

FASTA_OUT = QC_DIR / "FINAL_109_KJU_16S_QC_V2.fasta"
QC_OUT = QC_DIR / "FINAL_109_SEQUENCE_QC_V2.csv"
SUMMARY = QC_DIR / "FINAL_109_SEQUENCE_QC_SUMMARY_V2.txt"
DUP_OUT = QC_DIR / "DUPLICATE_SEQUENCE_GROUPS_V2.csv"

print("=" * 110)
print(" DODDAGGUBBI - IUPAC-AWARE 109 SEQUENCE QC")
print("=" * 110)

df = pd.read_csv(FINAL, dtype=str).fillna("")

# Standard IUPAC DNA symbols
IUPAC_DNA = set("ACGTRYSWKMBDHVN")

# Unambiguous DNA
ACGT = set("ACGT")

rows = []
fasta_records = []

for _, row in df.iterrows():

    kju = row["kju_id"].strip()
    accession = row["accession"].strip()
    organism = row["organism"].strip()

    gb_file = Path(row["genbank_file"].strip())

    if not gb_file.is_absolute():
        gb_file = MASTER / gb_file

    if not gb_file.exists():
        raise SystemExit(
            f"ERROR: Missing GenBank file for {kju}: {gb_file}"
        )

    text = gb_file.read_text(
        encoding="utf-8",
        errors="ignore"
    )

    # Extract sequence section from GenBank
    m = re.search(
        r"ORIGIN(.*?)(?=//)",
        text,
        flags=re.DOTALL
    )

    if not m:
        raise SystemExit(
            f"ERROR: No ORIGIN sequence found for {kju}"
        )

    origin = m.group(1)

    seq = "".join(
        re.findall(
            r"[A-Za-z]",
            origin
        )
    ).upper()

    if not seq:
        raise SystemExit(
            f"ERROR: Empty sequence for {kju}"
        )

    invalid = sorted(
        set(seq) - IUPAC_DNA
    )

    ambiguity_positions = [
        i + 1
        for i, base in enumerate(seq)
        if base not in ACGT
    ]

    ambiguity_counts = Counter(
        base
        for base in seq
        if base not in ACGT
    )

    ambiguity_count = sum(
        ambiguity_counts.values()
    )

    gc = (
        (seq.count("G") + seq.count("C"))
        / len(seq)
        * 100
    )

    problems = []

    if invalid:
        problems.append(
            "INVALID_NON_IUPAC_CHARACTER"
        )

    qc_status = (
        "PASS"
        if not problems
        else "FAIL"
    )

    rows.append({
        "kju_id": kju,
        "accession": accession,
        "organism": organism,
        "sequence_length": len(seq),
        "gc_percent": round(gc, 4),
        "iupac_ambiguity_count": ambiguity_count,
        "iupac_ambiguity_percent": round(
            ambiguity_count / len(seq) * 100,
            4
        ),
        "ambiguity_symbols":
            ";".join(
                f"{k}:{v}"
                for k, v in sorted(
                    ambiguity_counts.items()
                )
            ),
        "ambiguity_positions":
            ";".join(
                map(str, ambiguity_positions[:100])
            ),
        "invalid_characters":
            "".join(invalid),
        "qc_status": qc_status,
        "qc_problems":
            ";".join(problems),
        "genbank_file":
            str(gb_file),
    })

    fasta_records.append(
        (
            f">{kju}|{accession}|{organism.replace(' ', '_')}",
            seq
        )
    )

qc = pd.DataFrame(rows)

# ------------------------------------------------------------
# Duplicate sequence detection
# ------------------------------------------------------------

qc["_seq_key"] = [
    seq
    for _, seq in fasta_records
]

duplicate_rows = []

for seq, g in qc.groupby("_seq_key"):

    if len(g) > 1:

        duplicate_rows.append({
            "sequence_length": len(seq),
            "kju_ids": ";".join(g["kju_id"]),
            "accessions": ";".join(g["accession"]),
            "count": len(g),
        })

duplicates = pd.DataFrame(duplicate_rows)

if duplicates.empty:
    duplicates = pd.DataFrame(
        columns=[
            "sequence_length",
            "kju_ids",
            "accessions",
            "count",
        ]
    )

duplicates.to_csv(
    DUP_OUT,
    index=False
)

qc = qc.drop(columns="_seq_key")

qc.to_csv(
    QC_OUT,
    index=False
)

# ------------------------------------------------------------
# FASTA
# ------------------------------------------------------------

with open(
    FASTA_OUT,
    "w",
    encoding="utf-8"
) as fh:

    for header, seq in fasta_records:

        fh.write(header + "\n")

        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

pass_count = (
    qc["qc_status"] == "PASS"
).sum()

fail_count = (
    qc["qc_status"] == "FAIL"
).sum()

ambiguous_records = (
    qc["iupac_ambiguity_count"] > 0
).sum()

summary = []

summary.append(
    "DODDAGGUBBI - IUPAC-AWARE 109 SEQUENCE QC"
)

summary.append("=" * 110)
summary.append("")

summary.append(
    f"Records checked: {len(qc)}"
)

summary.append(
    f"PASS: {pass_count}"
)

summary.append(
    f"FAIL: {fail_count}"
)

summary.append("")

summary.append("LENGTH")
summary.append("-" * 110)

summary.append(
    f"Minimum: {qc['sequence_length'].min()}"
)

summary.append(
    f"Maximum: {qc['sequence_length'].max()}"
)

summary.append(
    f"Mean: {qc['sequence_length'].mean():.2f}"
)

summary.append(
    f"Median: {qc['sequence_length'].median():.2f}"
)

summary.append("")

summary.append("IUPAC AMBIGUITY")
summary.append("-" * 110)

summary.append(
    f"Sequences with IUPAC ambiguity: "
    f"{ambiguous_records}"
)

summary.append(
    f"Total ambiguity symbols: "
    f"{qc['iupac_ambiguity_count'].sum()}"
)

summary.append("")

summary.append("DUPLICATE SEQUENCES")
summary.append("-" * 110)

summary.append(
    f"Duplicate groups: {len(duplicates)}"
)

summary.append(
    f"Records in duplicate groups: "
    f"{sum(duplicates['count']) if len(duplicates) else 0}"
)

summary.append("")

summary.append(
    "INTERPRETATION"
)

summary.append("-" * 110)

summary.append(
    "IUPAC ambiguity symbols are valid nucleotide notation "
    "and are reported separately from invalid characters."
)

summary.append(
    "Identical 16S sequences across different accessions are "
    "retained unless an independent identity conflict is demonstrated."
)

summary.append("")

summary.append("KJU0095")
summary.append("-" * 110)

k95 = qc[qc["kju_id"] == "KJU0095"]

if len(k95):
    r = k95.iloc[0]

    summary.append(
        f"Length: {r['sequence_length']}"
    )

    summary.append(
        f"Ambiguity: {r['ambiguity_symbols']}"
    )

    summary.append(
        f"Status: {r['qc_status']}"
    )

summary.append("")

summary.append(
    "No source GenBank files were modified."
)

summary.append(
    "Original master was not modified."
)

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

# ------------------------------------------------------------
# Console
# ------------------------------------------------------------

print("")
print("=" * 110)
print("IUPAC-AWARE QC RESULT")
print("=" * 110)

print("Records:", len(qc))
print("PASS:", pass_count)
print("FAIL:", fail_count)

print("")
print(
    "Sequences with IUPAC ambiguity:",
    ambiguous_records
)

print(
    "Total ambiguity symbols:",
    qc["iupac_ambiguity_count"].sum()
)

print(
    "Duplicate groups:",
    len(duplicates)
)

print("")
print("KJU0095:")
print(
    qc[
        qc["kju_id"] == "KJU0095"
    ].to_string(index=False)
)

print("")
print("FASTA:", FASTA_OUT)
print("QC:", QC_OUT)
print("SUMMARY:", SUMMARY)
print("DUPLICATES:", DUP_OUT)

if fail_count:
    raise SystemExit(
        "\nQC STILL HAS FAILURES."
    )

print("")
print("ALL 109 SEQUENCES PASS IUPAC-AWARE QC.")
