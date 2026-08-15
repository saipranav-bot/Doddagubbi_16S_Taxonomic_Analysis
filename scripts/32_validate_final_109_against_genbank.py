from pathlib import Path
import pandas as pd
import re
import sys

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"
GB_DIR = MASTER / "01_NCBI" / "genbank"

FINAL = META / "FINAL_KJU_16S_DATASET_VALIDATED_V1.csv"

OUT = META / "FINAL_109_GENBANK_VALIDATION_V1.csv"
SUMMARY = META / "FINAL_109_GENBANK_VALIDATION_SUMMARY_V1.txt"

print("=" * 110)
print(" DODDAGGUBBI - FINAL 109 DATASET ↔ GENBANK VALIDATION")
print("=" * 110)

if not FINAL.exists():
    raise SystemExit(f"ERROR: Missing:\n{FINAL}")

df = pd.read_csv(FINAL, dtype=str).fillna("")

required = {
    "kju_id",
    "accession",
    "organism",
    "sequence_type",
    "length",
    "genbank_file",
}

missing = required - set(df.columns)

if missing:
    raise SystemExit(
        f"ERROR: Missing columns: {sorted(missing)}"
    )

print("Final dataset rows:", len(df))

if len(df) != 109:
    raise SystemExit(
        f"ERROR: Expected 109 records, found {len(df)}"
    )

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def norm_acc(x):
    x = str(x).strip().upper()
    if not x:
        return ""
    return x.split(".")[0]


def norm_org(x):
    x = str(x).strip().lower()
    x = re.sub(r"\s+", " ", x)
    return x


def parse_genbank(path):
    text = path.read_text(
        encoding="utf-8",
        errors="ignore"
    )

    version_match = re.search(
        r"^VERSION\s+(\S+)",
        text,
        flags=re.MULTILINE
    )

    accession_match = re.search(
        r"^ACCESSION\s+(\S+)",
        text,
        flags=re.MULTILINE
    )

    locus_match = re.search(
        r"^LOCUS\s+\S+\s+(\d+)\s+bp",
        text,
        flags=re.MULTILINE
    )

    definition_match = re.search(
        r"^DEFINITION\s+(.+?)(?=\n[A-Z]|\Z)",
        text,
        flags=re.MULTILINE | re.DOTALL
    )

    organism_match = re.search(
        r"^\s+ORGANISM\s+(.+)$",
        text,
        flags=re.MULTILINE
    )

    kju_hits = sorted(set(
        x.upper()
        for x in re.findall(
            r"\bKJU\d{4}\b",
            text,
            flags=re.IGNORECASE
        )
    ))

    return {
        "accession": (
            version_match.group(1)
            if version_match
            else (
                accession_match.group(1)
                if accession_match
                else ""
            )
        ),
        "length": (
            locus_match.group(1)
            if locus_match
            else ""
        ),
        "definition": (
            " ".join(
                definition_match.group(1)
                .split()
            )
            if definition_match
            else ""
        ),
        "organism": (
            organism_match.group(1).strip()
            if organism_match
            else ""
        ),
        "kju_hits": kju_hits,
    }


# ------------------------------------------------------------
# Validate every record
# ------------------------------------------------------------

rows = []

for _, r in df.iterrows():

    kju = r["kju_id"].strip()
    accession = r["accession"].strip()
    organism = r["organism"].strip()
    expected_length = r["length"].strip()
    gb_file = r["genbank_file"].strip()

    problems = []

    if not gb_file:
        problems.append("NO_GENBANK_FILE")

        rows.append({
            "kju_id": kju,
            "accession": accession,
            "organism": organism,
            "genbank_file": gb_file,
            "genbank_exists": False,
            "genbank_accession": "",
            "genbank_length": "",
            "genbank_organism": "",
            "kju_in_genbank": False,
            "accession_match": False,
            "organism_match": False,
            "length_match": False,
            "validation_status": "FAIL",
            "problems": ";".join(problems),
        })

        continue

    gb_path = MASTER / gb_file

    # also permit direct path if stored absolute
    if not gb_path.exists():
        gb_path = Path(gb_file)

    if not gb_path.exists():
        problems.append("GENBANK_FILE_MISSING")

        rows.append({
            "kju_id": kju,
            "accession": accession,
            "organism": organism,
            "genbank_file": gb_file,
            "genbank_exists": False,
            "genbank_accession": "",
            "genbank_length": "",
            "genbank_organism": "",
            "kju_in_genbank": False,
            "accession_match": False,
            "organism_match": False,
            "length_match": False,
            "validation_status": "FAIL",
            "problems": ";".join(problems),
        })

        continue

    gb = parse_genbank(gb_path)

    expected_acc_base = norm_acc(accession)
    gb_acc_base = norm_acc(gb["accession"])

    accession_match = (
        expected_acc_base != ""
        and expected_acc_base == gb_acc_base
    )

    kju_match = kju in gb["kju_hits"]

    # direct organism comparison
    organism_match = (
        norm_org(organism)
        == norm_org(gb["organism"])
    )

    length_match = (
        expected_length != ""
        and gb["length"] != ""
        and expected_length == gb["length"]
    )

    if not accession_match:
        problems.append("ACCESSION_MISMATCH")

    if not kju_match:
        problems.append("KJU_NOT_IN_GENBANK")

    if not organism_match:
        problems.append("ORGANISM_MISMATCH")

    if not length_match:
        problems.append("LENGTH_MISMATCH")

    status = (
        "PASS"
        if not problems
        else "FAIL"
    )

    rows.append({
        "kju_id": kju,
        "accession": accession,
        "organism": organism,
        "genbank_file": str(gb_path),
        "genbank_exists": True,
        "genbank_accession": gb["accession"],
        "genbank_length": gb["length"],
        "genbank_organism": gb["organism"],
        "kju_in_genbank": kju_match,
        "accession_match": accession_match,
        "organism_match": organism_match,
        "length_match": length_match,
        "validation_status": status,
        "problems": ";".join(problems),
    })


result = pd.DataFrame(rows)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

result.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

passes = (
    result["validation_status"] == "PASS"
).sum()

fails = len(result) - passes

summary = []

summary.append(
    "DODDAGGUBBI - FINAL 109 DATASET ↔ GENBANK VALIDATION"
)

summary.append("=" * 110)
summary.append("")

summary.append(f"Records checked: {len(result)}")
summary.append(f"PASS: {passes}")
summary.append(f"FAIL: {fails}")

summary.append("")
summary.append("FIELD CHECKS")
summary.append("-" * 110)

for col in [
    "accession_match",
    "kju_in_genbank",
    "organism_match",
    "length_match",
]:
    summary.append(
        f"{col}: {result[col].sum()} / {len(result)}"
    )

summary.append("")
summary.append("FAILED RECORDS")
summary.append("-" * 110)

failed = result[
    result["validation_status"] == "FAIL"
]

if failed.empty:
    summary.append("NONE")
else:
    for _, r in failed.iterrows():
        summary.append(
            f"{r['kju_id']} | "
            f"{r['accession']} | "
            f"{r['problems']}"
        )

summary.append("")
summary.append("CRITICAL KJU0038 CHECK")
summary.append("-" * 110)

k38 = result[
    result["kju_id"] == "KJU0038"
]

if len(k38) == 1:
    r = k38.iloc[0]

    summary.append(
        f"accession={r['accession']}"
    )
    summary.append(
        f"genbank_accession={r['genbank_accession']}"
    )
    summary.append(
        f"kju_in_genbank={r['kju_in_genbank']}"
    )
    summary.append(
        f"organism_match={r['organism_match']}"
    )
    summary.append(
        f"status={r['validation_status']}"
    )

summary.append("")
summary.append("KJU0039 CHECK")
summary.append("-" * 110)

k39 = result[
    result["kju_id"] == "KJU0039"
]

if len(k39) == 0:
    summary.append(
        "KJU0039 is correctly absent from the 109 sequence dataset."
    )
else:
    summary.append(
        "WARNING: KJU0039 unexpectedly appears in the 109 sequence dataset."
    )

summary.append("")
summary.append("IMPORTANT")
summary.append("-" * 110)
summary.append(
    "No sequence or metadata files were modified."
)
summary.append(
    "This validation uses only the actual GenBank file associated "
    "with each final record."
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
print("GENBANK VALIDATION RESULT")
print("=" * 110)

print("Records:", len(result))
print("PASS:", passes)
print("FAIL:", fails)

print("")
print("Field checks:")
for col in [
    "accession_match",
    "kju_in_genbank",
    "organism_match",
    "length_match",
]:
    print(
        f"{col}: {result[col].sum()} / {len(result)}"
    )

print("")
print("Output:", OUT)
print("Summary:", SUMMARY)

if fails:
    print("")
    print("FAILED RECORDS:")
    print(
        failed[
            [
                "kju_id",
                "accession",
                "problems"
            ]
        ].to_string(index=False)
    )

    print("")
    print("VALIDATION DID NOT PASS COMPLETELY.")
    print("DO NOT FREEZE THE 109 DATASET YET.")
    sys.exit(1)

print("")
print("VALIDATION PASSED FOR ALL 109 RECORDS.")
print("DO NOT MODIFY THE ORIGINAL MASTER.")
