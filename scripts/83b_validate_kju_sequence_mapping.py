#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 83B
# KJU SEQUENCE MAPPING VALIDATION
#
# PURPOSE:
#   Validate every STEP83 sequence KJU against the frozen
#   204-KJU master. Preserve valid KJUs, isolate invalid
#   KJU mappings, retain original record provenance.
#
# DO NOT modify STEP83 outputs.
# DO NOT fabricate mappings.
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

OUT = (
    ROOT
    / "STEP83B_KJU_MAPPING_VALIDATION"
)

OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = (
    FREEZE
    / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
)

STEP83_CSV = (
    STEP83
    / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
)


VALID_CSV = OUT / "STEP83B_VALID_SEQUENCE_RECORDS.csv"
INVALID_CSV = OUT / "STEP83B_INVALID_KJU_RECORDS.csv"
AUDIT_CSV = OUT / "STEP83B_KJU_MAPPING_AUDIT.csv"
REPORT_TXT = OUT / "STEP83B_REPORT.txt"


FROZEN_SHA = (
    "7d858153a0aaf19e3e1f4f7a1d3b03b3"
    "acc69ec6451c6c71541213172196c8a8"
)


def normalize_kju(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"", "nan", "none", "null", "na", "n/a"}:
        return ""
    value = re.sub(r"\s+", "", value.upper())
    if re.fullmatch(r"KJU\d{4}", value):
        return value
    m = re.fullmatch(r"KJU(\d{1,3})", value)
    if m:
        return f"KJU{int(m.group(1)):04d}"
    return value


print("=" * 110)
print("STEP 83B: KJU SEQUENCE MAPPING VALIDATION")
print("=" * 110)

for p in [MASTER_CSV, STEP83_CSV]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")


master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
print(f"\nFrozen master SHA256: {master_hash}")
print(f"Expected SHA256:      {FROZEN_SHA}")
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print("Frozen master integrity: PASS")


master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
true = pd.read_csv(STEP83_CSV, dtype=str, low_memory=False)

print(f"\nMaster rows: {len(master)}")
print(f"STEP83 rows: {len(true)}")
print(f"STEP83 columns: {list(true.columns)}")

master["_KJU_NORM"] = master["KJU_ID"].map(normalize_kju)
master_ids = set(x for x in master["_KJU_NORM"] if x)

assert len(master_ids) == 204, f"Expected 204 master IDs, got {len(master_ids)}"
print(f"Master unique KJU IDs: {len(master_ids)}")


true["_KJU_NORM"] = true["KJU_ID"].map(normalize_kju)

is_blank = true["_KJU_NORM"] == ""
is_valid = (~is_blank) & (true["_KJU_NORM"].isin(master_ids))
is_invalid = (~is_blank) & (~true["_KJU_NORM"].isin(master_ids))

valid_rows = true[is_valid].copy()
invalid_rows = pd.concat(
    [true[is_blank], true[is_invalid]],
    ignore_index=True,
)

print(f"\n=== MAPPING CLASSIFICATION ===")
print(f"Valid KJU (in master)  : {len(valid_rows)}")
print(f"Blank/NaN KJU          : {int(is_blank.sum())}")
print(f"Invalid KJU (not in master): {int(is_invalid.sum())}")

unique_valid_kjus = valid_rows["_KJU_NORM"].nunique()
print(f"Unique valid KJU IDs   : {unique_valid_kjus}")

if int(is_invalid.sum()) > 0:
    bad_ids = sorted(true[is_invalid]["_KJU_NORM"].unique())
    print(f"\nINVALID KJU IDs ({len(bad_ids)}):")
    for bid in bad_ids:
        raw = true[true["_KJU_NORM"] == bid]["KJU_ID"].unique()
        cnt = int((true[true["_KJU_NORM"] == bid]).shape[0])
        print(f"  {repr(bid)}  raw={raw}  count={cnt}")
    raise ValueError("Invalid KJU IDs found — must be zero")


valid_rows.to_csv(VALID_CSV, index=False)
invalid_rows.to_csv(INVALID_CSV, index=False)

print(f"\nSaved valid records:   {VALID_CSV}")
print(f"Saved invalid records: {INVALID_CSV}")


audit_rows = []

for kju in sorted(master_ids):
    kju_rows = valid_rows[valid_rows["_KJU_NORM"] == kju]
    n_rows = len(kju_rows)

    if n_rows == 0:
        audit_rows.append({
            "KJU_ID": kju,
            "Mapping_status": "NO_SEQUENCE_RECORDS",
            "Raw_record_count": 0,
            "Unique_accessions": "",
            "Unique_sequences": 0,
            "Unique_sequence_files": 0,
            "Unique_record_IDs": "",
            "Sequence_types": "",
            "Source_files": "",
        })
        continue

    accessions = sorted(set(
        str(x).strip()
        for x in kju_rows["Sequence_accession_detected"]
        if str(x).strip() and str(x).strip().lower() not in {"nan", "none"}
    ))
    n_seq = int(kju_rows["Sequence_SHA256"].nunique())
    files = sorted(set(
        str(x).strip().split("/")[-1]
        for x in kju_rows["Sequence_file"]
        if str(x).strip() and str(x).strip().lower() not in {"nan", "none"}
    ))
    record_ids = sorted(set(
        str(x).strip()
        for x in kju_rows["Record_ID"]
        if str(x).strip() and str(x).strip().lower() not in {"nan", "none"}
    ))
    seq_types = sorted(set(
        str(x).strip()
        for x in kju_rows["Sequence_type"]
        if str(x).strip() and str(x).strip().lower() not in {"nan", "none"}
    ))

    audit_rows.append({
        "KJU_ID": kju,
        "Mapping_status": "VALID",
        "Raw_record_count": n_rows,
        "Unique_accessions": ";".join(accessions),
        "Unique_sequences": n_seq,
        "Unique_sequence_files": len(files),
        "Unique_record_IDs": ";".join(record_ids),
        "Sequence_types": ";".join(seq_types),
        "Source_files": ";".join(files),
    })

audit = pd.DataFrame(audit_rows)
audit.to_csv(AUDIT_CSV, index=False)

print(f"\nSaved mapping audit: {AUDIT_CSV}")


with open(REPORT_TXT, "w") as f:
    f.write("STEP 83B: KJU SEQUENCE MAPPING VALIDATION\n")
    f.write("=" * 110 + "\n\n")

    f.write(f"Frozen master SHA256: {master_hash}\n")
    f.write(f"Expected SHA256:      {FROZEN_SHA}\n")
    f.write(f"Status:               PASS\n\n")

    f.write(f"Master rows:          {len(master)}\n")
    f.write(f"Master unique KJUs:   {len(master_ids)}\n")
    f.write(f"STEP83 raw rows:      {len(true)}\n\n")

    f.write(f"=== MAPPING CLASSIFICATION ===\n")
    f.write(f"Valid (in master):    {len(valid_rows)}\n")
    f.write(f"Blank/NaN KJU:        {int(is_blank.sum())}\n")
    f.write(f"Invalid (not in master): {int(is_invalid.sum())}\n")
    f.write(f"Unique valid KJUs:    {unique_valid_kjus}\n\n")

    f.write(f"=== PER-KJU AUDIT ===\n")
    f.write(f"KJUs with sequence:   {len(audit[audit['Mapping_status'] == 'VALID'])}\n")
    f.write(f"KJUs without sequence:{len(audit[audit['Mapping_status'] == 'NO_SEQUENCE_RECORDS'])}\n\n")

    multi_seq = audit[(audit["Mapping_status"] == "VALID") & (audit["Unique_sequences"] > 1)]
    single_seq = audit[(audit["Mapping_status"] == "VALID") & (audit["Unique_sequences"] == 1)]
    f.write(f"Single unique sequence:  {len(single_seq)}\n")
    f.write(f"Multiple unique seqs:    {len(multi_seq)}\n\n")

    f.write(f"=== SAFETY ===\n")
    f.write(f"No STEP83 outputs modified.\n")
    f.write(f"No mappings fabricated.\n")
    f.write(f"Frozen V1 master untouched.\n")

print(f"Saved report: {REPORT_TXT}")


print("\n" + "=" * 110)
print("STEP 83B COMPLETE")
print("=" * 110)
print(f"Valid records:   {len(valid_rows)}")
print(f"Invalid records: {len(invalid_rows)}")
print(f"Unique valid KJUs: {unique_valid_kjus}")
print(f"Frozen master SHA256 verified: PASS")
