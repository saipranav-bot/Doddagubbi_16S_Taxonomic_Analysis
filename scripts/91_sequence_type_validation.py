#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 91
# 16S/ITS SEQUENCE-TYPE VALIDATION
# ============================================================


BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP83 = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"
STEP87 = ROOT / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
OUT = ROOT / "STEP91_SEQUENCE_TYPE_VALIDATION"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP83_CSV = STEP83 / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
STEP87_AUDIT = STEP87 / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"

FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

SIXS_MIN = 1200
SIXS_MAX = 1550
ITS_MIN = 200
ITS_MAX = 900


def numeric_length(value):
    if pd.isna(value):
        return None
    value = str(value).strip()
    if not value or value.lower() in {"nan", "none", "null", "na", "n/a"}:
        return None
    try:
        n = float(value)
        if n != n or n <= 0:
            return None
        return int(n)
    except Exception:
        return None


print("=" * 110)
print("STEP 91: 16S/ITS SEQUENCE-TYPE VALIDATION")
print("=" * 110)

for p in [MASTER_CSV, STEP83_CSV, STEP87_CANONICAL, STEP87_AUDIT]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master SHA256: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
true = pd.read_csv(STEP83_CSV, dtype=str, low_memory=False)
canonical = pd.read_csv(STEP87_CANONICAL, dtype=str, low_memory=False)
audit87 = pd.read_csv(STEP87_AUDIT, dtype=str, low_memory=False)

auto = canonical[canonical["Selection_status"].str.startswith("AUTO_SELECTED_", na=False)].copy()

print(f"Master: {len(master)}, Canonical: {len(auto)}, STEP83: {len(true)}")

type_rows = []
anomaly_rows = []
length_dist_rows = []

type_counts = {}
length_counts = {}

for _, crow in auto.iterrows():
    kju = crow["KJU_ID"]
    seq_len = numeric_length(crow["Canonical_sequence_length"])
    accepted_len = numeric_length(crow["Accepted_length"])

    kju_true = true[true["KJU_ID"] == kju]
    seq_types_in_source = sorted(set(
        str(x).strip()
        for x in kju_true["Sequence_type"]
        if str(x).strip() and str(x).strip().lower() not in {"nan", "none"}
    ))
    primary_type = seq_types_in_source[0] if seq_types_in_source else "UNKNOWN"

    type_counts[primary_type] = type_counts.get(primary_type, 0) + 1

    if seq_len is not None:
        bucket = f"{(seq_len // 50) * 50}-{(seq_len // 50) * 50 + 49}"
        length_counts[bucket] = length_counts.get(bucket, 0) + 1

    in_range = True
    range_detail = ""
    if primary_type == "16S":
        in_range = SIXS_MIN <= seq_len <= SIXS_MAX if seq_len else False
        range_detail = f"16S expected {SIXS_MIN}-{SIXS_MAX}"
    elif primary_type == "ITS":
        in_range = ITS_MIN <= seq_len <= ITS_MAX if seq_len else False
        range_detail = f"ITS expected {ITS_MIN}-{ITS_MAX}"

    type_rows.append({
        "KJU_ID": kju,
        "Canonical_sequence_type": primary_type,
        "All_source_types": ";".join(seq_types_in_source),
        "Canonical_sequence_length": str(seq_len) if seq_len else "",
        "Accepted_length": str(accepted_len) if accepted_len else "",
        "Length_in_expected_range": str(in_range),
        "Expected_range": range_detail,
        "Canonical_accession": str(crow["Canonical_accession"]),
    })

    if not in_range:
        anomaly_rows.append({
            "KJU_ID": kju,
            "Anomaly_type": "LENGTH_OUT_OF_RANGE",
            "Sequence_type": primary_type,
            "Sequence_length": str(seq_len),
            "Expected_range": range_detail,
            "Canonical_accession": str(crow["Canonical_accession"]),
        })

    if accepted_len is not None and seq_len is not None and seq_len != accepted_len:
        anomaly_rows.append({
            "KJU_ID": kju,
            "Anomaly_type": "CANONICAL_LENGTH_DIFFERS_FROM_ACCEPTED",
            "Sequence_type": primary_type,
            "Canonical_length": str(seq_len),
            "Accepted_length": str(accepted_len),
            "Expected_range": "",
            "Canonical_accession": str(crow["Canonical_accession"]),
        })

for kju_count, count in sorted(type_counts.items()):
    length_dist_rows.append({
        "Category": f"Type_{kju_count}",
        "Count": str(count),
    })

for bucket, count in sorted(length_counts.items()):
    length_dist_rows.append({
        "Category": f"Length_{bucket}",
        "Count": str(count),
    })

type_df = pd.DataFrame(type_rows)
anomaly_df = pd.DataFrame(anomaly_rows)
dist_df = pd.DataFrame(length_dist_rows)

type_df.to_csv(OUT / "STEP91_SEQUENCE_TYPE_AUDIT.csv", index=False)
anomaly_df.to_csv(OUT / "STEP91_ANOMALIES.csv", index=False)
dist_df.to_csv(OUT / "STEP91_LENGTH_DISTRIBUTION.csv", index=False)

print(f"\n=== SEQUENCE TYPE DISTRIBUTION ===")
for t, c in sorted(type_counts.items()):
    print(f"  {t}: {c}")

print(f"\n=== LENGTH RANGE ===")
all_lens = [numeric_length(r["Canonical_sequence_length"]) for _, r in auto.iterrows()]
all_lens = [l for l in all_lens if l is not None]
print(f"  Min: {min(all_lens)}")
print(f"  Max: {max(all_lens)}")
print(f"  Mean: {sum(all_lens)/len(all_lens):.0f}")

print(f"\n=== ANOMALIES ===")
print(f"  Total anomalies: {len(anomaly_df)}")
for _, a in anomaly_df.iterrows():
    print(f"  {a['KJU_ID']}: {a['Anomaly_type']}")

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(OUT / "STEP91_REPORT.txt", "w") as f:
    f.write("STEP 91: 16S/ITS SEQUENCE-TYPE VALIDATION\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS\n\n")
    f.write(f"Total canonical: {len(auto)}\n\n")
    f.write("SEQUENCE TYPE DISTRIBUTION:\n")
    for t, c in sorted(type_counts.items()):
        f.write(f"  {t}: {c}\n")
    f.write(f"\nLENGTH RANGE: {min(all_lens)}-{max(all_lens)}\n")
    f.write(f"MEAN: {sum(all_lens)/len(all_lens):.0f}\n\n")
    f.write(f"ANOMALIES: {len(anomaly_df)}\n")
    for _, a in anomaly_df.iterrows():
        f.write(f"  {a['KJU_ID']}: {a['Anomaly_type']}\n")
    f.write(f"\nFrozen V1 master untouched.\n")

print("\n" + "=" * 110)
print("STEP 91 COMPLETE")
print("=" * 110)
