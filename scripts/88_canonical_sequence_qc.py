#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
import os


# ============================================================
# STEP 88
# CANONICAL SEQUENCE QC + PROVENANCE VALIDATION
#
# PURPOSE:
#   Deep validation of every auto-selected canonical
#   sequence from STEP87 against STEP83 source data
#   and the frozen master.
# ============================================================


BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "final_integrated_master"
)

FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP83 = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"
STEP87 = ROOT / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
OUT = ROOT / "STEP88_CANONICAL_SEQUENCE_QC"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP83_CSV = STEP83 / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
STEP87_AUDIT = STEP87 / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"

QC_DETAIL = OUT / "STEP88_CANONICAL_QC_DETAIL.csv"
PROVENANCE = OUT / "STEP88_PROVENANCE_AUDIT.csv"
REPORT = OUT / "STEP88_QC_REPORT.txt"

FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

VALID_IUPAC = set("ACGTNRYWSMKBDHVU")


def normalize_accession(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"", "nan", "none", "null", "na", "n/a"}:
        return ""
    value = value.upper()
    value = re.sub(r"\s+", "", value)
    value = re.sub(r"\.\d+$", "", value)
    return value


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
print("STEP 88: CANONICAL SEQUENCE QC + PROVENANCE VALIDATION")
print("=" * 110)

for p in [MASTER_CSV, STEP83_CSV, STEP87_CANONICAL, STEP87_AUDIT]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
print(f"\nFrozen master SHA256: {master_hash}")
print(f"Expected SHA256:      {FROZEN_SHA}")
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print("Frozen master integrity: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
true = pd.read_csv(STEP83_CSV, dtype=str, low_memory=False)
canonical = pd.read_csv(STEP87_CANONICAL, dtype=str, low_memory=False)
audit87 = pd.read_csv(STEP87_AUDIT, dtype=str, low_memory=False)

print(f"\nMaster rows:   {len(master)}")
print(f"STEP83 rows:   {len(true)}")
print(f"Canonical rows: {len(canonical)}")
print(f"Audit87 rows:  {len(audit87)}")

assert len(audit87) == 204
assert audit87["KJU_ID"].nunique() == 204

auto_selected = canonical[
    canonical["Selection_status"].str.startswith("AUTO_SELECTED_", na=False)
].copy()

print(f"Auto-selected canonical: {len(auto_selected)}")

qc_rows = []
prov_rows = []

for _, row in auto_selected.iterrows():
    kju = row["KJU_ID"]
    canon_sha = str(row["Canonical_sequence_SHA256"]).strip()
    canon_seq = str(row["Canonical_sequence"]).strip()
    canon_acc = str(row["Canonical_accession"]).strip()
    canon_len_raw = str(row["Canonical_sequence_length"]).strip()
    canon_len = numeric_length(canon_len_raw)
    accepted_acc = str(row["Accepted_accession"]).strip()
    accepted_len = numeric_length(row["Accepted_length"])
    selection_status = str(row["Selection_status"]).strip()
    selection_evidence = str(row["Selection_evidence"]).strip()

    checks = {}

    import hashlib as hl
    calc_sha = hl.sha256(canon_seq.encode("utf-8")).hexdigest()
    checks["sha256_match"] = calc_sha == canon_sha

    checks["sequence_valid_iupac"] = all(c in VALID_IUPAC for c in canon_seq)

    checks["sequence_not_empty"] = len(canon_seq) > 0

    checks["sha256_not_empty"] = len(canon_sha) > 0

    checks["accession_not_empty"] = len(canon_acc) > 0

    checks["length_not_empty"] = canon_len is not None

    canon_acc_root = normalize_accession(canon_acc)
    accepted_acc_root = normalize_accession(accepted_acc)
    checks["accession_root_match"] = canon_acc_root == accepted_acc_root

    if accepted_len is not None and canon_len is not None:
        checks["length_match"] = canon_len == accepted_len
    else:
        checks["length_match"] = True

    source_match = true[
        (true["KJU_ID"] == kju)
        & (true["Sequence_SHA256"] == canon_sha)
    ]
    checks["exists_in_step83"] = len(source_match) > 0

    if checks["exists_in_step83"]:
        seq_match = source_match[
            source_match["Sequence"].apply(
                lambda x: str(x).strip() == canon_seq
            )
        ]
        checks["sequence_matches_source"] = len(seq_match) > 0
    else:
        checks["sequence_matches_source"] = False

    source_files = set()
    source_record_ids = set()
    if checks["exists_in_step83"]:
        for _, sr in source_match.iterrows():
            sf = str(sr.get("Sequence_file", "")).strip()
            if sf and sf.lower() not in {"nan", "none"}:
                source_files.add(sf)
                fname = sf.split("/")[-1]
                source_files.add(fname)
            rid = str(sr.get("Record_ID", "")).strip()
            if rid and rid.lower() not in {"nan", "none"}:
                source_record_ids.add(rid)

    source_files_exist = True
    for sf in source_files:
        if sf.startswith("/") and not os.path.exists(sf):
            source_files_exist = False

    checks["source_files_exist"] = source_files_exist

    all_pass = all(checks.values())
    failed = [k for k, v in checks.items() if not v]

    qc_rows.append({
        "KJU_ID": kju,
        "Selection_status": selection_status,
        "Selection_evidence": selection_evidence,
        "Accepted_accession": accepted_acc,
        "Accepted_length": str(accepted_len) if accepted_len else "",
        "Canonical_accession": canon_acc,
        "Canonical_accession_root": canon_acc_root,
        "Accepted_accession_root": accepted_acc_root,
        "Canonical_sequence_length": str(canon_len) if canon_len else "",
        "Canonical_sequence_SHA256": canon_sha,
        "Calculated_SHA256": calc_sha,
        "sha256_match": str(checks["sha256_match"]),
        "sequence_valid_iupac": str(checks["sequence_valid_iupac"]),
        "sequence_not_empty": str(checks["sequence_not_empty"]),
        "accession_root_match": str(checks["accession_root_match"]),
        "length_match": str(checks["length_match"]),
        "exists_in_step83": str(checks["exists_in_step83"]),
        "sequence_matches_source": str(checks["sequence_matches_source"]),
        "source_files_exist": str(checks["source_files_exist"]),
        "all_checks_pass": str(all_pass),
        "failed_checks": ";".join(failed) if failed else "",
    })

    prov_rows.append({
        "KJU_ID": kju,
        "Canonical_accession": canon_acc,
        "Canonical_sequence_SHA256": canon_sha,
        "Canonical_sequence_length": str(canon_len) if canon_len else "",
        "Source_record_count": len(source_match),
        "Source_record_IDs": ";".join(sorted(source_record_ids)),
        "Source_files": ";".join(sorted(source_files)),
        "Accepted_accession": accepted_acc,
        "Accepted_length": str(accepted_len) if accepted_len else "",
        "Accepted_taxonomy": str(row.get("Accepted_taxonomy", "")),
        "Biological_identity": str(row.get("Biological_identity", "")),
    })

qc_df = pd.DataFrame(qc_rows)
qc_df.to_csv(QC_DETAIL, index=False)

prov_df = pd.DataFrame(prov_rows)
prov_df.to_csv(PROVENANCE, index=False)

print(f"\nSaved QC detail: {QC_DETAIL}")
print(f"Saved provenance: {PROVENANCE}")


all_pass_count = (qc_df["all_checks_pass"] == "True").sum()
any_fail = (qc_df["all_checks_pass"] == "False").sum()
sha_mismatches = (qc_df["sha256_match"] == "False").sum()
invalid_seq = (qc_df["sequence_valid_iupac"] == "False").sum()
acc_mismatches = (qc_df["accession_root_match"] == "False").sum()
len_mismatches = (qc_df["length_match"] == "False").sum()
missing_step83 = (qc_df["exists_in_step83"] == "False").sum()
missing_source = (qc_df["source_files_exist"] == "False").sum()

print(f"\n=== QC SUMMARY ===")
print(f"Total canonical: {len(auto_selected)}")
print(f"All checks pass: {all_pass_count}")
print(f"Any fail:        {any_fail}")
print(f"SHA256 mismatch: {sha_mismatches}")
print(f"Invalid IUPAC:   {invalid_seq}")
print(f"Accession mismatch: {acc_mismatches}")
print(f"Length mismatch: {len_mismatches}")
print(f"Missing from STEP83: {missing_step83}")
print(f"Source files missing: {missing_source}")

if any_fail > 0:
    failed_rows = qc_df[qc_df["all_checks_pass"] == "False"]
    print(f"\nFAILED KJUs:")
    for _, fr in failed_rows.iterrows():
        print(f"  {fr['KJU_ID']}: {fr['failed_checks']}")


no_seq_audit = audit87[audit87["Selection_status"] == "NO_SEQUENCE"]
unresolved_audit = audit87[audit87["Selection_status"] == "UNRESOLVED"]

print(f"\n=== NON-CANONICAL CASES ===")
print(f"NO_SEQUENCE:   {len(no_seq_audit)}")
print(f"UNRESOLVED:    {len(unresolved_audit)}")

for _, nr in no_seq_audit.iterrows():
    print(f"  {nr['KJU_ID']}: accession={nr['Accepted_accession']}, taxonomy={nr['Accepted_taxonomy']}")

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA, "FROZEN MASTER CHANGED"

with open(REPORT, "w") as f:
    f.write("STEP 88: CANONICAL SEQUENCE QC + PROVENANCE VALIDATION\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: {master_hash}\n")
    f.write(f"Expected SHA256:      {FROZEN_SHA}\n")
    f.write(f"Status:               PASS\n\n")
    f.write(f"Total canonical sequences: {len(auto_selected)}\n")
    f.write(f"All checks pass:          {all_pass_count}\n")
    f.write(f"Any fail:                 {any_fail}\n\n")
    f.write(f"SHA256 mismatch:          {sha_mismatches}\n")
    f.write(f"Invalid IUPAC:            {invalid_seq}\n")
    f.write(f"Accession root mismatch:  {acc_mismatches}\n")
    f.write(f"Length mismatch:          {len_mismatches}\n")
    f.write(f"Missing from STEP83:      {missing_step83}\n")
    f.write(f"Source files missing:     {missing_source}\n\n")
    f.write(f"NO_SEQUENCE KJUs:         {len(no_seq_audit)}\n")
    f.write(f"UNRESOLVED KJUs:          {len(unresolved_audit)}\n\n")
    f.write(f"Frozen V1 master untouched.\n")

print(f"\nSaved report: {REPORT}")

print("\n" + "=" * 110)
print("STEP 88 COMPLETE")
print("=" * 110)
