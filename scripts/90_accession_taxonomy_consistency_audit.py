#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 90
# CANONICAL SEQUENCE ACCESSION/TAXONOMY CONSISTENCY AUDIT
# ============================================================


BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP83 = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"
STEP87 = ROOT / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
OUT = ROOT / "STEP90_ACCESSION_TAXONOMY_AUDIT"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP83_CSV = STEP83 / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
STEP87_AUDIT = STEP87 / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"

FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"


def normalize_accession(value):
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"", "nan", "none", "null", "na", "n/a"}:
        return ""
    value = value.upper()
    value = re.sub(r"\s+", "", value)
    return value


def accession_root(value):
    v = normalize_accession(value)
    if not v:
        return ""
    return re.sub(r"\.\d+$", "", v)


print("=" * 110)
print("STEP 90: ACCESSION/TAXONOMY CONSISTENCY AUDIT")
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

acc_rows = []
tax_rows = []
disc_rows = []

for _, crow in auto.iterrows():
    kju = crow["KJU_ID"]
    canon_acc = str(crow["Canonical_accession"]).strip()
    canon_acc_root = accession_root(canon_acc)
    accepted_acc = str(crow["Accepted_accession"]).strip()
    accepted_acc_root = accession_root(accepted_acc)

    master_row = master[master["KJU_ID"] == kju].iloc[0]
    accepted_tax = str(master_row.get("Accepted_taxonomy", "")).strip()
    bio_id = str(master_row.get("Biological_identity", "")).strip()

    acc_match = canon_acc_root == accepted_acc_root

    acc_rows.append({
        "KJU_ID": kju,
        "Canonical_accession": canon_acc,
        "Canonical_accession_root": canon_acc_root,
        "Accepted_accession": accepted_acc,
        "Accepted_accession_root": accepted_acc_root,
        "Accession_root_match": str(acc_match),
        "Accepted_taxonomy": accepted_tax,
        "Biological_identity": bio_id,
        "Selection_status": str(crow["Selection_status"]),
    })

    kju_true = true[true["KJU_ID"] == kju]
    source_accs = sorted(set(
        accession_root(x)
        for x in kju_true["Sequence_accession_detected"]
        if accession_root(x)
    ))

    tax_from_source = True
    source_tax_evidence = ""

    if kju_true.empty:
        tax_from_source = False
        source_tax_evidence = "NO_SOURCE_RECORDS"
    elif not acc_match:
        source_tax_evidence = f"ACCESSION_MISMATCH: source={source_accs}, canon={canon_acc_root}"
        tax_from_source = False

    tax_rows.append({
        "KJU_ID": kju,
        "Accepted_taxonomy": accepted_tax,
        "Canonical_accession": canon_acc,
        "Source_accessions": ";".join(source_accs),
        "Taxonomy_from_source": str(tax_from_source),
        "Taxonomy_evidence": source_tax_evidence,
    })

    if not acc_match:
        disc_rows.append({
            "KJU_ID": kju,
            "Discrepancy_type": "ACCESSION_ROOT_MISMATCH",
            "Canonical_accession": canon_acc,
            "Canonical_accession_root": canon_acc_root,
            "Accepted_accession": accepted_acc,
            "Accepted_accession_root": accepted_acc_root,
            "Accepted_taxonomy": accepted_tax,
        })

    if not tax_from_source and source_tax_evidence:
        disc_rows.append({
            "KJU_ID": kju,
            "Discrepancy_type": "TAXONOMY_SOURCE_UNVERIFIED",
            "Canonical_accession": canon_acc,
            "Canonical_accession_root": canon_acc_root,
            "Accepted_accession": accepted_acc,
            "Accepted_accession_root": accepted_acc_root,
            "Accepted_taxonomy": accepted_tax,
        })

unresolved = audit87[audit87["Selection_status"] == "UNRESOLVED"]
for _, urow in unresolved.iterrows():
    kju = urow["KJU_ID"]
    acc = str(urow.get("Accepted_accession", "")).strip()
    tax = str(urow.get("Accepted_taxonomy", "")).strip()
    acc_rows.append({
        "KJU_ID": kju,
        "Canonical_accession": "",
        "Canonical_accession_root": "",
        "Accepted_accession": acc,
        "Accepted_accession_root": accession_root(acc),
        "Accession_root_match": "N/A_NO_SEQUENCE",
        "Accepted_taxonomy": tax,
        "Biological_identity": str(urow.get("Biological_identity", "")),
        "Selection_status": "UNRESOLVED",
    })
    tax_rows.append({
        "KJU_ID": kju,
        "Accepted_taxonomy": tax,
        "Canonical_accession": "",
        "Source_accessions": "",
        "Taxonomy_from_source": "N/A_UNRESOLVED",
        "Taxonomy_evidence": "UNRESOLVED_NO_SEQUENCE",
    })

no_seq = audit87[audit87["Selection_status"] == "NO_SEQUENCE"]
for _, nsrow in no_seq.iterrows():
    kju = nsrow["KJU_ID"]
    acc = str(nsrow.get("Accepted_accession", "")).strip()
    tax = str(nsrow.get("Accepted_taxonomy", "")).strip()
    acc_rows.append({
        "KJU_ID": kju,
        "Canonical_accession": "",
        "Canonical_accession_root": "",
        "Accepted_accession": acc,
        "Accepted_accession_root": accession_root(acc),
        "Accession_root_match": "N/A_NO_SEQUENCE",
        "Accepted_taxonomy": tax,
        "Biological_identity": str(nsrow.get("Biological_identity", "")),
        "Selection_status": "NO_SEQUENCE",
    })
    tax_rows.append({
        "KJU_ID": kju,
        "Accepted_taxonomy": tax,
        "Canonical_accession": "",
        "Source_accessions": "",
        "Taxonomy_from_source": "N/A_NO_SEQUENCE",
        "Taxonomy_evidence": f"NO_SEQUENCE: accession={acc}",
    })

acc_df = pd.DataFrame(acc_rows)
tax_df = pd.DataFrame(tax_rows)
disc_df = pd.DataFrame(disc_rows)

if disc_df.empty:
    disc_df = pd.DataFrame(columns=["KJU_ID", "Discrepancy_type", "Canonical_accession",
                                     "Canonical_accession_root", "Accepted_accession",
                                     "Accepted_accession_root", "Accepted_taxonomy"])

acc_df.to_csv(OUT / "STEP90_ACCESSION_CONSISTENCY.csv", index=False)
tax_df.to_csv(OUT / "STEP90_TAXONOMY_CONSISTENCY.csv", index=False)
disc_df.to_csv(OUT / "STEP90_DISCREPANCIES.csv", index=False)

n_auto = len(auto)
acc_match_count = int(acc_df[acc_df["Accession_root_match"] == "True"].shape[0])
acc_mismatch_count = int(disc_df[disc_df["Discrepancy_type"] == "ACCESSION_ROOT_MISMATCH"].shape[0]) if "Discrepancy_type" in disc_df.columns else 0
tax_unverified = int(disc_df[disc_df["Discrepancy_type"] == "TAXONOMY_SOURCE_UNVERIFIED"].shape[0]) if "Discrepancy_type" in disc_df.columns else 0

print(f"\n=== ACCESSION CONSISTENCY ===")
print(f"Auto-selected canonical: {n_auto}")
print(f"Accession root match:    {acc_match_count}")
print(f"Accession root mismatch: {acc_mismatch_count}")
print(f"Taxonomy unverified:     {tax_unverified}")
print(f"Total discrepancies:     {len(disc_df)}")

if len(disc_df) > 0:
    print(f"\nDiscrepancies:")
    for _, d in disc_df.iterrows():
        print(f"  {d['KJU_ID']}: {d['Discrepancy_type']}")

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(OUT / "STEP90_AUDIT_REPORT.txt", "w") as f:
    f.write("STEP 90: ACCESSION/TAXONOMY CONSISTENCY AUDIT\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS ({master_hash})\n\n")
    f.write(f"Auto-selected canonical: {n_auto}\n")
    f.write(f"Accession root match:    {acc_match_count}\n")
    f.write(f"Accession root mismatch: {acc_mismatch_count}\n")
    f.write(f"Taxonomy unverified:     {tax_unverified}\n")
    f.write(f"Total discrepancies:     {len(disc_df)}\n\n")
    f.write(f"Frozen V1 master untouched.\n")

print("\n" + "=" * 110)
print("STEP 90 COMPLETE")
print("=" * 110)
