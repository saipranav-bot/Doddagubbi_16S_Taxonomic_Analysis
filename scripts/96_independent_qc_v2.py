#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re

BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP95 = ROOT / "STEP95_FINAL_INTEGRATED_MASTER_V2"
STEP83 = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"
OUT = ROOT / "STEP96_INDEPENDENT_QC_V2"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
V2_CSV = STEP95 / "STEP95_INTEGRATED_MASTER_V2.csv"
STEP83_CSV = STEP83 / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

VALID_IUPAC = set("ACGTNRYWSMKBDHVU")

print("=" * 110)
print("STEP 96: INDEPENDENT QC OF V2")
print("=" * 110)

for p in [MASTER_CSV, V2_CSV, STEP83_CSV]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master SHA256: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
v2 = pd.read_csv(V2_CSV, dtype=str, low_memory=False)
true = pd.read_csv(STEP83_CSV, dtype=str, low_memory=False)

qc_results = []

# QC1: Row count
qc_results.append(("V2 row count == 204", len(v2) == 204, f"got {len(v2)}"))

# QC2: Unique KJU count
qc_results.append(("V2 unique KJUs == 204", v2["KJU_ID"].nunique() == 204, f"got {v2['KJU_ID'].nunique()}"))

# QC3: KJU IDs match master
master_ids = set(master["KJU_ID"])
v2_ids = set(v2["KJU_ID"])
qc_results.append(("V2 KJU IDs == master KJU IDs", v2_ids == master_ids, f"diff={v2_ids.symmetric_difference(master_ids)}"))

# QC4: Status counts sum to 204
status_counts = v2["Final_status"].value_counts().to_dict()
total_status = sum(status_counts.values())
qc_results.append(("Status counts sum to 204", total_status == 204, f"got {total_status}"))

canon_count = sum(v for k, v in status_counts.items() if "CANONICAL" in k)
noseq_count = status_counts.get("ADJUDICATED_NO_SEQUENCE", 0)
unresolved_count = status_counts.get("UNRESOLVED", 0)
qc_results.append(("canonical + noseq + unresolved == 204", canon_count + noseq_count + unresolved_count == 204, f"{canon_count}+{noseq_count}+{unresolved_count}"))

# QC5: Frozen master SHA256
master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
qc_results.append(("Frozen master SHA256 unchanged", master_hash_after == FROZEN_SHA, f"got {master_hash_after[:16]}..."))

# QC6: All canonical sequences have SHA
canon_rows = v2[v2["Final_status"].str.contains("CANONICAL", na=False)]
blank_sha = canon_rows[canon_rows["Canonical_sequence_SHA256"].apply(lambda x: str(x).strip() in {"", "nan", "none"})]
qc_results.append(("All canonical have SHA256", len(blank_sha) == 0, f"{len(blank_sha)} blank"))

# QC7: All canonical sequences have sequence
blank_seq = canon_rows[canon_rows["Canonical_sequence"].apply(lambda x: str(x).strip() in {"", "nan", "none"})]
qc_results.append(("All canonical have sequence", len(blank_seq) == 0, f"{len(blank_seq)} blank"))

# QC8: SHA256 integrity
sha_fail = 0
for _, row in canon_rows.iterrows():
    seq = str(row["Canonical_sequence"]).strip().upper()
    sha = str(row["Canonical_sequence_SHA256"]).strip()
    calc = hashlib.sha256(seq.encode("utf-8")).hexdigest()
    if calc != sha:
        sha_fail += 1
qc_results.append(("SHA256 integrity for all canonical", sha_fail == 0, f"{sha_fail} mismatches"))

# QC9: IUPAC validity
iupac_fail = 0
for _, row in canon_rows.iterrows():
    seq = str(row["Canonical_sequence"]).strip().upper()
    if not all(c in VALID_IUPAC for c in seq):
        iupac_fail += 1
qc_results.append(("IUPAC validity for all canonical", iupac_fail == 0, f"{iupac_fail} invalid"))

# QC10: No unresolved KJU has sequence
unresolved_rows = v2[v2["Final_status"] == "UNRESOLVED"]
unresolved_with_seq = unresolved_rows[unresolved_rows["Canonical_sequence"].apply(lambda x: str(x).strip() not in {"", "nan", "none"})]
qc_results.append(("No unresolved has sequence", len(unresolved_with_seq) == 0, f"{len(unresolved_with_seq)}"))

# QC11: No fabricated accession
no_fabricated_acc = True
for _, row in v2.iterrows():
    canon_acc = str(row.get("Canonical_accession", "")).strip()
    accepted = str(row.get("Accepted_accession", "")).strip()
    if canon_acc and canon_acc.lower() not in {"nan", "none"}:
        canon_root = re.sub(r"\.\d+$", "", canon_acc.upper())
        acc_root = re.sub(r"\.\d+$", "", accepted.upper())
        if canon_root != acc_root:
            no_fabricated_acc = False
qc_results.append(("No fabricated accession", no_fabricated_acc, ""))

# QC12: No fabricated taxonomy
no_fabricated_tax = True
for _, row in v2.iterrows():
    bio = str(row.get("Final_biological_identity", "")).strip()
    accepted = str(row.get("Biological_identity", "")).strip()
    if bio and bio.lower() not in {"", "nan", "none"} and bio != accepted:
        no_fabricated_tax = False
qc_results.append(("No fabricated biological identity", no_fabricated_tax, ""))

# QC13: KJU0039 is NO_SEQUENCE
k0039 = v2[v2["KJU_ID"] == "KJU0039"]
k0039_ok = len(k0039) == 1 and k0039.iloc[0]["Final_status"] == "ADJUDICATED_NO_SEQUENCE"
qc_results.append(("KJU0039 is NO_SEQUENCE", k0039_ok, ""))

# QC14: Special KJUs
for kju, expected_status in [("KJU0038", "ADJUDICATED_CANONICAL"), ("KJU0105", "ADJUDICATED_CANONICAL"),
                               ("KJU0224", "ADJUDICATED_CANONICAL"), ("KJU0225", "ADJUDICATED_CANONICAL")]:
    row = v2[v2["KJU_ID"] == kju]
    ok = len(row) == 1 and row.iloc[0]["Final_status"] == expected_status
    qc_results.append((f"{kju} is {expected_status}", ok, ""))

# Print results
print(f"\n=== QC RESULTS ({sum(1 for _, p, _ in qc_results if p)}/{len(qc_results)} pass) ===")
all_pass = True
for desc, passed, detail in qc_results:
    status = "PASS" if passed else "FAIL"
    if not passed:
        all_pass = False
    print(f"  [{status}] {desc}" + (f" ({detail})" if detail else ""))

qc_df = pd.DataFrame([{"Check": d, "Pass": str(p), "Detail": det} for d, p, det in qc_results])
qc_df.to_csv(OUT / "STEP96_QC_DETAIL.csv", index=False)

with open(OUT / "STEP96_QC_REPORT.txt", "w") as f:
    f.write("STEP 96: INDEPENDENT QC OF V2\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS\n\n")
    f.write(f"Total checks: {len(qc_results)}\n")
    f.write(f"Passed: {sum(1 for _, p, _ in qc_results if p)}\n")
    f.write(f"Failed: {sum(1 for _, p, _ in qc_results if not p)}\n\n")
    for desc, passed, detail in qc_results:
        status = "PASS" if passed else "FAIL"
        f.write(f"  [{status}] {desc}" + (f" ({detail})" if detail else "") + "\n")
    f.write(f"\nFrozen V1 master untouched.\n")

if not all_pass:
    print("\nWARNING: Some QC checks failed")
else:
    print("\nALL QC CHECKS PASSED")

print("\n" + "=" * 110)
print("STEP 96 COMPLETE")
print("=" * 110)
