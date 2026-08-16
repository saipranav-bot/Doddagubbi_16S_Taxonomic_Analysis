#!/usr/bin/env python3
"""
STEP 101: Release/Package Integrity Check

Verifies the Step 100 release package is internally consistent
before any downstream work. Creates a clean analysis workspace.

READ-ONLY on V1, V2, and STEP100_RELEASE_PACKAGE.
"""

import os
import re
import hashlib
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

# ──────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────
HOME = Path.home()
BASE = HOME / "Doddagubbi_Master"
RESULTS = BASE / "github_results/reconciliation_111_vs_109/final_integrated_master"

REL_DIR = RESULTS / "STEP100_RELEASE_PACKAGE"
V2_DIR = RESULTS / "STEP98_FINAL_FREEZE_V2/FINAL_KJU_MASTER_FREEZE_V2"
V1_DIR = RESULTS / "FINAL_KJU_MASTER_FREEZE_V1"

OUT = BASE / "downstream_analysis" / "STEP101_RELEASE_INTEGRITY"
OUT.mkdir(parents=True, exist_ok=True)

# Create all downstream directories
DOWNSTREAM = BASE / "downstream_analysis"
for d in [
    "STEP101_RELEASE_INTEGRITY",
    "STEP102_ANALYSIS_METADATA",
    "STEP103_TAXONOMIC_COMPOSITION",
    "STEP104_SEQUENCE_DIVERSITY",
    "STEP105_REPRODUCIBLE_ALIGNMENT",
    "STEP106_REPRODUCIBLE_PHYLOGENY",
    "STEP107_PHYLOGENETIC_TAXONOMIC_DIVERSITY",
]:
    (DOWNSTREAM / d).mkdir(parents=True, exist_ok=True)

V1_SHA256_EXPECTED = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

# ──────────────────────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────────────────────
def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_seq(seq_str):
    return hashlib.sha256(seq_str.upper().encode("utf-8")).hexdigest()


IUPAC = set("ACGTURYSWKMBDHVN")


def iupac_valid(seq):
    return all(c in IUPAC for c in seq.upper().strip())


print("STEP 101: RELEASE/PACKAGE INTEGRITY CHECK")
print("=" * 70)
print()

# ──────────────────────────────────────────────────────────
# CHECKS
# ──────────────────────────────────────────────────────────
checks = []  # (check_num, name, status, detail)


def check(num, name, condition, detail):
    status = "PASS" if condition else "FAIL"
    checks.append((num, name, status, detail))
    tag = "PASS" if condition else "FAIL"
    print(f"  [{tag}] Check {num:02d} — {name}: {detail}")
    return condition


all_pass = True

# ──────────────────────────────────────────────────────────
# CHECK 1-6: FILE EXISTENCE
# ──────────────────────────────────────────────────────────
print("\nFILE EXISTENCE")
print("-" * 70)

rel_package = REL_DIR
all_pass &= check(1, "Release package exists", rel_package.is_dir(), str(rel_package))

# Discover files dynamically
v2_csv_path = None
v2_bio_path = None
v2_fasta_path = None
v2_manifest_path = None
checksums_path = None

for f in REL_DIR.rglob("*"):
    if f.is_file():
        name = f.name
        if name == "FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv":
            v2_csv_path = f
        elif name == "FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER_V2.csv":
            v2_bio_path = f
        elif name == "FINAL_CANONICAL_KJU_16S_V2.fasta":
            v2_fasta_path = f
        elif name == "FINAL_KJU_MASTER_FREEZE_MANIFEST_V2.txt" and "MANIFESTS" in str(f):
            v2_manifest_path = f
        elif name == "STEP100_RELEASE_CHECKSUMS.txt":
            checksums_path = f

# Also check V2_DIR directly
if v2_csv_path is None:
    v2_csv_path = V2_DIR / "FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv"
if v2_bio_path is None:
    v2_bio_path = V2_DIR / "FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER_V2.csv"
if v2_fasta_path is None:
    v2_fasta_path = V2_DIR / "FINAL_CANONICAL_KJU_16S_V2.fasta"
if v2_manifest_path is None:
    v2_manifest_path = V2_DIR / "FINAL_KJU_MASTER_FREEZE_MANIFEST_V2.txt"

all_pass &= check(2, "V2 CSV identified", v2_csv_path is not None and v2_csv_path.is_file(),
                  str(v2_csv_path) if v2_csv_path else "NOT FOUND")
all_pass &= check(3, "V2 Bio identified", v2_bio_path is not None and v2_bio_path.is_file(),
                  str(v2_bio_path) if v2_bio_path else "NOT FOUND")
all_pass &= check(4, "V2 FASTA identified", v2_fasta_path is not None and v2_fasta_path.is_file(),
                  str(v2_fasta_path) if v2_fasta_path else "NOT FOUND")
all_pass &= check(5, "V2 manifest identified", v2_manifest_path is not None and v2_manifest_path.is_file(),
                  str(v2_manifest_path) if v2_manifest_path else "NOT FOUND")
all_pass &= check(6, "Checksums identified", checksums_path is not None and checksums_path.is_file(),
                  str(checksums_path) if checksums_path else "NOT FOUND")

# ──────────────────────────────────────────────────────────
# CHECK 7-8: SHA256 VERIFICATION
# ──────────────────────────────────────────────────────────
print("\nSHA256 VERIFICATION")
print("-" * 70)

v1_csv_path = V1_DIR / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"

v2_sha = sha256(v2_csv_path)
v1_sha = sha256(v1_csv_path)

# Parse manifest for expected V2 SHA
manifest_v2_sha = None
if v2_manifest_path:
    current_file = None
    with open(v2_manifest_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.endswith(".csv"):
                current_file = stripped
            elif stripped.startswith("SHA256:") and current_file:
                sha_val = stripped.split("SHA256:", 1)[1].strip()
                # Match the exact V2 master CSV filename
                if current_file == "FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv":
                    manifest_v2_sha = sha_val

v2_match = v2_sha == manifest_v2_sha if manifest_v2_sha else False
v1_match = v1_sha == V1_SHA256_EXPECTED

all_pass &= check(7, "V2 SHA256 matches manifest", v2_match,
                  f"calculated={v2_sha[:16]}... manifest={manifest_v2_sha[:16] if manifest_v2_sha else 'N/A'}...")
all_pass &= check(8, "V1 SHA256 matches expected", v1_match,
                  f"calculated={v1_sha[:16]}... expected={V1_SHA256_EXPECTED[:16]}...")

# ──────────────────────────────────────────────────────────
# CHECK 9-10: V2 MASTER COUNTS
# ──────────────────────────────────────────────────────────
print("\nV2 MASTER COUNTS")
print("-" * 70)

v2 = pd.read_csv(v2_csv_path, dtype=str)

all_pass &= check(9, "V2 row count", len(v2) == 204, f"{len(v2)} rows")
all_pass &= check(10, "V2 unique KJUs", v2["KJU_ID"].nunique() == 204,
                  f"{v2['KJU_ID'].nunique()} unique KJU IDs")

# ──────────────────────────────────────────────────────────
# CHECK 11-14: FASTA VERIFICATION
# ──────────────────────────────────────────────────────────
print("\nFASTA VERIFICATION")
print("-" * 70)

fasta_records = list(SeqIO.parse(str(v2_fasta_path), "fasta"))
# Extract KJU IDs from pipe-delimited headers
fasta_ids = [r.id.split("|")[0] if "|" in r.id else r.id for r in fasta_records]
fasta_seqs = [str(r.seq) for r in fasta_records]

all_pass &= check(11, "FASTA entry count", len(fasta_records) == 111,
                  f"{len(fasta_records)} entries")

# Check FASTA → V2 mapping
canon_kjus = set(v2[v2["Final_status"] == "ADJUDICATED_CANONICAL"]["KJU_ID"].values)
fasta_kjus = set(fasta_ids)
fasta_not_in_v2 = fasta_kjus - canon_kjus
v2_not_in_fasta = canon_kjus - fasta_kjus

all_pass &= check(12, "FASTA → V2 mapping", not fasta_not_in_v2 and not v2_not_in_fasta,
                  f"fasta_not_in_v2={len(fasta_not_in_v2)}, v2_not_in_fasta={len(v2_not_in_fasta)}")

all_pass &= check(13, "FASTA duplicate IDs", len(fasta_ids) == len(set(fasta_ids)),
                  f"duplicates={len(fasta_ids) - len(set(fasta_ids))}")

# ──────────────────────────────────────────────────────────
# CHECK 14: BIO UNIQUE COUNT
# ──────────────────────────────────────────────────────────
print("\nBIOLOGICAL UNIQUE MASTER")
print("-" * 70)

v2_bio = pd.read_csv(v2_bio_path, dtype=str)
all_pass &= check(14, "Bio unique count", len(v2_bio) == 111,
                  f"{len(v2_bio)} rows")

# ──────────────────────────────────────────────────────────
# CHECK 15: RELEASE CHECKSUMS
# ──────────────────────────────────────────────────────────
print("\nRELEASE CHECKSUMS")
print("-" * 70)

checksum_verified = 0
checksum_total = 0
checksum_failures = []
if checksums_path:
    with open(checksums_path) as fh:
        for line in fh:
            line = line.strip()
            if "  " in line and not line.startswith("=") and not line.startswith("-") and not line.startswith("V") and not line.startswith("Generated") and not line.startswith("Total"):
                parts = line.split("  ", 1)
                if len(parts) == 2:
                    expected_sha, rel_file = parts
                    # Skip self-referential checksum (circular reference)
                    if "STEP100_RELEASE_CHECKSUMS" in rel_file:
                        continue
                    fpath = REL_DIR / rel_file
                    if fpath.is_file():
                        actual_sha = sha256(fpath)
                        checksum_total += 1
                        if actual_sha == expected_sha:
                            checksum_verified += 1
                        else:
                            checksum_failures.append(rel_file)
                    else:
                        checksum_failures.append(f"{rel_file} (NOT FOUND)")

all_pass &= check(15, "Release checksums", checksum_verified == checksum_total,
                  f"{checksum_verified}/{checksum_total} verified" + (f", failures: {checksum_failures[:3]}" if checksum_failures else ""))

# ──────────────────────────────────────────────────────────
# CHECK 16: STATUS COUNTS
# ──────────────────────────────────────────────────────────
print("\nSTATUS COUNTS")
print("-" * 70)

status_counts = v2["Final_status"].value_counts()
n_canonical = int(status_counts.get("ADJUDICATED_CANONICAL", 0))
n_no_seq = int(status_counts.get("ADJUDICATED_NO_SEQUENCE", 0))
n_unresolved = int(status_counts.get("UNRESOLVED", 0))

status_ok = (n_canonical == 111 and n_no_seq == 1 and n_unresolved == 92)
all_pass &= check(16, "Status counts", status_ok,
                  f"canonical={n_canonical}, no_seq={n_no_seq}, unresolved={n_unresolved}")

# ──────────────────────────────────────────────────────────
# CHECK 17: FASTA SEQUENCE VALIDITY
# ──────────────────────────────────────────────────────────
print("\nFASTA SEQUENCE VALIDITY")
print("-" * 70)

invalid_seqs = [fasta_ids[i] for i, seq in enumerate(fasta_seqs) if not iupac_valid(seq)]
all_pass &= check(17, "FASTA IUPAC valid", not invalid_seqs,
                  f"invalid={len(invalid_seqs)} ({invalid_seqs[:5] if invalid_seqs else 'none'})")

# ──────────────────────────────────────────────────────────
# CHECK 18: FASTA LENGTH CONSISTENCY
# ──────────────────────────────────────────────────────────
print("\nFASTA LENGTH CONSISTENCY")
print("-" * 70)

length_mismatches = []
for i, kju in enumerate(fasta_ids):
    row = v2[v2["KJU_ID"] == kju]
    if len(row) > 0:
        row = row.iloc[0]
        fasta_len = len(fasta_seqs[i])
        try:
            v2_len = int(float(str(row.get("Canonical_sequence_length", "")).strip()))
            if fasta_len != v2_len:
                length_mismatches.append(kju)
        except (ValueError, TypeError):
            length_mismatches.append(kju)

all_pass &= check(18, "FASTA length consistency", not length_mismatches,
                  f"mismatches={len(length_mismatches)} ({length_mismatches[:5] if length_mismatches else 'none'})")

# ──────────────────────────────────────────────────────────
# CHECK 19: FASTA SHA256 CONSISTENCY
# ──────────────────────────────────────────────────────────
print("\nFASTA SHA256 CONSISTENCY")
print("-" * 70)

sha_mismatches = []
for i, kju in enumerate(fasta_ids):
    row = v2[v2["KJU_ID"] == kju]
    if len(row) > 0:
        row = row.iloc[0]
        fasta_sha = sha256_seq(fasta_seqs[i])
        v2_sha_val = str(row.get("Canonical_sequence_SHA256", "")).strip()
        if fasta_sha != v2_sha_val:
            sha_mismatches.append(kju)

all_pass &= check(19, "FASTA SHA256 consistency", not sha_mismatches,
                  f"mismatches={len(sha_mismatches)} ({sha_mismatches[:5] if sha_mismatches else 'none'})")

# ──────────────────────────────────────────────────────────
# CHECK 20: WORKSPACE CREATION
# ──────────────────────────────────────────────────────────
print("\nWORKSPACE CREATION")
print("-" * 70)

dirs_expected = [
    "STEP101_RELEASE_INTEGRITY",
    "STEP102_ANALYSIS_METADATA",
    "STEP103_TAXONOMIC_COMPOSITION",
    "STEP104_SEQUENCE_DIVERSITY",
    "STEP105_REPRODUCIBLE_ALIGNMENT",
    "STEP106_REPRODUCIBLE_PHYLOGENY",
    "STEP107_PHYLOGENETIC_TAXONOMIC_DIVERSITY",
]
dirs_exist = [(DOWNSTREAM / d).is_dir() for d in dirs_expected]
all_pass &= check(20, "Workspace directories", all(dirs_exist),
                  f"{sum(dirs_exist)}/{len(dirs_expected)} exist")

# ──────────────────────────────────────────────────────────
# FINAL REPORT
# ──────────────────────────────────────────────────────────
print()
print("=" * 70)
n_pass = sum(1 for _, _, s, _ in checks if s == "PASS")
n_fail = sum(1 for _, _, s, _ in checks if s == "FAIL")
overall = "PASS" if all_pass else "FAIL"

print(f"\nSTEP 101: {overall}")
print(f"  Checks: {n_pass}/{len(checks)} PASS, {n_fail} FAIL")
print()

# ──────────────────────────────────────────────────────────
# MACHINE-READABLE REPORT
# ──────────────────────────────────────────────────────────
report_rows = []
for num, name, status, detail in checks:
    report_rows.append({
        "Check_Number": num,
        "Check_Name": name,
        "Status": status,
        "Detail": detail,
    })
report_df = pd.DataFrame(report_rows)
report_df.to_csv(OUT / "STEP101_INTEGRITY_REPORT.csv", index=False)

# ──────────────────────────────────────────────────────────
# HUMAN-READABLE REPORT
# ──────────────────────────────────────────────────────────
with open(OUT / "STEP101_REPORT.txt", "w") as fh:
    fh.write("STEP 101: RELEASE/PACKAGE INTEGRITY CHECK\n")
    fh.write("=" * 70 + "\n\n")

    fh.write("OVERALL STATUS: " + overall + "\n")
    fh.write(f"Checks: {n_pass}/{len(checks)} PASS, {n_fail} FAIL\n\n")

    fh.write("CHECK RESULTS\n")
    fh.write("-" * 70 + "\n")
    for num, name, status, detail in checks:
        fh.write(f"  [{status}] Check {num:02d} — {name}: {detail}\n")

    fh.write("\nFILES VERIFIED\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  V2 CSV: {v2_csv_path}\n")
    fh.write(f"  V2 Bio: {v2_bio_path}\n")
    fh.write(f"  V2 FASTA: {v2_fasta_path}\n")
    fh.write(f"  V2 Manifest: {v2_manifest_path}\n")
    fh.write(f"  V1 CSV: {v1_csv_path}\n")
    fh.write(f"  Checksums: {checksums_path}\n")

    fh.write("\nSHA256 VALUES\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  V1: {v1_sha}\n")
    fh.write(f"  V2: {v2_sha}\n")
    fh.write(f"  V2 FASTA: {sha256(v2_fasta_path)}\n")
    fh.write(f"  V2 Bio: {sha256(v2_bio_path)}\n")

    fh.write("\nWORKSPACE CREATED\n")
    fh.write("-" * 70 + "\n")
    fh.write(f"  Root: {DOWNSTREAM}\n")
    for d in dirs_expected:
        fh.write(f"  {d}/\n")

    fh.write(f"\nGenerated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# ──────────────────────────────────────────────────────────
# DOWNSTREAM MANIFEST (SHA256 of all files steps 102-107 will read)
# ──────────────────────────────────────────────────────────
downstream_files = {
    "V2_MASTER": v2_csv_path,
    "V2_BIO": v2_bio_path,
    "V2_FASTA": v2_fasta_path,
    "V2_MANIFEST": v2_manifest_path,
    "V1_MASTER": v1_csv_path,
    "BLAST_RESULTS": RESULTS / "STEP92_BLAST_VALIDATION" / "STEP92_BLAST_RESULTS.csv",
    "BLAST_TAXONOMY_CORRECTED": RESULTS / "STEP92_BLAST_VALIDATION" / "STEP92b_TAXONOMY_CORRECTION" / "STEP92b_TAXONOMY_COMPARISON_CORRECTED.csv",
    "PHYLOGENY_TREE": RESULTS / "STEP93_PHYLOGENETIC_VALIDATION" / "STEP93_PHYLOGENY.treefile",
    "PHYLOGENY_ALIGNMENT": RESULTS / "STEP93_PHYLOGENETIC_VALIDATION" / "STEP93_ALIGNED_CANONICAL.fasta",
    "CLUSTERING_AUDIT": RESULTS / "STEP93_PHYLOGENETIC_VALIDATION" / "STEP93_CLUSTERING_AUDIT.csv",
    "CROSS_KJU": RESULTS / "STEP97_BIOLOGICAL_UNIQUENESS" / "STEP97_CROSS_KJU_RELATIONSHIPS.csv",
    "GENBANK_METADATA": BASE / "07_Master_Metadata" / "direct_genbank_kju_metadata.csv",
    "STEP83_RECORDS": RESULTS / "STEP83_TRUE_SEQUENCE_INVENTORY" / "STEP83_TRUE_SEQUENCE_RECORDS.csv",
    "STEP87_CANONICAL": RESULTS / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE" / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv",
}

with open(OUT / "STEP101_DOWNSTREAM_MANIFEST.txt", "w") as fh:
    fh.write("STEP 101 DOWNSTREAM MANIFEST\n")
    fh.write("=" * 70 + "\n\n")
    fh.write(f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write(f"V1 SHA256: {V1_SHA256_EXPECTED}\n")
    fh.write(f"V2 SHA256: {v2_sha}\n\n")

    fh.write("FILES FOR DOWNSTREAM ANALYSIS\n")
    fh.write("-" * 70 + "\n")
    for label, path in downstream_files.items():
        if path.is_file():
            fsha = sha256(path)
            fh.write(f"{fsha}  {label}: {path}\n")
        else:
            fh.write(f"MISSING  {label}: {path}\n")
    fh.write(f"\nTotal: {sum(1 for p in downstream_files.values() if p.is_file())}/{len(downstream_files)} files found\n")

# ──────────────────────────────────────────────────────────
# FINAL SUMMARY
# ──────────────────────────────────────────────────────────
print()
print("OUTPUTS:")
for f in sorted(OUT.glob("*")):
    print(f"  {f}")
print()

print("WORKSPACE:")
for d in dirs_expected:
    print(f"  {DOWNSTREAM / d}")
print()

print("=" * 70)
print(f"V2 integrity: {'PASS' if v2_match else 'FAIL'}")
print(f"V1 integrity: {'PASS' if v1_match else 'FAIL'}")
print(f"FASTA integrity: {'PASS' if not invalid_seqs and not length_mismatches and not sha_mismatches else 'FAIL'}")
print(f"Bio unique: {'PASS' if len(v2_bio) == 111 else 'FAIL'}")
print(f"Checksums: {'PASS' if checksum_verified == checksum_total else 'FAIL'}")
print(f"OVERALL: {overall}")
print()
print("STEP 101 COMPLETE")
