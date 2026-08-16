#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import shutil

BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP95 = ROOT / "STEP95_FINAL_INTEGRATED_MASTER_V2"
STEP89 = ROOT / "STEP89_CANONICAL_FASTA"
OUT = ROOT / "STEP98_FINAL_FREEZE_V2"
FREEZE_V2 = OUT / "FINAL_KJU_MASTER_FREEZE_V2"
FREEZE_V2.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
V2_CSV = STEP95 / "STEP95_INTEGRATED_MASTER_V2.csv"
CANONICAL_FASTA = STEP89 / "STEP89_CANONICAL_KJU_16S.fasta"
FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

V2_MASTER = FREEZE_V2 / "FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv"
V2_BIO = FREEZE_V2 / "FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER_V2.csv"
V2_FASTA = FREEZE_V2 / "FINAL_CANONICAL_KJU_16S_V2.fasta"
V2_MANIFEST = FREEZE_V2 / "FINAL_KJU_MASTER_FREEZE_MANIFEST_V2.txt"

print("=" * 110)
print("STEP 98: FINAL FREEZE V2")
print("=" * 110)

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master V1 SHA256: PASS")

v2 = pd.read_csv(V2_CSV, dtype=str, low_memory=False)
print(f"V2 rows: {len(v2)}")
assert len(v2) == 204
assert v2["KJU_ID"].nunique() == 204

shutil.copy2(V2_CSV, V2_MASTER)
print(f"V2 master copied: {V2_MASTER}")

bio_rows = []
seen_bio = set()
for _, row in v2.iterrows():
    bio = str(row.get("Final_biological_identity", "")).strip()
    if bio in {"", "nan", "none"}:
        continue
    if bio in seen_bio:
        continue
    seen_bio.add(bio)
    bio_rows.append(row.to_dict())

bio_df = pd.DataFrame(bio_rows)
bio_df.to_csv(V2_BIO, index=False)
print(f"Biologically unique V2: {len(bio_df)} rows")

shutil.copy2(CANONICAL_FASTA, V2_FASTA)
print(f"Canonical FASTA copied: {V2_FASTA}")

v2_master_sha = hashlib.sha256(V2_MASTER.read_bytes()).hexdigest()
v2_bio_sha = hashlib.sha256(V2_BIO.read_bytes()).hexdigest()
v2_fasta_sha = hashlib.sha256(V2_FASTA.read_bytes()).hexdigest()

master_hash_final = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_final == FROZEN_SHA, "V1 MASTER CHANGED"

status_counts = v2["Final_status"].value_counts()
canon_count = sum(v for k, v in status_counts.items() if "CANONICAL" in k)
noseq_count = status_counts.get("ADJUDICATED_NO_SEQUENCE", 0)
unresolved_count = status_counts.get("UNRESOLVED", 0)

with open(V2_MANIFEST, "w") as f:
    f.write("FINAL KJU MASTER FREEZE V2 MANIFEST\n")
    f.write("=" * 90 + "\n\n")
    f.write(f"Freeze date: 2026-08-16\n\n")
    f.write("FROZEN V1 INTEGRITY:\n")
    f.write(f"  SHA256: {master_hash_final}\n")
    f.write(f"  Status: PASS\n\n")
    f.write("V2 ARTIFACTS:\n")
    f.write(f"  FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv\n")
    f.write(f"    SHA256: {v2_master_sha}\n")
    f.write(f"    Rows: {len(v2)}\n")
    f.write(f"    Unique KJUs: {v2['KJU_ID'].nunique()}\n\n")
    f.write(f"  FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER_V2.csv\n")
    f.write(f"    SHA256: {v2_bio_sha}\n")
    f.write(f"    Unique biological identities: {len(bio_df)}\n\n")
    f.write(f"  FINAL_CANONICAL_KJU_16S_V2.fasta\n")
    f.write(f"    SHA256: {v2_fasta_sha}\n\n")
    f.write("V2 STATISTICS:\n")
    f.write(f"  Total KJUs: 204\n")
    f.write(f"  Canonical sequences: {canon_count}\n")
    f.write(f"  No sequence: {noseq_count}\n")
    f.write(f"  Unresolved: {unresolved_count}\n\n")
    f.write("V1 FROZEN MASTER DIRECTORY: NEVER MODIFIED\n")
    f.write("SOURCE FASTA FILES: NEVER MODIFIED\n")
    f.write("NO SEQUENCE FABRICATED\n")
    f.write("NO ACCESSION FABRICATED\n")
    f.write("NO TAXONOMY FABRICATED\n")

print(f"\nManifest: {V2_MANIFEST}")

print(f"\n=== V2 FREEZE SUMMARY ===")
print(f"V2 master:       {V2_MASTER} ({v2_master_sha[:16]}...)")
print(f"V2 biologically unique: {V2_BIO} ({v2_bio_sha[:16]}...)")
print(f"V2 canonical FASTA: {V2_FASTA} ({v2_fasta_sha[:16]}...)")
print(f"V1 frozen SHA256: {master_hash_final[:16]}... (unchanged)")

print("\n" + "=" * 110)
print("STEP 98 COMPLETE — V2 FROZEN")
print("=" * 110)
