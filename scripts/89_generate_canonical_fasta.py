#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re


# ============================================================
# STEP 89
# CANONICAL SEQUENCE FASTA GENERATION
#
# PURPOSE:
#   Generate a multi-FASTA file with one entry per KJU
#   that has a canonical sequence.
# ============================================================


BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "final_integrated_master"
)

FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP87 = ROOT / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
OUT = ROOT / "STEP89_CANONICAL_FASTA"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"

FASTA_OUT = OUT / "STEP89_CANONICAL_KJU_16S.fasta"
MANIFEST = OUT / "STEP89_FASTA_MANIFEST.csv"
REPORT = OUT / "STEP89_FASTA_QC_REPORT.txt"

FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

VALID_IUPAC = set("ACGTNRYWSMKBDHVU")


def format_fasta_seq(seq, width=80):
    seq = seq.strip().upper()
    lines = []
    for i in range(0, len(seq), width):
        lines.append(seq[i:i + width])
    return "\n".join(lines)


print("=" * 110)
print("STEP 89: CANONICAL SEQUENCE FASTA GENERATION")
print("=" * 110)

for p in [MASTER_CSV, STEP87_CANONICAL]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
print(f"\nFrozen master SHA256: {master_hash}")
print(f"Expected SHA256:      {FROZEN_SHA}")
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print("Frozen master integrity: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
canonical = pd.read_csv(STEP87_CANONICAL, dtype=str, low_memory=False)

auto = canonical[
    canonical["Selection_status"].str.startswith("AUTO_SELECTED_", na=False)
].copy()

print(f"\nMaster rows:    {len(master)}")
print(f"Auto-selected:  {len(auto)}")

taxonomy_lookup = {}
for _, r in master.iterrows():
    taxonomy_lookup[r["KJU_ID"]] = str(r.get("Accepted_taxonomy", "")).strip()

manifest_rows = []
headers_seen = set()
fasta_errors = []

with open(FASTA_OUT, "w") as f:
    for _, row in auto.iterrows():
        kju = row["KJU_ID"]
        accession = str(row["Canonical_accession"]).strip()
        seq = str(row["Canonical_sequence"]).strip().upper()
        sha = str(row["Canonical_sequence_SHA256"]).strip()
        seq_len = int(row["Canonical_sequence_length"])
        tax = taxonomy_lookup.get(kju, "")
        tax_clean = re.sub(r"[^A-Za-z0-9]", "_", tax)

        header = f"{kju}|{accession}|{tax_clean}|length={seq_len}"

        if header in headers_seen:
            fasta_errors.append(f"Duplicate header: {header}")
        headers_seen.add(header)

        invalid_chars = [c for c in seq if c not in VALID_IUPAC]
        if invalid_chars:
            fasta_errors.append(
                f"{kju}: invalid IUPAC chars: {set(invalid_chars)}"
            )

        if len(seq) != seq_len:
            fasta_errors.append(
                f"{kju}: declared length {seq_len} != actual {len(seq)}"
            )

        calc_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
        if calc_sha != sha:
            fasta_errors.append(
                f"{kju}: SHA256 mismatch"
            )

        f.write(f">{header}\n")
        f.write(format_fasta_seq(seq) + "\n")

        manifest_rows.append({
            "KJU_ID": kju,
            "FASTA_header": header,
            "Accession": accession,
            "Taxonomy": tax,
            "Sequence_length": len(seq),
            "Sequence_SHA256": sha,
            "Selection_status": str(row["Selection_status"]),
            "Selection_evidence": str(row["Selection_evidence"]),
        })

manifest_df = pd.DataFrame(manifest_rows)
manifest_df.to_csv(MANIFEST, index=False)

print(f"\nFASTA written:  {FASTA_OUT}")
print(f"Manifest:       {MANIFEST}")


entry_count = len(auto)
header_count = len(headers_seen)
unique_kjus = auto["KJU_ID"].nunique()
total_bases = sum(len(str(r["Canonical_sequence"]).strip()) for _, r in auto.iterrows())
gc_count = sum(
    sum(1 for c in str(r["Canonical_sequence"]).strip().upper() if c in "GC")
    for _, r in auto.iterrows()
)
gc_pct = (gc_count / total_bases * 100) if total_bases > 0 else 0
min_len = min(len(str(r["Canonical_sequence"]).strip()) for _, r in auto.iterrows())
max_len = max(len(str(r["Canonical_sequence"]).strip()) for _, r in auto.iterrows())

print(f"\n=== FASTA QC ===")
print(f"Entries written:  {entry_count}")
print(f"Unique headers:   {header_count}")
print(f"Unique KJUs:      {unique_kjus}")
print(f"Total bases:      {total_bases}")
print(f"GC%:              {gc_pct:.1f}%")
print(f"Length range:     {min_len}-{max_len}")
print(f"Errors:           {len(fasta_errors)}")

if fasta_errors:
    for e in fasta_errors:
        print(f"  ERROR: {e}")
    raise ValueError(f"FASTA generation errors: {len(fasta_errors)}")

assert entry_count == header_count, "Header count mismatch"
assert entry_count == unique_kjus, "KJU count mismatch"

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(REPORT, "w") as f:
    f.write("STEP 89: CANONICAL SEQUENCE FASTA GENERATION\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: {master_hash}\n")
    f.write(f"Expected SHA256:      {FROZEN_SHA}\n")
    f.write(f"Status:               PASS\n\n")
    f.write(f"Entries:              {entry_count}\n")
    f.write(f"Unique headers:       {header_count}\n")
    f.write(f"Unique KJUs:          {unique_kjus}\n")
    f.write(f"Total bases:          {total_bases}\n")
    f.write(f"GC%:                  {gc_pct:.1f}%\n")
    f.write(f"Length range:         {min_len}-{max_len}\n")
    f.write(f"Errors:               {len(fasta_errors)}\n\n")
    f.write(f"FASTA: {FASTA_OUT}\n")
    f.write(f"Manifest: {MANIFEST}\n\n")
    f.write(f"Frozen V1 master untouched.\n")

print(f"Saved report: {REPORT}")

print("\n" + "=" * 110)
print("STEP 89 COMPLETE")
print("=" * 110)
