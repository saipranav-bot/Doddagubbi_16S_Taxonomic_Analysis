#!/usr/bin/env python3
"""STEP 105: Reproducible MAFFT Alignment"""

import hashlib
import subprocess
import time
from pathlib import Path

BASE = Path.home() / "Doddagubbi_Master"
OUT = BASE / "downstream_analysis/STEP105_REPRODUCIBLE_ALIGNMENT"
OUT.mkdir(parents=True, exist_ok=True)

V2_FASTA = BASE / "github_results/reconciliation_111_vs_109/final_integrated_master/STEP98_FINAL_FREEZE_V2/FINAL_KJU_MASTER_FREEZE_V2/FINAL_CANONICAL_KJU_16S_V2.fasta"
OUTPUT_FASTA = OUT / "STEP105_MAFFT_ALIGNED.fasta"

def sha256(p): return hashlib.sha256(p.read_bytes()).hexdigest()

print("STEP 105: REPRODUCIBLE MAFFT ALIGNMENT")
print("=" * 60)
print(f"Input: {V2_FASTA}")
print(f"Output: {OUTPUT_FASTA}")
print()

# Run MAFFT
print("Running MAFFT (FFT-NS-2)...")
start = time.time()

result = subprocess.run(
    ["mafft", "--auto", "--quiet", "--thread", "-1", str(V2_FASTA)],
    capture_output=True, text=True, timeout=600
)

elapsed = time.time() - start
print(f"MAFFT completed in {elapsed:.1f}s")

if result.returncode != 0:
    print(f"MAFFT STDERR: {result.stderr[:500]}")
    raise RuntimeError(f"MAFFT failed with return code {result.returncode}")

# Write alignment
OUTPUT_FASTA.write_text(result.stdout)

# Verify
from Bio import SeqIO
aln_records = list(SeqIO.parse(str(OUTPUT_FASTA), "fasta"))
print(f"Alignment sequences: {len(aln_records)}")
aln_len = len(aln_records[0].seq) if aln_records else 0
print(f"Alignment length: {aln_len} bp")

# Alignment statistics
gap_counts = [str(r.seq).count("-") for r in aln_records]
conserved = sum(1 for i in range(aln_len) if len(set(str(r.seq)[i] for r in aln_records if str(r.seq)[i] != "-")) == 1)
print(f"Mean gaps per sequence: {sum(gap_counts)/len(gap_counts):.1f}")
print(f"Conserved positions: {conserved}/{aln_len} ({conserved/aln_len*100:.1f}%)")

# Save stats
import pandas as pd
stats = pd.DataFrame([{
    "Input_sequences": 111,
    "Output_sequences": len(aln_records),
    "Alignment_length": aln_len,
    "Mean_gaps_per_seq": round(sum(gap_counts)/len(gap_counts), 1),
    "Conserved_positions": conserved,
    "Conserved_percent": round(conserved/aln_len*100, 1),
    "Runtime_seconds": round(elapsed, 1),
}])
stats.to_csv(OUT / "STEP105_ALIGNMENT_STATS.csv", index=False)

# Checksums
for f in sorted(OUT.glob("STEP105_*")):
    if f.is_file():
        with open(OUT / "STEP105_CHECKSUMS.txt", "a") as fh:
            fh.write(f"{sha256(f)}  {f.name}\n")

with open(OUT / "STEP105_REPORT.txt", "w") as fh:
    fh.write("STEP 105: REPRODUCIBLE MAFFT ALIGNMENT\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Input: {V2_FASTA}\n")
    fh.write(f"Method: MAFFT --auto (FFT-NS-2)\n")
    fh.write(f"Runtime: {elapsed:.1f}s\n\n")
    fh.write(f"Sequences: {len(aln_records)}\n")
    fh.write(f"Alignment length: {aln_len} bp\n")
    fh.write(f"Mean gaps per sequence: {sum(gap_counts)/len(gap_counts):.1f}\n")
    fh.write(f"Conserved positions: {conserved}/{aln_len} ({conserved/aln_len*100:.1f}%)\n")

print("\nOUTPUTS:")
for f in sorted(OUT.glob("STEP105_*")):
    print(f"  {f.name}")
print("\nSTEP 105 COMPLETE")
