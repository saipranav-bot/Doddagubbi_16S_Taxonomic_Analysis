#!/usr/bin/env python3
"""STEP 104: Sequence-Level Diversity Analysis"""

import hashlib
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from Bio import SeqIO

BASE = Path.home() / "Doddagubbi_Master"
OUT = BASE / "downstream_analysis/STEP104_SEQUENCE_DIVERSITY"
OUT.mkdir(parents=True, exist_ok=True)

V2_FASTA = BASE / "github_results/reconciliation_111_vs_109/final_integrated_master/STEP98_FINAL_FREEZE_V2/FINAL_KJU_MASTER_FREEZE_V2/FINAL_CANONICAL_KJU_16S_V2.fasta"

def sha256(p): return hashlib.sha256(p.read_bytes()).hexdigest()

print("STEP 104: SEQUENCE-LEVEL DIVERSITY")
print("=" * 60)

# Load sequences
records = []
for rec in SeqIO.parse(str(V2_FASTA), "fasta"):
    kju = rec.id.split("|")[0] if "|" in rec.id else rec.id
    seq = str(rec.seq).upper()
    records.append({"KJU_ID": kju, "Sequence": seq, "Length": len(seq)})

df = pd.DataFrame(records)

# GC content
df["GC"] = df["Sequence"].apply(lambda s: round((s.count("G") + s.count("C")) / len(s) * 100, 2))

# Nucleotide composition
def nuc_comp(seq):
    c = Counter(seq)
    total = len(seq)
    return {n: round(c.get(n, 0) / total * 100, 2) for n in "ATGC"}

nuc_df = pd.DataFrame([nuc_comp(s) for s in df["Sequence"]])
df = pd.concat([df, nuc_df], axis=1)

# Summary statistics
print(f"Sequences: {len(df)}")
print(f"Length range: {df['Length'].min()} – {df['Length'].max()} bp")
print(f"Length mean: {df['Length'].mean():.0f} bp")
print(f"Length std: {df['Length'].std():.1f} bp")
print(f"GC range: {df['GC'].min():.1f}% – {df['GC'].max():.1f}%")
print(f"GC mean: {df['GC'].mean():.1f}%")
print(f"GC std: {df['GC'].std():.2f}%")
print()
print("Nucleotide composition (mean):")
for n in "ATGC":
    print(f"  {n}: {df[n].mean():.1f}%")

# Save per-sequence stats
df[["KJU_ID", "Length", "GC", "A", "T", "G", "C"]].to_csv(OUT / "STEP104_SEQUENCE_STATS.csv", index=False)

# Pairwise identity matrix (111x111)
print("\nComputing pairwise identity matrix...")
seqs = df["Sequence"].tolist()
kjus = df["KJU_ID"].tolist()
n = len(seqs)
identity_matrix = np.zeros((n, n), dtype=float)

for i in range(n):
    identity_matrix[i, i] = 100.0
    for j in range(i + 1, n):
        matches = sum(a == b for a, b in zip(seqs[i], seqs[j]))
        ident = round(matches / min(len(seqs[i]), len(seqs[j])) * 100, 2)
        identity_matrix[i, j] = ident
        identity_matrix[j, i] = ident

# Save identity matrix
id_df = pd.DataFrame(identity_matrix, index=kjus, columns=kjus)
id_df.to_csv(OUT / "STEP104_PAIRWISE_IDENTITY_MATRIX.csv")

# Save summary of pairwise identities
upper_tri = identity_matrix[np.triu_indices(n, k=1)]
print(f"Pairwise identity range: {upper_tri.min():.2f}% – {upper_tri.max():.2f}%")
print(f"Pairwise identity mean: {upper_tri.mean():.2f}%")
print(f"Pairwise identity median: {np.median(upper_tri):.2f}%")

# Distance matrix (100 - identity)
dist_matrix = 100.0 - identity_matrix
np.fill_diagonal(dist_matrix, 0)
dist_df = pd.DataFrame(dist_matrix, index=kjus, columns=kjus)
dist_df.to_csv(OUT / "STEP104_DISTANCE_MATRIX.csv")

# Save summary
summary = pd.DataFrame([{
    "Metric": "Sequences", "Value": len(df)
}, {
    "Metric": "Length_min", "Value": int(df["Length"].min())
}, {
    "Metric": "Length_max", "Value": int(df["Length"].max())
}, {
    "Metric": "Length_mean", "Value": round(df["Length"].mean(), 1)
}, {
    "Metric": "Length_std", "Value": round(df["Length"].std(), 2)
}, {
    "Metric": "GC_min", "Value": df["GC"].min()
}, {
    "Metric": "GC_max", "Value": df["GC"].max()
}, {
    "Metric": "GC_mean", "Value": round(df["GC"].mean(), 2)
}, {
    "Metric": "GC_std", "Value": round(df["GC"].std(), 2)
}, {
    "Metric": "Pairwise_identity_min", "Value": round(upper_tri.min(), 2)
}, {
    "Metric": "Pairwise_identity_max", "Value": round(upper_tri.max(), 2)
}, {
    "Metric": "Pairwise_identity_mean", "Value": round(upper_tri.mean(), 2)
}, {
    "Metric": "Pairwise_identity_median", "Value": round(np.median(upper_tri), 2)
}])
summary.to_csv(OUT / "STEP104_SUMMARY.csv", index=False)

# Checksums
for f in sorted(OUT.glob("STEP104_*")):
    if f.is_file():
        with open(OUT / "STEP104_CHECKSUMS.txt", "a") as fh:
            fh.write(f"{sha256(f)}  {f.name}\n")

with open(OUT / "STEP104_REPORT.txt", "w") as fh:
    fh.write("STEP 104: SEQUENCE-LEVEL DIVERSITY\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Input: {V2_FASTA}\n")
    fh.write(f"Sequences: {len(df)}\n\n")
    fh.write("LENGTH DISTRIBUTION:\n")
    fh.write(f"  Min: {int(df['Length'].min())} bp\n")
    fh.write(f"  Max: {int(df['Length'].max())} bp\n")
    fh.write(f"  Mean: {df['Length'].mean():.0f} bp\n")
    fh.write(f"  Std: {df['Length'].std():.1f} bp\n\n")
    fh.write("GC CONTENT:\n")
    fh.write(f"  Min: {df['GC'].min():.1f}%\n")
    fh.write(f"  Max: {df['GC'].max():.1f}%\n")
    fh.write(f"  Mean: {df['GC'].mean():.1f}%\n")
    fh.write(f"  Std: {df['GC'].std():.2f}%\n\n")
    fh.write("NUCLEOTIDE COMPOSITION (mean):\n")
    for n in "ATGC":
        fh.write(f"  {n}: {df[n].mean():.1f}%\n")
    fh.write(f"\nPAIRWISE IDENTITY:\n")
    fh.write(f"  Min: {upper_tri.min():.2f}%\n")
    fh.write(f"  Max: {upper_tri.max():.2f}%\n")
    fh.write(f"  Mean: {upper_tri.mean():.2f}%\n")
    fh.write(f"  Median: {np.median(upper_tri):.2f}%\n")

print("\nOUTPUTS:")
for f in sorted(OUT.glob("STEP104_*")):
    print(f"  {f.name}")
print("\nSTEP 104 COMPLETE")
