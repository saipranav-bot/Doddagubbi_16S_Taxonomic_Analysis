#!/usr/bin/env python3
"""STEP 106: Reproducible IQ-TREE Phylogeny"""

import hashlib
import subprocess
import time
from pathlib import Path

BASE = Path.home() / "Doddagubbi_Master"
OUT = BASE / "downstream_analysis/STEP106_REPRODUCIBLE_PHYLOGENY"
OUT.mkdir(parents=True, exist_ok=True)

ALIGNMENT = BASE / "downstream_analysis/STEP105_REPRODUCIBLE_ALIGNMENT/STEP105_MAFFT_ALIGNED.fasta"

def sha256(p): return hashlib.sha256(p.read_bytes()).hexdigest()

print("STEP 106: REPRODUCIBLE IQ-TREE PHYLOGENY")
print("=" * 60)
print(f"Input: {ALIGNMENT}")
print()

# Run IQ-TREE
print("Running IQ-TREE (GTR+G)...")
start = time.time()

result = subprocess.run(
    [
        "iqtree",
        "-s", str(ALIGNMENT),
        "-m", "GTR+G",
        "-pre", str(OUT / "STEP106_PHYLOGENY"),
        "-nt", "AUTO",
        "-seed", "12345",
        "-redo",
        "-quiet",
    ],
    capture_output=True, text=True, timeout=1800
)

elapsed = time.time() - start
print(f"IQ-TREE completed in {elapsed:.1f}s")

if result.returncode != 0:
    print(f"IQ-TREE STDERR: {result.stderr[:1000]}")
    # Try without -redo
    print("Retrying without -redo...")
    result = subprocess.run(
        [
            "iqtree",
            "-s", str(ALIGNMENT),
            "-m", "GTR+G",
            "-pre", str(OUT / "STEP106_PHYLOGENY"),
            "-nt", "AUTO",
            "-seed", "12345",
            "-quiet",
        ],
        capture_output=True, text=True, timeout=1800
    )
    elapsed = time.time() - start
    if result.returncode != 0:
        print(f"Retry STDERR: {result.stderr[:1000]}")
        raise RuntimeError(f"IQ-TREE failed with return code {result.returncode}")

# Check outputs
treefile = OUT / "STEP106_PHYLOGENY.treefile"
if not treefile.exists():
    # Look for .contree
    contree = OUT / "STEP106_PHYLOGENY.contree"
    if contree.exists():
        treefile = contree

if treefile.exists():
    tree_text = treefile.read_text()
    print(f"Tree file: {treefile.name} ({len(tree_text)} chars)")
    # Count tips
    import re
    tips = re.findall(r"[A-Z0-9_]+\|[A-Z0-9]+\|[A-Za-z_]+\|length_\d+", tree_text)
    print(f"Tree tips: {len(tips)}")
else:
    print("WARNING: No tree file found")
    tree_text = ""

# Save stats
import pandas as pd
stats = pd.DataFrame([{
    "Input_sequences": 111,
    "Model": "GTR+G",
    "Runtime_seconds": round(elapsed, 1),
    "Tree_file": str(treefile.name) if treefile.exists() else "NONE",
}])
stats.to_csv(OUT / "STEP106_PHYLOGENY_STATS.csv", index=False)

# Checksums
for f in sorted(OUT.glob("STEP106_*")):
    if f.is_file():
        with open(OUT / "STEP106_CHECKSUMS.txt", "a") as fh:
            fh.write(f"{sha256(f)}  {f.name}\n")

with open(OUT / "STEP106_REPORT.txt", "w") as fh:
    fh.write("STEP 106: REPRODUCIBLE IQ-TREE PHYLOGENY\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Input: {ALIGNMENT}\n")
    fh.write(f"Model: GTR+G\n")
    fh.write(f"Seed: 12345\n")
    fh.write(f"Runtime: {elapsed:.1f}s\n\n")
    fh.write(f"Tree tips: {len(tips)}\n")
    fh.write(f"Tree file: {treefile.name}\n")
    fh.write(f"IQ-TREE log: STEP106_PHYLOGENY.iqtree\n")

print("\nOUTPUTS:")
for f in sorted(OUT.glob("STEP106_*")):
    print(f"  {f.name}")
print("\nSTEP 106 COMPLETE")
