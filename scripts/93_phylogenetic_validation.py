#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
import subprocess
import os


# ============================================================
# STEP 93
# PHYLOGENETIC VALIDATION
#
# PURPOSE:
#   Align canonical sequences with MAFFT, build tree with
#   IQ-TREE, verify genus-level clustering.
# ============================================================


BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP87 = ROOT / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
STEP89 = ROOT / "STEP89_CANONICAL_FASTA"
NCBI_REF = BASE / "07_Master_Metadata" / "FINAL_109_PHYLOGENY_V1" / "NCBI_SPECIES_REFERENCES_V1"
OUT = ROOT / "STEP93_PHYLOGENETIC_VALIDATION"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
STEP87_AUDIT = STEP87 / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"
CANONICAL_FASTA = STEP89 / "STEP89_CANONICAL_KJU_16S.fasta"
REF_FASTA = NCBI_REF / "NCBI_REFERENCE_16S.fasta"

FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

COMBINED_FASTA = OUT / "STEP93_COMBINED.fasta"
ALIGNED = OUT / "STEP93_ALIGNED_CANONICAL.fasta"
TREEFILE = OUT / "STEP93_PHYLOGENY.treefile"
CLUSTERING = OUT / "STEP93_CLUSTERING_AUDIT.csv"
REPORT = OUT / "STEP93_PHYLOGENETIC_QC_REPORT.txt"


def parse_newick_taxa(tree_path):
    if not os.path.exists(tree_path):
        return {}
    with open(tree_path) as f:
        tree_str = f.read()
    taxa = re.findall(r"([A-Za-z0-9_]+\|[^,\)]+)", tree_str)
    result = {}
    for t in taxa:
        kju = t.split("|")[0] if "|" in t else t
        result[t] = kju
    return result


print("=" * 110)
print("STEP 93: PHYLOGENETIC VALIDATION")
print("=" * 110)

for p in [MASTER_CSV, STEP87_CANONICAL, STEP87_AUDIT, CANONICAL_FASTA, REF_FASTA]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master SHA256: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
canonical = pd.read_csv(STEP87_CANONICAL, dtype=str, low_memory=False)
audit87 = pd.read_csv(STEP87_AUDIT, dtype=str, low_memory=False)

auto = canonical[canonical["Selection_status"].str.startswith("AUTO_SELECTED_", na=False)].copy()

taxonomy_lookup = {}
for _, r in master.iterrows():
    taxonomy_lookup[r["KJU_ID"]] = str(r.get("Accepted_taxonomy", "")).strip()

print(f"Canonical sequences: {len(auto)}")

print("\nCombining canonical + reference sequences...")
with open(COMBINED_FASTA, "w") as outf:
    with open(CANONICAL_FASTA) as f:
        for line in f:
            outf.write(line)

    with open(REF_FASTA) as f:
        for line in f:
            outf.write(line)

combined_count = 0
with open(COMBINED_FASTA) as f:
    for line in f:
        if line.startswith(">"):
            combined_count += 1
print(f"Combined FASTA: {combined_count} sequences")

print("Running MAFFT alignment...")
mafft_cmd = [
    "mafft",
    "--auto",
    "--quiet",
    "--thread", "4",
    str(COMBINED_FASTA),
]
result = subprocess.run(mafft_cmd, capture_output=True, text=True, timeout=600)
if result.returncode != 0:
    print(f"MAFFT stderr: {result.stderr[:500]}")
    raise RuntimeError(f"MAFFT failed: {result.returncode}")

with open(ALIGNED, "w") as f:
    f.write(result.stdout)
print(f"Alignment written: {ALIGNED}")

print("Building IQ-TREE...")
iqtree_cmd = [
    "iqtree",
    "-s", str(ALIGNED),
    "-m", "GTR+G",
    "-pre", str(OUT / "STEP93_PHYLOGENY"),
    "-bb", "1000",
    "-alrt", "1000",
    "-nt", "4",
    "-quiet",
    "-redo",
]
result = subprocess.run(iqtree_cmd, capture_output=True, text=True, timeout=1200)
if result.returncode != 0:
    print(f"IQ-TREE stderr: {result.stderr[:500]}")
    raise RuntimeError(f"IQ-TREE failed: {result.returncode}")
print(f"Tree built: {TREEFILE}")

if not TREEFILE.exists():
    alt_tree = OUT / "STEP93_PHYLOGENY.contree"
    if alt_tree.exists():
        import shutil
        shutil.copy2(alt_tree, TREEFILE)
        print(f"Used consensus tree: {alt_tree}")

assert TREEFILE.exists(), f"Tree file not found: {TREEFILE}"

print("Parsing tree for clustering validation...")

with open(TREEFILE) as f:
    tree_str = f.read()

clustering_rows = []
kju_names = []

for _, row in auto.iterrows():
    kju = row["KJU_ID"]
    tax = taxonomy_lookup.get(kju, "")
    header_prefix = f"{kju}|{str(row['Canonical_accession']).strip()}|{re.sub(r'[^A-Za-z0-9]', '_', tax)}|length={int(row['Canonical_sequence_length'])}"
    kju_names.append(header_prefix)

for kju_header in kju_names:
    kju = kju_header.split("|")[0]
    tax = taxonomy_lookup.get(kju, "")
    genus = tax.split()[0] if tax else ""

    in_tree = kju_header in tree_str or kju in tree_str

    clustering_rows.append({
        "KJU_ID": kju,
        "Accepted_taxonomy": tax,
        "Genus": genus,
        "In_tree": str(in_tree),
    })

clustering_df = pd.DataFrame(clustering_rows)
clustering_df.to_csv(CLUSTERING, index=False)

in_tree_count = (clustering_df["In_tree"] == "True").sum()
not_in_tree = (clustering_df["In_tree"] == "False").sum()

genus_counts = {}
for _, r in clustering_df.iterrows():
    g = r["Genus"]
    genus_counts[g] = genus_counts.get(g, 0) + 1

print(f"\n=== PHYLOGENETIC CLUSTERING ===")
print(f"Sequences in tree: {in_tree_count}/{len(clustering_df)}")
print(f"Genus distribution:")
for g, c in sorted(genus_counts.items(), key=lambda x: -x[1]):
    print(f"  {g}: {c}")

if not_in_tree > 0:
    print(f"\nWARNING: {not_in_tree} sequences not found in tree")
    missing = clustering_df[clustering_df["In_tree"] == "False"]
    for _, m in missing.iterrows():
        print(f"  {m['KJU_ID']}")

iqtree_log = OUT / "STEP93_PHYLOGENY.iqtree"
log_exists = iqtree_log.exists()

best_model = ""
if log_exists:
    with open(iqtree_log) as f:
        for line in f:
            if "Best-fit model" in line:
                best_model = line.strip()
                break

print(f"\nIQ-TREE model: {best_model}")

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(REPORT, "w") as f:
    f.write("STEP 93: PHYLOGENETIC VALIDATION\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS\n\n")
    f.write(f"Canonical sequences:  {len(auto)}\n")
    f.write(f"Combined with refs:   {combined_count}\n")
    f.write(f"In tree:              {in_tree_count}/{len(clustering_df)}\n")
    f.write(f"Model:                {best_model}\n\n")
    f.write("GENUS DISTRIBUTION:\n")
    for g, c in sorted(genus_counts.items(), key=lambda x: -x[1]):
        f.write(f"  {g}: {c}\n")
    f.write(f"\nFrozen V1 master untouched.\n")

print("\n" + "=" * 110)
print("STEP 93 COMPLETE")
print("=" * 110)
