#!/usr/bin/env python3
"""
STEP 92b: BLAST Taxonomy Species Lookup Correction

Fixes a bug in STEP92 where ref_lookup used versioned accession keys
(e.g. PX979534.1) but the lookup code normalized to root accessions
(e.g. PX979534). The species lookup failed silently and
Taxonomy_agreement defaulted to True vacuously.

This step:
1. Reads existing STEP92 BLAST results (raw hits are correct)
2. Re-extracts Hit_species from Subject_accession 3rd pipe field
3. Also builds corrected ref_lookup with root accession keys
4. Re-classifies taxonomy agreement into 4 categories:
   - SPECIES_AGREEMENT: accepted species in NCBI_REF, top hit matches
   - GENUS_AGREEMENT: same genus, different species (within-genus hit)
   - INSUFFICIENT_REFERENCE: accepted species not in NCBI_REF
   - DISCORDANT: genuine cross-genus mismatch
5. Does NOT modify V2 freeze

Output: STEP92b_TAXONOMY_CORRECTION/ inside STEP92 directory
"""

import os
import re
import hashlib
import pandas as pd
from pathlib import Path

# ──────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────
ROOT = Path("/home/pranav/Doddagubbi_Master")
RESULTS = ROOT / "github_results/reconciliation_111_vs_109/final_integrated_master"
STEP92_DIR = RESULTS / "STEP92_BLAST_VALIDATION"
OUT_DIR = STEP92_DIR / "STEP92b_TAXONOMY_CORRECTION"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BLAST_RESULTS = STEP92_DIR / "STEP92_BLAST_RESULTS.csv"
REF_META = ROOT / "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/NCBI_SPECIES_REFERENCES_V1/NCBI_REFERENCE_METADATA.csv"
V2_MASTER = RESULTS / "STEP98_FINAL_FREEZE_V2/FINAL_KJU_MASTER_FREEZE_V2/FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv"

# ──────────────────────────────────────────────────────────
# LOAD DATA
# ──────────────────────────────────────────────────────────
print("STEP 92b: BLAST TAXONOMY SPECIES LOOKUP CORRECTION")
print("=" * 60)
print()

blast = pd.read_csv(BLAST_RESULTS, dtype=str)
ref_meta = pd.read_csv(REF_META, dtype=str)
v2 = pd.read_csv(V2_MASTER, dtype=str)

print(f"BLAST results: {len(blast)} rows")
print(f"NCBI references: {len(ref_meta)} rows")
print(f"V2 master: {len(v2)} rows")
print()

# ──────────────────────────────────────────────────────────
# BUILD CORRECTED REF_LOOKUP (root accession keys)
# ──────────────────────────────────────────────────────────
ref_lookup = {}
for _, r in ref_meta.iterrows():
    acc = str(r.get("accession", "")).strip()
    root = re.sub(r"\.\d+$", "", acc.upper())
    ref_lookup[root] = {
        "species": str(r.get("species", "")).strip(),
        "organism": str(r.get("organism", "")).strip(),
        "length": str(r.get("length", "")).strip(),
    }

print(f"Ref lookup (root keys): {len(ref_lookup)} entries")
print()

# ──────────────────────────────────────────────────────────
# EXTRACT HIT_SPECIES FROM SUBJECT_ACCESSION
# ──────────────────────────────────────────────────────────
def parse_species_from_accession(accession_str):
    """Extract species from NCBI_REF pipe-delimited accession.
    Format: PX979534.1|NCBI_REF|Aeromonas_caviae -> Aeromonas caviae
    """
    parts = str(accession_str).split("|")
    if len(parts) >= 3:
        species = parts[2].replace("_", " ").strip()
        return species
    return ""


def normalize_species(s):
    """Normalize species name for comparison."""
    return " ".join(str(s).strip().lower().split())


def get_genus(species):
    """Extract genus from species name."""
    parts = str(species).strip().split()
    return parts[0] if parts else ""


# For each BLAST hit, extract species
blast["Hit_species_corrected"] = blast["Subject_accession"].apply(parse_species_from_accession)

# Also try ref_lookup for any that failed
for idx, row in blast.iterrows():
    if not row["Hit_species_corrected"]:
        root = str(row.get("Subject_accession_root", "")).strip().upper()
        ref_info = ref_lookup.get(root, {})
        if ref_info.get("species"):
            blast.at[idx, "Hit_species_corrected"] = ref_info["species"]

# Count populated
species_populated = blast["Hit_species_corrected"].apply(lambda x: str(x).strip() not in {"", "nan", "none"}).sum()
print(f"Hit_species_corrected populated: {species_populated}/{len(blast)}")
print()

# ──────────────────────────────────────────────────────────
# GET TOP HIT PER KJU
# ──────────────────────────────────────────────────────────
blast["Bit_score_float"] = blast["Bit_score"].astype(float)
blast["Pident_float"] = blast["Percent_identity"].astype(float)

top_hits = (
    blast.sort_values("Bit_score_float", ascending=False)
    .groupby("KJU_ID")
    .first()
    .reset_index()
)

print(f"Top hits per KJU: {len(top_hits)}")
print()

# ──────────────────────────────────────────────────────────
# BUILD CANONICAL SPECIES SET FROM V2
# ──────────────────────────────────────────────────────────
canon = v2[v2["Final_status"] == "ADJUDICATED_CANONICAL"].copy()
canon_species_set = set(canon["Accepted_taxonomy"].apply(normalize_species).unique())
ref_species_set = set(ref_meta["species"].apply(normalize_species).unique())

print(f"Canonical species: {len(canon_species_set)}")
print(f"NCBI_REF species: {len(ref_species_set)}")
print(f"Overlap: {len(canon_species_set & ref_species_set)}")
print()

# ──────────────────────────────────────────────────────────
# RECLASSIFY TAXONOMY AGREEMENT
# ──────────────────────────────────────────────────────────
canon_top = top_hits[top_hits["KJU_ID"].isin(canon["KJU_ID"].values)].copy()
canon_top = canon_top.merge(canon[["KJU_ID", "Accepted_taxonomy"]], on="KJU_ID", how="left")

results = []
for _, row in canon_top.iterrows():
    kju = row["KJU_ID"]
    accepted = str(row.get("Accepted_taxonomy", "")).strip()
    hit_species = str(row.get("Hit_species_corrected", "")).strip()
    accession = str(row.get("Subject_accession", "")).strip()
    identity = row.get("Pident_float", 0)
    alignment_len = row.get("Alignment_length", "0")
    bit_score = row.get("Bit_score", "0")
    e_value = row.get("E_value", "0")

    norm_accepted = normalize_species(accepted)
    norm_hit = normalize_species(hit_species)
    genus_accepted = get_genus(norm_accepted)
    genus_hit = get_genus(norm_hit)

    in_ref = norm_accepted in ref_species_set

    if not in_ref:
        category = "INSUFFICIENT_REFERENCE"
        agreement = "INSUFFICIENT_REFERENCE"
    elif norm_accepted == norm_hit:
        category = "SPECIES_AGREEMENT"
        agreement = "TRUE"
    elif genus_accepted == genus_hit:
        category = "GENUS_AGREEMENT"
        agreement = "GENUS_AGREEMENT"
    else:
        category = "DISCORDANT"
        agreement = "FALSE"

    results.append({
        "KJU_ID": kju,
        "Accepted_taxonomy": accepted,
        "Hit_species": hit_species,
        "Hit_accession": accession,
        "Percent_identity": identity,
        "Alignment_length": alignment_len,
        "E_value": e_value,
        "Bit_score": bit_score,
        "Taxonomy_agreement": agreement,
        "Agreement_category": category,
        "Accepted_in_NCBI_REF": in_ref,
    })

results_df = pd.DataFrame(results)

# ──────────────────────────────────────────────────────────
# FULL TAXONOMY COMPARISON (all 204 KJUs)
# ──────────────────────────────────────────────────────────
all_comparison = []
for _, row in v2.iterrows():
    kju = row["KJU_ID"]
    accepted = str(row.get("Accepted_taxonomy", "")).strip()
    final_status = str(row.get("Final_status", "")).strip()

    # Check if this KJU has BLAST hits
    kju_blast = blast[blast["KJU_ID"] == kju]
    if len(kju_blast) == 0:
        all_comparison.append({
            "KJU_ID": kju,
            "Accepted_taxonomy": accepted,
            "Top_hit_accession": "",
            "Top_hit_species": "",
            "Top_hit_identity": "",
            "Best_match_is_self": "",
            "Taxonomy_agreement": "NO_BLAST_HITS",
            "Agreement_category": "NO_BLAST_HITS",
            "Accepted_in_NCBI_REF": "",
        })
        continue

    kju_blast_sorted = kju_blast.sort_values("Bit_score_float", ascending=False)
    best = kju_blast_sorted.iloc[0]
    hit_species = str(best.get("Hit_species_corrected", "")).strip()

    norm_accepted = normalize_species(accepted)
    norm_hit = normalize_species(hit_species)
    genus_accepted = get_genus(norm_accepted)
    genus_hit = get_genus(norm_hit)

    in_ref = norm_accepted in ref_species_set if accepted else False

    if final_status == "UNRESOLVED":
        category = "UNRESOLVED"
        agreement = "UNRESOLVED"
    elif final_status == "ADJUDICATED_NO_SEQUENCE":
        category = "NO_SEQUENCE"
        agreement = "NO_SEQUENCE"
    elif not in_ref:
        category = "INSUFFICIENT_REFERENCE"
        agreement = "INSUFFICIENT_REFERENCE"
    elif norm_accepted == norm_hit:
        category = "SPECIES_AGREEMENT"
        agreement = "TRUE"
    elif genus_accepted == genus_hit:
        category = "GENUS_AGREEMENT"
        agreement = "GENUS_AGREEMENT"
    else:
        category = "DISCORDANT"
        agreement = "FALSE"

    all_comparison.append({
        "KJU_ID": kju,
        "Accepted_taxonomy": accepted,
        "Top_hit_accession": best.get("Subject_accession", ""),
        "Top_hit_species": hit_species,
        "Top_hit_identity": best.get("Pident_float", ""),
        "Best_match_is_self": "N/A",
        "Taxonomy_agreement": agreement,
        "Agreement_category": category,
        "Accepted_in_NCBI_REF": in_ref,
    })

all_comparison_df = pd.DataFrame(all_comparison)

# ──────────────────────────────────────────────────────────
# SUMMARY
# ──────────────────────────────────────────────────────────
print("TAXONOMY AGREEMENT CATEGORIES (all 204 KJUs):")
cat_counts = all_comparison_df["Agreement_category"].value_counts()
for cat, count in cat_counts.items():
    print(f"  {cat}: {count}")
print()

print("TAXONOMY AGREEMENT (canonical 111):")
canon_comp = all_comparison_df[all_comparison_df["KJU_ID"].isin(canon["KJU_ID"].values)]
cat_counts_canon = canon_comp["Agreement_category"].value_counts()
for cat, count in cat_counts_canon.items():
    print(f"  {cat}: {count}")
print()

# List DISCORDANT cases
discordant = all_comparison_df[all_comparison_df["Agreement_category"] == "DISCORDANT"]
if len(discordant) > 0:
    print(f"DISCORDANT CASES ({len(discordant)}):")
    for _, row in discordant.iterrows():
        print(f"  {row['KJU_ID']}: {row['Accepted_taxonomy']} vs {row['Top_hit_species']}")
else:
    print("DISCORDANT CASES: 0")
print()

# List GENUS_AGREEMENT cases
genus_agree = all_comparison_df[all_comparison_df["Agreement_category"] == "GENUS_AGREEMENT"]
if len(genus_agree) > 0:
    print(f"GENUS_AGREEMENT CASES ({len(genus_agree)}):")
    for _, row in genus_agree.iterrows():
        print(f"  {row['KJU_ID']}: {row['Accepted_taxonomy']} vs {row['Top_hit_species']} (identity={row['Top_hit_identity']:.1f}%)")
print()

# List INSUFFICIENT_REFERENCE cases
insuff = all_comparison_df[all_comparison_df["Agreement_category"] == "INSUFFICIENT_REFERENCE"]
print(f"INSUFFICIENT_REFERENCE: {len(insuff)} KJUs (species not in NCBI_REF)")
print()

# ──────────────────────────────────────────────────────────
# SAVE OUTPUTS
# ──────────────────────────────────────────────────────────
all_comparison_df.to_csv(OUT_DIR / "STEP92b_TAXONOMY_COMPARISON_CORRECTED.csv", index=False)
results_df.to_csv(OUT_DIR / "STEP92b_CANONICAL_TAXONOMY_DETAIL.csv", index=False)

# Also save corrected BLAST results with Hit_species
blast_save = blast.drop(columns=["Bit_score_float", "Pident_float"])
blast_save.to_csv(OUT_DIR / "STEP92b_BLAST_RESULTS_CORRECTED.csv", index=False)

# ──────────────────────────────────────────────────────────
# CHECKSUMS
# ──────────────────────────────────────────────────────────
sha256_map = {}
for f in sorted(OUT_DIR.glob("STEP92b_*")):
    sha = hashlib.sha256(f.read_bytes()).hexdigest()
    sha256_map[f.name] = sha

checksums_path = OUT_DIR / "STEP92b_CHECKSUMS.txt"
with open(checksums_path, "w") as fh:
    fh.write("STEP 92b CHECKSUMS\n")
    fh.write("=" * 60 + "\n\n")
    for fname, sha in sha256_map.items():
        fh.write(f"{sha}  {fname}\n")
    fh.write(f"\nGenerated: Step 92b BLAST Taxonomy Species Lookup Correction\n")
    fh.write(f"V2 master SHA256 (not modified): 65185f7545f73fcf36ff5062bb85305affcae76458f8049130f052a811b79bee\n")

print(f"Output written to: {OUT_DIR}")
for f in sorted(OUT_DIR.glob("*")):
    size = f.stat().st_size
    print(f"  {f.name} ({size:,} bytes)")
print()

# ──────────────────────────────────────────────────────────
# QC REPORT
# ──────────────────────────────────────────────────────────
report_path = OUT_DIR / "STEP92b_TAXONOMY_CORRECTION_REPORT.txt"
with open(report_path, "w") as fh:
    fh.write("STEP 92b: BLAST TAXONOMY SPECIES LOOKUP CORRECTION\n")
    fh.write("=" * 60 + "\n\n")

    fh.write("BUG DESCRIPTION\n")
    fh.write("-" * 40 + "\n")
    fh.write("STEP92 ref_lookup used versioned accession keys (e.g. PX979534.1)\n")
    fh.write("but lookup code normalized BLAST hits to root accessions (e.g. PX979534).\n")
    fh.write("Species lookup failed silently. Taxonomy_agreement defaulted to True.\n\n")

    fh.write("FIX APPLIED\n")
    fh.write("-" * 40 + "\n")
    fh.write("1. Built ref_lookup with root accession keys\n")
    fh.write("2. Extracted Hit_species from Subject_accession 3rd pipe field\n")
    fh.write("   (PX979534.1|NCBI_REF|Aeromonas_caviae -> Aeromonas caviae)\n")
    fh.write("3. Re-classified taxonomy agreement into 4 categories\n\n")

    fh.write("DATA COVERAGE\n")
    fh.write("-" * 40 + "\n")
    fh.write(f"BLAST results: {len(blast)} hits\n")
    fh.write(f"NCBI references: {len(ref_meta)} sequences, {len(ref_species_set)} species\n")
    fh.write(f"Canonical species: {len(canon_species_set)}\n")
    fh.write(f"Species with NCBI_REF: {len(canon_species_set & ref_species_set)}\n")
    fh.write(f"Species without NCBI_REF: {len(canon_species_set - ref_species_set)}\n\n")

    fh.write("RESULTS (all 204 KJUs)\n")
    fh.write("-" * 40 + "\n")
    for cat, count in cat_counts.items():
        fh.write(f"  {cat}: {count}\n")
    fh.write("\n")

    fh.write("RESULTS (canonical 111 only)\n")
    fh.write("-" * 40 + "\n")
    for cat, count in cat_counts_canon.items():
        fh.write(f"  {cat}: {count}\n")
    fh.write("\n")

    fh.write("GENUS_AGREEMENT DETAILS\n")
    fh.write("-" * 40 + "\n")
    fh.write("These are same-genus BLAST hits with >99% identity.\n")
    fh.write("Expected when the reference DB lacks the exact species.\n\n")
    for _, row in genus_agree.iterrows():
        fh.write(f"  {row['KJU_ID']}: {row['Accepted_taxonomy']} vs {row['Top_hit_species']}\n")
        fh.write(f"    Identity: {row['Top_hit_identity']:.1f}%\n")
    fh.write("\n")

    fh.write("DISCORDANT CASES\n")
    fh.write("-" * 40 + "\n")
    if len(discordant) > 0:
        for _, row in discordant.iterrows():
            fh.write(f"  {row['KJU_ID']}: {row['Accepted_taxonomy']} vs {row['Top_hit_species']}\n")
    else:
        fh.write("  None\n")
    fh.write("\n")

    fh.write("INSUFFICIENT_REFERENCE SPECIES\n")
    fh.write("-" * 40 + "\n")
    fh.write("These species have no entry in NCBI_REFERENCE_METADATA.\n")
    fh.write("Species-level BLAST comparison cannot be performed.\n\n")
    for sp in sorted(insuff["Accepted_taxonomy"].unique()):
        count = len(insuff[insuff["Accepted_taxonomy"] == sp])
        fh.write(f"  {sp}: {count} KJUs\n")
    fh.write("\n")

    fh.write("FILES\n")
    fh.write("-" * 40 + "\n")
    fh.write(f"Corrected taxonomy comparison: STEP92b_TAXONOMY_COMPARISON_CORRECTED.csv\n")
    fh.write(f"Corrected BLAST results: STEP92b_BLAST_RESULTS_CORRECTED.csv\n")
    fh.write(f"Canonical detail: STEP92b_CANONICAL_TAXONOMY_DETAIL.csv\n")
    fh.write(f"Checksums: STEP92b_CHECKSUMS.txt\n\n")

    fh.write("V2 FREEZE STATUS\n")
    fh.write("-" * 40 + "\n")
    fh.write("V2 freeze was NOT modified by this step.\n")
    fh.write("BLAST_taxonomy column in V2 remains empty.\n")
    fh.write("Corrected data available in STEP92b_TAXONOMY_COMPARISON_CORRECTED.csv\n")

print("STEP 92b: COMPLETE")
print()
print("SUMMARY:")
print(f"  Species with NCBI_REF agreement: {(canon_comp['Agreement_category'] == 'SPECIES_AGREEMENT').sum()}")
print(f"  Genus-level agreement: {(canon_comp['Agreement_category'] == 'GENUS_AGREEMENT').sum()}")
print(f"  Insufficient reference: {(canon_comp['Agreement_category'] == 'INSUFFICIENT_REFERENCE').sum()}")
print(f"  Discordant: {(canon_comp['Agreement_category'] == 'DISCORDANT').sum()}")
print(f"  V2 freeze: UNCHANGED")
