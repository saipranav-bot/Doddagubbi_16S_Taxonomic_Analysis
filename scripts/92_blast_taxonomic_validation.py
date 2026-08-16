#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re
import subprocess
import tempfile
import os


# ============================================================
# STEP 92
# BLAST-BASED TAXONOMIC VALIDATION
#
# PURPOSE:
#   Local BLAST of canonical sequences against NCBI reference
#   sequences. Compare top-hit taxonomy with master taxonomy.
# ============================================================


BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP87 = ROOT / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
STEP89 = ROOT / "STEP89_CANONICAL_FASTA"
NCBI_REF = BASE / "07_Master_Metadata" / "FINAL_109_PHYLOGENY_V1" / "NCBI_SPECIES_REFERENCES_V1"
OUT = ROOT / "STEP92_BLAST_VALIDATION"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
STEP87_AUDIT = STEP87 / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"
REF_FASTA = NCBI_REF / "NCBI_REFERENCE_16S.fasta"
REF_META = NCBI_REF / "NCBI_REFERENCE_METADATA.csv"
CANONICAL_FASTA = STEP89 / "STEP89_CANONICAL_KJU_16S.fasta"

FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

BLAST_DB = OUT / "STEP92_REF_DB"
BLAST_OUT = OUT / "STEP92_BLAST_RESULTS.csv"
TAX_COMP = OUT / "STEP92_TAXONOMY_COMPARISON.csv"
DISCORDANT = OUT / "STEP92_DISCORDANT_KJUS.csv"
REPORT = OUT / "STEP92_BLAST_QC_REPORT.txt"


def normalize_species(value):
    if pd.isna(value):
        return ""
    value = str(value).strip().lower()
    value = re.sub(r"\s+", " ", value)
    return value


print("=" * 110)
print("STEP 92: BLAST-BASED TAXONOMIC VALIDATION")
print("=" * 110)

for p in [MASTER_CSV, STEP87_CANONICAL, STEP87_AUDIT, REF_FASTA, REF_META, CANONICAL_FASTA]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master SHA256: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
canonical = pd.read_csv(STEP87_CANONICAL, dtype=str, low_memory=False)
audit87 = pd.read_csv(STEP87_AUDIT, dtype=str, low_memory=False)
ref_meta = pd.read_csv(REF_META, dtype=str, low_memory=False)

auto = canonical[canonical["Selection_status"].str.startswith("AUTO_SELECTED_", na=False)].copy()

print(f"Canonical queries: {len(auto)}")
print(f"Reference seqs:    {len(ref_meta)}")


ref_lookup = {}
for _, r in ref_meta.iterrows():
    acc = str(r.get("accession", "")).strip()
    ref_lookup[acc] = {
        "species": str(r.get("species", "")).strip(),
        "organism": str(r.get("organism", "")).strip(),
    }


print("\nBuilding local BLAST database...")
db_cmd = [
    "makeblastdb",
    "-in", str(REF_FASTA),
    "-dbtype", "nucl",
    "-out", str(BLAST_DB),
    "-parse_seqids",
]
result = subprocess.run(db_cmd, capture_output=True, text=True, timeout=120)
if result.returncode != 0:
    print(f"makeblastdb stderr: {result.stderr[:500]}")
    raise RuntimeError(f"makeblastdb failed: {result.returncode}")
print("BLAST database built.")

print("Running BLASTN...")
blast_cmd = [
    "blastn",
    "-query", str(CANONICAL_FASTA),
    "-db", str(BLAST_DB),
    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
    "-max_target_seqs", "5",
    "-evalue", "1e-50",
    "-num_threads", "4",
]
result = subprocess.run(blast_cmd, capture_output=True, text=True, timeout=300)
if result.returncode != 0:
    print(f"blastn stderr: {result.stderr[:500]}")
    raise RuntimeError(f"blastn failed: {result.returncode}")

blast_lines = result.stdout.strip().split("\n")
print(f"BLAST hits: {len(blast_lines)}")

blast_rows = []
for line in blast_lines:
    if not line.strip():
        continue
    parts = line.split("\t")
    if len(parts) < 12:
        continue
    qseqid = parts[0]
    sseqid = parts[1]
    pident = parts[2]
    length = parts[3]
    evalue = parts[10]
    bitscore = parts[11]
    stitle = parts[12] if len(parts) > 12 else ""

    kju = qseqid.split("|")[0] if "|" in qseqid else qseqid

    ref_acc = sseqid.split("|")[0] if "|" in sseqid else sseqid
    ref_acc = re.sub(r"\.\d+$", "", ref_acc.upper())

    ref_info = ref_lookup.get(ref_acc, {})
    hit_species = ref_info.get("species", "")
    hit_organism = ref_info.get("organism", "")

    blast_rows.append({
        "Query_ID": qseqid,
        "KJU_ID": kju,
        "Subject_accession": sseqid,
        "Subject_accession_root": ref_acc,
        "Percent_identity": pident,
        "Alignment_length": length,
        "E_value": evalue,
        "Bit_score": bitscore,
        "Hit_species": hit_species,
        "Hit_organism": hit_organism,
        "Subject_title": stitle,
    })

blast_df = pd.DataFrame(blast_rows)
blast_df.to_csv(BLAST_OUT, index=False)

print(f"BLAST results saved: {BLAST_OUT}")

taxonomy_rows = []
discordant_rows = []

for _, crow in auto.iterrows():
    kju = crow["KJU_ID"]
    accepted_acc = str(crow["Accepted_accession"]).strip()
    accepted_acc_root = re.sub(r"\.\d+$", "", accepted_acc.upper())

    master_row = master[master["KJU_ID"] == kju].iloc[0]
    accepted_species = normalize_species(master_row.get("Accepted_taxonomy", ""))

    kju_blast = blast_df[blast_df["KJU_ID"] == kju].copy()
    kju_blast["Percent_identity"] = kju_blast["Percent_identity"].astype(float)

    if kju_blast.empty:
        taxonomy_rows.append({
            "KJU_ID": kju,
            "Accepted_taxonomy": str(master_row.get("Accepted_taxonomy", "")),
            "Top_hit_accession": "",
            "Top_hit_species": "",
            "Top_hit_identity": "",
            "Best_match_is_self": "N/A",
            "Taxonomy_agreement": "NO_BLAST_HITS",
        })
        continue

    best = kju_blast.sort_values("Bit_score", ascending=False).iloc[0]
    top_species = normalize_species(best["Hit_species"])
    top_identity = float(best["Percent_identity"])

    self_hit = accepted_acc_root in str(best["Subject_accession_root"])

    tax_agree = True
    if not self_hit and top_species and accepted_species:
        tax_agree = normalize_species(accepted_species) == top_species

    taxonomy_rows.append({
        "KJU_ID": kju,
        "Accepted_taxonomy": str(master_row.get("Accepted_taxonomy", "")),
        "Top_hit_accession": str(best["Subject_accession"]),
        "Top_hit_species": str(best["Hit_species"]),
        "Top_hit_identity": str(top_identity),
        "Best_match_is_self": str(self_hit),
        "Taxonomy_agreement": str(tax_agree),
    })

    if not tax_agree and not self_hit:
        discordant_rows.append({
            "KJU_ID": kju,
            "Accepted_taxonomy": str(master_row.get("Accepted_taxonomy", "")),
            "Top_hit_accession": str(best["Subject_accession"]),
            "Top_hit_species": str(best["Hit_species"]),
            "Top_hit_identity": str(top_identity),
            "Accepted_accession": accepted_acc,
        })

unresolved = audit87[audit87["Selection_status"] == "UNRESOLVED"]
for _, urow in unresolved.iterrows():
    taxonomy_rows.append({
        "KJU_ID": urow["KJU_ID"],
        "Accepted_taxonomy": str(urow.get("Accepted_taxonomy", "")),
        "Top_hit_accession": "",
        "Top_hit_species": "",
        "Top_hit_identity": "",
        "Best_match_is_self": "N/A",
        "Taxonomy_agreement": "UNRESOLVED",
    })

no_seq = audit87[audit87["Selection_status"] == "NO_SEQUENCE"]
for _, nsrow in no_seq.iterrows():
    taxonomy_rows.append({
        "KJU_ID": nsrow["KJU_ID"],
        "Accepted_taxonomy": str(nsrow.get("Accepted_taxonomy", "")),
        "Top_hit_accession": "",
        "Top_hit_species": "",
        "Top_hit_identity": "",
        "Best_match_is_self": "N/A",
        "Taxonomy_agreement": "NO_SEQUENCE",
    })

tax_df = pd.DataFrame(taxonomy_rows)
tax_df.to_csv(TAX_COMP, index=False)

disc_df = pd.DataFrame(discordant_rows)
disc_df.to_csv(DISCORDANT, index=False)

agree_count = int(tax_df[tax_df["Taxonomy_agreement"] == "True"].shape[0])
disagree_count = len(discordant_rows)
self_match = int(tax_df[tax_df["Best_match_is_self"] == "True"].shape[0])
no_hit = int(tax_df[tax_df["Taxonomy_agreement"] == "NO_BLAST_HITS"].shape[0])

print(f"\n=== BLAST TAXONOMY COMPARISON ===")
print(f"Total canonical: {len(auto)}")
print(f"Self-matching (best hit = own accession): {self_match}")
print(f"Taxonomy agreement: {agree_count}")
print(f"Taxonomy discordance: {disagree_count}")
print(f"No BLAST hits: {no_hit}")

if discordant_rows:
    print(f"\nDiscordant KJUs:")
    for d in discordant_rows:
        print(f"  {d['KJU_ID']}: accepted={d['Accepted_taxonomy']}, blast={d['Top_hit_species']} ({d['Top_hit_identity']}%)")

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(REPORT, "w") as f:
    f.write("STEP 92: BLAST-BASED TAXONOMIC VALIDATION\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS\n\n")
    f.write(f"Reference sequences: {len(ref_meta)}\n")
    f.write(f"Canonical queries:   {len(auto)}\n\n")
    f.write(f"BLAST HITS:          {len(blast_df)}\n")
    f.write(f"Self-matching:       {self_match}\n")
    f.write(f"Taxonomy agreement:  {agree_count}\n")
    f.write(f"Taxonomy discordance:{disagree_count}\n")
    f.write(f"No BLAST hits:       {no_hit}\n\n")
    if discordant_rows:
        f.write("DISCORDANT KJUs:\n")
        for d in discordant_rows:
            f.write(f"  {d['KJU_ID']}: accepted={d['Accepted_taxonomy']}, blast={d['Top_hit_species']} ({d['Top_hit_identity']}%)\n")
    f.write(f"\nFrozen V1 master untouched.\n")

print("\n" + "=" * 110)
print("STEP 92 COMPLETE")
print("=" * 110)
