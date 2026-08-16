#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib
import re

BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP83 = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"
STEP87 = ROOT / "STEP87_DETERMINISTIC_CANONICAL_SEQUENCE"
STEP88 = ROOT / "STEP88_CANONICAL_SEQUENCE_QC"
STEP89 = ROOT / "STEP89_CANONICAL_FASTA"
STEP90 = ROOT / "STEP90_ACCESSION_TAXONOMY_AUDIT"
STEP91 = ROOT / "STEP91_SEQUENCE_TYPE_VALIDATION"
STEP92 = ROOT / "STEP92_BLAST_VALIDATION"
STEP93 = ROOT / "STEP93_PHYLOGENETIC_VALIDATION"
STEP94 = ROOT / "STEP94_BIOLOGICAL_IDENTITY_ADJUDICATION"
OUT = ROOT / "STEP95_FINAL_INTEGRATED_MASTER_V2"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
STEP87_AUDIT = STEP87 / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"
STEP94_IDENTITY = STEP94 / "STEP94_BIOLOGICAL_IDENTITY_FINAL.csv"
STEP94_EVIDENCE = STEP94 / "STEP94_EVIDENCE_ASSEMBLY.csv"
STEP88_PROV = STEP88 / "STEP88_PROVENANCE_AUDIT.csv"
STEP92_TAXCOMP = STEP92 / "STEP92_TAXONOMY_COMPARISON.csv"
STEP91_TYPE = STEP91 / "STEP91_SEQUENCE_TYPE_AUDIT.csv"
STEP93_CLUST = STEP93 / "STEP93_CLUSTERING_AUDIT.csv"

V2_OUT = OUT / "STEP95_INTEGRATED_MASTER_V2.csv"
SUMMARY = OUT / "STEP95_SUMMARY.txt"
FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

print("=" * 110)
print("STEP 95: FINAL INTEGRATED MASTER CANDIDATE V2")
print("=" * 110)

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master SHA256: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
audit87 = pd.read_csv(STEP87_AUDIT, dtype=str, low_memory=False)
identity = pd.read_csv(STEP94_IDENTITY, dtype=str, low_memory=False) if STEP94_IDENTITY.exists() else pd.DataFrame()
evidence = pd.read_csv(STEP94_EVIDENCE, dtype=str, low_memory=False) if STEP94_EVIDENCE.exists() else pd.DataFrame()
provenance = pd.read_csv(STEP88_PROV, dtype=str, low_memory=False) if STEP88_PROV.exists() else pd.DataFrame()
taxcomp = pd.read_csv(STEP92_TAXCOMP, dtype=str, low_memory=False) if STEP92_TAXCOMP.exists() else pd.DataFrame()
typeaudit = pd.read_csv(STEP91_TYPE, dtype=str, low_memory=False) if STEP91_TYPE.exists() else pd.DataFrame()
clustering = pd.read_csv(STEP93_CLUST, dtype=str, low_memory=False) if STEP93_CLUST.exists() else pd.DataFrame()

print(f"Master rows: {len(master)}")
assert len(master) == 204

def safe_lookup(df, kju, col):
    if df.empty or "KJU_ID" not in df.columns:
        return ""
    rows = df[df["KJU_ID"] == kju]
    if len(rows) == 0:
        return ""
    val = rows.iloc[0].get(col, "")
    return str(val) if pd.notna(val) else ""

v2_rows = []
for _, mrow in master.iterrows():
    kju = mrow["KJU_ID"]
    accepted_acc = str(mrow.get("Accepted_accession", "")).strip()
    accepted_tax = str(mrow.get("Accepted_taxonomy", "")).strip()
    accepted_len = str(mrow.get("Accepted_length", "")).strip()
    accepted_seq_type = str(mrow.get("Accepted_sequence_type", "")).strip()
    bio_id = str(mrow.get("Biological_identity", "")).strip()
    reconcile = str(mrow.get("Reconciliation_status", "")).strip()
    adj_type = str(mrow.get("Adjudication_type", "")).strip()
    confidence = str(mrow.get("Confidence", "")).strip()
    adj_reason = str(mrow.get("Adjudication_reason", "")).strip()
    adj_evidence_str = str(mrow.get("Adjudication_evidence", "")).strip()
    tax_evidence_source = str(mrow.get("Taxonomic_evidence_source", "")).strip()
    norm_acc = str(mrow.get("Normalized_accession", "")).strip()

    a = audit87[audit87["KJU_ID"] == kju]
    sel_status = str(a.iloc[0]["Selection_status"]) if len(a) > 0 else ""
    sel_evidence = str(a.iloc[0]["Selection_evidence"]) if len(a) > 0 else ""
    canon_sha = str(a.iloc[0].get("Canonical_sequence_SHA256", "")) if len(a) > 0 else ""
    canon_len = str(a.iloc[0].get("Canonical_sequence_length", "")) if len(a) > 0 else ""
    canon_acc = str(a.iloc[0].get("Canonical_accession", "")) if len(a) > 0 else ""
    canon_seq = str(a.iloc[0].get("Canonical_sequence", "")) if len(a) > 0 else ""
    n_cand = str(a.iloc[0].get("Unique_candidate_count", "")) if len(a) > 0 else ""

    prov_row = provenance[provenance["KJU_ID"] == kju] if not provenance.empty else pd.DataFrame()
    src_file = str(prov_row.iloc[0].get("Source_files", "")) if len(prov_row) > 0 else ""
    src_record = str(prov_row.iloc[0].get("Source_record_IDs", "")) if len(prov_row) > 0 else ""

    final_id_row = identity[identity["KJU_ID"] == kju] if not identity.empty else pd.DataFrame()
    final_bio = str(final_id_row.iloc[0].get("Final_biological_identity", "")) if len(final_id_row) > 0 else bio_id
    final_status = str(final_id_row.iloc[0].get("Final_status", "")) if len(final_id_row) > 0 else ""

    blast_tax = safe_lookup(taxcomp, kju, "Top_hit_species")
    blast_agree = safe_lookup(taxcomp, kju, "Taxonomy_agreement")
    seq_type_valid = safe_lookup(typeaudit, kju, "Length_in_expected_range")
    in_phylo = safe_lookup(clustering, kju, "In_tree")

    ev_row = evidence[evidence["KJU_ID"] == kju] if not evidence.empty else pd.DataFrame()
    ev_conflicts = str(ev_row.iloc[0].get("Evidence_conflicts", "")) if len(ev_row) > 0 else ""

    v2_rows.append({
        "KJU_ID": kju,
        "Accepted_accession": accepted_acc,
        "Accepted_taxonomy": accepted_tax,
        "Accepted_length": accepted_len,
        "Accepted_sequence_type": accepted_seq_type,
        "Biological_identity": bio_id,
        "Reconciliation_status": reconcile,
        "Adjudication_type": adj_type,
        "Confidence": confidence,
        "Adjudication_reason": adj_reason,
        "Adjudication_evidence": adj_evidence_str,
        "Taxonomic_evidence_source": tax_evidence_source,
        "Normalized_accession": norm_acc,
        "Canonical_sequence": canon_seq,
        "Canonical_sequence_SHA256": canon_sha,
        "Canonical_sequence_length": canon_len,
        "Canonical_accession": canon_acc,
        "Canonical_sequence_source_files": src_file,
        "Canonical_sequence_source_record_IDs": src_record,
        "Selection_status": sel_status,
        "Selection_evidence": sel_evidence,
        "Unique_candidate_count": n_cand,
        "Sequence_type_validated": seq_type_valid,
        "BLAST_taxonomy": blast_tax,
        "BLAST_taxonomy_agreement": blast_agree,
        "In_phylogenetic_tree": in_phylo,
        "Evidence_conflicts": ev_conflicts,
        "Final_biological_identity": final_bio,
        "Final_status": final_status,
    })

v2 = pd.DataFrame(v2_rows)
v2.to_csv(V2_OUT, index=False)

print(f"\nV2 rows: {len(v2)}")
print(f"V2 unique KJUs: {v2['KJU_ID'].nunique()}")
assert len(v2) == 204
assert v2["KJU_ID"].nunique() == 204

status_counts = v2["Final_status"].value_counts()
print(f"\nFinal status distribution:")
for s, c in status_counts.items():
    print(f"  {s}: {c}")

has_canon_seq = (v2["Canonical_sequence"].apply(lambda x: str(x).strip() != "" and str(x).strip().lower() != "nan")).sum()
print(f"\nKJUs with canonical sequence: {has_canon_seq}")
print(f"KJUs without canonical sequence: {204 - has_canon_seq}")

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(SUMMARY, "w") as f:
    f.write("STEP 95: FINAL INTEGRATED MASTER CANDIDATE V2\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS\n\n")
    f.write(f"Total rows: {len(v2)}\n")
    f.write(f"Unique KJUs: {v2['KJU_ID'].nunique()}\n\n")
    f.write("Final status distribution:\n")
    for s, c in status_counts.items():
        f.write(f"  {s}: {c}\n")
    f.write(f"\nWith canonical sequence: {has_canon_seq}\n")
    f.write(f"Without canonical sequence: {204 - has_canon_seq}\n")
    f.write(f"\nV2 output: {V2_OUT}\n")
    f.write(f"\nFrozen V1 master untouched.\n")

print("\n" + "=" * 110)
print("STEP 95 COMPLETE")
print("=" * 110)
