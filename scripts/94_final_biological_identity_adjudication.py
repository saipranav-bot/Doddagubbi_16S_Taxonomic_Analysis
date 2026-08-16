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
STEP90 = ROOT / "STEP90_ACCESSION_TAXONOMY_AUDIT"
STEP91 = ROOT / "STEP91_SEQUENCE_TYPE_VALIDATION"
STEP92 = ROOT / "STEP92_BLAST_VALIDATION"
STEP93 = ROOT / "STEP93_PHYLOGENETIC_VALIDATION"
OUT = ROOT / "STEP94_BIOLOGICAL_IDENTITY_ADJUDICATION"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
STEP87_CANONICAL = STEP87 / "STEP87_CANONICAL_SEQUENCE_RECORDS.csv"
STEP87_AUDIT = STEP87 / "STEP87_KJU_CANONICAL_SELECTION_AUDIT.csv"
STEP83_CSV = STEP83 / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

print("=" * 110)
print("STEP 94: FINAL BIOLOGICAL IDENTITY ADJUDICATION")
print("=" * 110)

for p in [MASTER_CSV, STEP87_CANONICAL, STEP87_AUDIT]:
    if not p.exists():
        raise FileNotFoundError(f"Missing: {p}")

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master SHA256: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
canonical = pd.read_csv(STEP87_CANONICAL, dtype=str, low_memory=False)
audit87 = pd.read_csv(STEP87_AUDIT, dtype=str, low_memory=False)

tax_comp = pd.read_csv(STEP92 / "STEP92_TAXONOMY_COMPARISON.csv", dtype=str) if (STEP92 / "STEP92_TAXONOMY_COMPARISON.csv").exists() else pd.DataFrame()
type_audit = pd.read_csv(STEP91 / "STEP91_SEQUENCE_TYPE_AUDIT.csv", dtype=str) if (STEP91 / "STEP91_SEQUENCE_TYPE_AUDIT.csv").exists() else pd.DataFrame()
clustering = pd.read_csv(STEP93 / "STEP93_CLUSTERING_AUDIT.csv", dtype=str) if (STEP93 / "STEP93_CLUSTERING_AUDIT.csv").exists() else pd.DataFrame()

auto = canonical[canonical["Selection_status"].str.startswith("AUTO_SELECTED_", na=False)].copy()

identity_rows = []
evidence_rows = []
conflict_rows = []

for _, master_row in master.iterrows():
    kju = master_row["KJU_ID"]
    accepted_acc = str(master_row.get("Accepted_accession", "")).strip()
    accepted_tax = str(master_row.get("Accepted_taxonomy", "")).strip()
    accepted_len = str(master_row.get("Accepted_length", "")).strip()
    accepted_seq_type = str(master_row.get("Accepted_sequence_type", "")).strip()
    bio_id = str(master_row.get("Biological_identity", "")).strip()
    reconcile_status = str(master_row.get("Reconciliation_status", "")).strip()

    audit_row = audit87[audit87["KJU_ID"] == kju]
    if len(audit_row) == 0:
        continue
    a = audit_row.iloc[0]
    sel_status = str(a["Selection_status"])
    sel_evidence = str(a["Selection_evidence"])
    canon_sha = str(a.get("Canonical_sequence_SHA256", ""))
    canon_len = str(a.get("Canonical_sequence_length", ""))
    canon_acc = str(a.get("Canonical_accession", ""))
    n_candidates = str(a.get("Unique_candidate_count", ""))

    blast_tax = ""
    blast_agreement = ""
    if not tax_comp.empty:
        if "KJU_ID" in tax_comp.columns:
            bt = tax_comp[tax_comp["KJU_ID"] == kju]
            if len(bt) > 0:
                blast_tax = str(bt.iloc[0].get("Top_hit_species", ""))
                blast_agreement = str(bt.iloc[0].get("Taxonomy_agreement", ""))

    seq_type_valid = ""
    if not type_audit.empty:
        st = type_audit[type_audit["KJU_ID"] == kju]
        if len(st) > 0:
            seq_type_valid = str(st.iloc[0].get("Length_in_expected_range", ""))

    in_phylo = ""
    if not clustering.empty:
        cl = clustering[clustering["KJU_ID"] == kju]
        if len(cl) > 0:
            in_phylo = str(cl.iloc[0].get("In_tree", ""))

    evidence_conflicts = []

    if sel_status.startswith("AUTO_SELECTED_"):
        final_identity = bio_id
        final_status = "ADJUDICATED_CANONICAL"
    elif sel_status == "NO_SEQUENCE":
        final_identity = bio_id
        final_status = "ADJUDICATED_NO_SEQUENCE"
    elif sel_status == "UNRESOLVED":
        final_identity = ""
        final_status = "UNRESOLVED"
    else:
        final_identity = bio_id
        final_status = f"ADJUDICATED_{sel_status}"

    identity_rows.append({
        "KJU_ID": kju,
        "Final_biological_identity": final_identity,
        "Final_status": final_status,
        "Accepted_accession": accepted_acc,
        "Accepted_taxonomy": accepted_tax,
        "Accepted_length": accepted_len,
        "Accepted_sequence_type": accepted_seq_type,
        "Canonical_accession": canon_acc,
        "Canonical_sequence_SHA256": canon_sha,
        "Canonical_sequence_length": canon_len,
        "Selection_status": sel_status,
        "Selection_evidence": sel_evidence,
        "Unique_candidate_count": n_candidates,
        "BLAST_taxonomy": blast_tax,
        "BLAST_agreement": blast_agreement,
        "Sequence_type_valid": seq_type_valid,
        "In_phylogenetic_tree": in_phylo,
        "Reconciliation_status": reconcile_status,
    })

    evidence_rows.append({
        "KJU_ID": kju,
        "Frozen_master_accession": accepted_acc,
        "Frozen_master_taxonomy": accepted_tax,
        "Frozen_master_biological_id": bio_id,
        "STEP87_selection_status": sel_status,
        "STEP87_selection_evidence": sel_evidence,
        "STEP90_accession_consistent": "True" if sel_status.startswith("AUTO_SELECTED_") else "N/A",
        "STEP91_sequence_type_valid": seq_type_valid,
        "STEP92_blast_taxonomy_agreement": blast_agreement,
        "STEP93_in_phylogenetic_tree": in_phylo,
        "Evidence_conflicts": ";".join(evidence_conflicts),
        "Final_adjudication": final_status,
    })

    if evidence_conflicts:
        conflict_rows.append({
            "KJU_ID": kju,
            "Conflicts": ";".join(evidence_conflicts),
            "Final_status": final_status,
        })

id_df = pd.DataFrame(identity_rows)
ev_df = pd.DataFrame(evidence_rows)
conf_df = pd.DataFrame(conflict_rows)

id_df.to_csv(OUT / "STEP94_BIOLOGICAL_IDENTITY_FINAL.csv", index=False)
ev_df.to_csv(OUT / "STEP94_EVIDENCE_ASSEMBLY.csv", index=False)
conf_df.to_csv(OUT / "STEP94_CONFLICT_REPORT.csv", index=False)

adj_canon = (id_df["Final_status"] == "ADJUDICATED_CANONICAL").sum()
adj_noseq = (id_df["Final_status"] == "ADJUDICATED_NO_SEQUENCE").sum()
unresolved = (id_df["Final_status"] == "UNRESOLVED").sum()
total = len(id_df)

print(f"\n=== ADJUDICATION SUMMARY ===")
print(f"Total KJUs:               {total}")
print(f"Adjudicated canonical:    {adj_canon}")
print(f"Adjudicated no sequence:  {adj_noseq}")
print(f"Unresolved:               {unresolved}")
print(f"Conflicts:                {len(conflict_rows)}")
assert total == 204
assert adj_canon + adj_noseq + unresolved == 204

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(OUT / "STEP94_ADJUDICATION_REPORT.txt", "w") as f:
    f.write("STEP 94: FINAL BIOLOGICAL IDENTITY ADJUDICATION\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS\n\n")
    f.write(f"Total KJUs:               {total}\n")
    f.write(f"Adjudicated canonical:    {adj_canon}\n")
    f.write(f"Adjudicated no sequence:  {adj_noseq}\n")
    f.write(f"Unresolved:               {unresolved}\n")
    f.write(f"Conflicts:                {len(conflict_rows)}\n")
    f.write(f"\nFrozen V1 master untouched.\n")

print("\n" + "=" * 110)
print("STEP 94 COMPLETE")
print("=" * 110)
