#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import hashlib

BASE = Path.home() / "Doddagubbi_Master"
ROOT = BASE / "github_results" / "reconciliation_111_vs_109" / "final_integrated_master"
FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
STEP83 = ROOT / "STEP83_TRUE_SEQUENCE_INVENTORY"
STEP95 = ROOT / "STEP95_FINAL_INTEGRATED_MASTER_V2"
OUT = ROOT / "STEP97_BIOLOGICAL_UNIQUENESS"
OUT.mkdir(parents=True, exist_ok=True)

MASTER_CSV = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
V2_CSV = STEP95 / "STEP95_INTEGRATED_MASTER_V2.csv"
STEP83_CSV = STEP83 / "STEP83_TRUE_SEQUENCE_RECORDS.csv"
FROZEN_SHA = "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"

print("=" * 110)
print("STEP 97: FINAL BIOLOGICAL UNIQUENESS RECONCILIATION")
print("=" * 110)

master_hash = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash == FROZEN_SHA, "FROZEN MASTER SHA256 MISMATCH"
print(f"Frozen master SHA256: PASS")

master = pd.read_csv(MASTER_CSV, dtype=str, low_memory=False)
v2 = pd.read_csv(V2_CSV, dtype=str, low_memory=False)
true = pd.read_csv(STEP83_CSV, dtype=str, low_memory=False)

bio_counts = v2["Final_biological_identity"].value_counts()
bio_nonempty = v2[v2["Final_biological_identity"].apply(lambda x: str(x).strip() not in {"", "nan", "none"})]

unique_bio = bio_nonempty["Final_biological_identity"].nunique()
total_bio = len(bio_nonempty)
duplicate_bio = bio_counts[bio_counts > 1]

print(f"\nTotal KJUs with biological identity: {total_bio}")
print(f"Unique biological identities: {unique_bio}")

bio_audit_rows = []
for bio_id, count in bio_counts.items():
    if str(bio_id).strip() in {"", "nan", "none"}:
        continue
    kjus = bio_nonempty[bio_nonempty["Final_biological_identity"] == bio_id]["KJU_ID"].tolist()
    bio_audit_rows.append({
        "Biological_identity": bio_id,
        "KJU_count": count,
        "KJU_IDs": ";".join(kjus),
        "Is_duplicate": str(count > 1),
    })

bio_audit = pd.DataFrame(bio_audit_rows)
bio_audit.to_csv(OUT / "STEP97_BIOLOGICAL_UNIQUENESS_AUDIT.csv", index=False)

print(f"\nDuplicate biological identities: {len(duplicate_bio)}")
if len(duplicate_bio) > 0:
    for bio_id, count in duplicate_bio.items():
        kjus = bio_nonempty[bio_nonempty["Final_biological_identity"] == bio_id]["KJU_ID"].tolist()
        print(f"  {bio_id}: {count} KJUs ({';'.join(kjus)})")

cross_kju = true[true["KJU_ID"].duplicated(keep=False)]
seq_sha_by_kju = true.groupby("KJU_ID")["Sequence_SHA256"].apply(set).to_dict()

cross_rows = []
checked = set()
for kju_a, shas_a in seq_sha_by_kju.items():
    for kju_b, shas_b in seq_sha_by_kju.items():
        if kju_a >= kju_b:
            continue
        pair = (kju_a, kju_b)
        if pair in checked:
            continue
        common = shas_a & shas_b
        if common:
            checked.add(pair)
            cross_rows.append({
                "KJU_ID_A": kju_a,
                "KJU_ID_B": kju_b,
                "Shared_sequence_SHA256_count": len(common),
                "Shared_SHA256s": ";".join(sorted(common)),
            })

cross_df = pd.DataFrame(cross_rows)
cross_df.to_csv(OUT / "STEP97_CROSS_KJU_RELATIONSHIPS.csv", index=False)

print(f"\nCross-KJU identical sequence pairs: {len(cross_df)}")
if len(cross_df) > 0:
    for _, cr in cross_df.head(10).iterrows():
        print(f"  {cr['KJU_ID_A']} <-> {cr['KJU_ID_B']}: {cr['Shared_sequence_SHA256_count']} shared seqs")

master_hash_after = hashlib.sha256(MASTER_CSV.read_bytes()).hexdigest()
assert master_hash_after == FROZEN_SHA

with open(OUT / "STEP97_REPORT.txt", "w") as f:
    f.write("STEP 97: FINAL BIOLOGICAL UNIQUENESS RECONCILIATION\n")
    f.write("=" * 110 + "\n\n")
    f.write(f"Frozen master SHA256: PASS\n\n")
    f.write(f"KJUs with biological identity: {total_bio}\n")
    f.write(f"Unique biological identities:   {unique_bio}\n")
    f.write(f"Duplicate biological identities:{len(duplicate_bio)}\n")
    f.write(f"Cross-KJU identical seq pairs:  {len(cross_df)}\n\n")
    if len(duplicate_bio) > 0:
        f.write("DUPLICATE BIOLOGICAL IDENTITIES:\n")
        for bio_id, count in duplicate_bio.items():
            kjus = bio_nonempty[bio_nonempty["Final_biological_identity"] == bio_id]["KJU_ID"].tolist()
            f.write(f"  {bio_id}: {count} KJUs ({';'.join(kjus)})\n")
    f.write(f"\nFrozen V1 master untouched.\n")

print("\n" + "=" * 110)
print("STEP 97 COMPLETE")
print("=" * 110)
