from pathlib import Path
import hashlib
import shutil
import pandas as pd
from datetime import datetime

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

FREEZE = META / "NCBI_RESOLUTION_FREEZE_V1"
FREEZE.mkdir(exist_ok=True)

FILES = [
    "MASTER_KJU_RECONCILIATION_001_183.csv",
    "NCBI_VALIDATION_DISPUTED_RECORDS.csv",
    "NCBI_RESOLUTION_TABLE_V1.csv",
    "NCBI_RESOLUTION_SUMMARY_V1.txt",
    "KJU0039_INVESTIGATION_V1.csv",
    "KJU0039_INVESTIGATION_SUMMARY_V1.txt",
    "KJU0039_FINAL_RESOLUTION_V1.csv",
    "KJU0039_FINAL_RESOLUTION_SUMMARY_V1.txt",
    "NCBI_MASTER_CHANGE_PROPOSAL_V1.csv",
    "NCBI_MASTER_CHANGE_PROPOSAL_SUMMARY_V1.txt",
]

print("=" * 90)
print(" DODDAGGUBBI - FREEZE NCBI RESOLUTION STATE")
print("=" * 90)

manifest = []

for filename in FILES:

    src = META / filename

    if not src.exists():
        print(f"WARNING: missing {filename}")
        continue

    dst = FREEZE / filename

    shutil.copy2(src, dst)

    sha256 = hashlib.sha256(
        dst.read_bytes()
    ).hexdigest()

    size = dst.stat().st_size

    manifest.append({
        "file": filename,
        "size_bytes": size,
        "sha256": sha256,
    })

    print(f"FROZEN: {filename}")

manifest_df = pd.DataFrame(manifest)

manifest_file = FREEZE / "SHA256_MANIFEST_V1.csv"
manifest_df.to_csv(manifest_file, index=False)

# ------------------------------------------------------------
# Additional master validation
# ------------------------------------------------------------

master_file = FREEZE / "MASTER_KJU_RECONCILIATION_001_183.csv"

df = pd.read_csv(
    master_file,
    dtype=str
).fillna("")

checks = {
    "master_rows": len(df),
    "master_unique_kju": df["kju_id"].nunique(),
    "duplicate_kju_ids": int(df["kju_id"].duplicated().sum()),
    "expected_kju_ids": 183,
    "kju0038_accession": (
        df.loc[df["kju_id"] == "KJU0038", "accession"].iloc[0]
        if (df["kju_id"] == "KJU0038").any()
        else ""
    ),
    "kju0038_organism": (
        df.loc[df["kju_id"] == "KJU0038", "organism"].iloc[0]
        if (df["kju_id"] == "KJU0038").any()
        else ""
    ),
    "kju0039_accession": (
        df.loc[df["kju_id"] == "KJU0039", "accession"].iloc[0]
        if (df["kju_id"] == "KJU0039").any()
        else ""
    ),
    "kju0039_status": (
        df.loc[df["kju_id"] == "KJU0039", "status"].iloc[0]
        if (df["kju_id"] == "KJU0039").any()
        else ""
    ),
}

check_file = FREEZE / "MASTER_INTEGRITY_CHECK_V1.txt"

with check_file.open("w", encoding="utf-8") as fh:

    fh.write("DODDAGGUBBI - FROZEN MASTER INTEGRITY CHECK\n")
    fh.write("=" * 80 + "\n\n")

    for key, value in checks.items():
        fh.write(f"{key}: {value}\n")

    fh.write("\n")

    if (
        checks["master_rows"] == 183
        and checks["master_unique_kju"] == 183
        and checks["duplicate_kju_ids"] == 0
        and checks["kju0038_accession"] == "PZ666011.1"
        and checks["kju0038_organism"] == "Shewanella algae"
        and checks["kju0039_accession"] == ""
        and checks["kju0039_status"] == "NOT_CURRENTLY_DOWNLOADED"
    ):
        fh.write("STATUS: PASS\n")
    else:
        fh.write("STATUS: FAIL\n")

print("")
print("=" * 90)
print(" FREEZE COMPLETE")
print("=" * 90)

print(f"Frozen directory:")
print(FREEZE)

print("")
print(f"SHA256 manifest:")
print(manifest_file)

print("")
print(f"Integrity report:")
print(check_file)

print("")
print("MASTER DATA WAS NOT MODIFIED.")
