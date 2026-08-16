#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import shutil
from datetime import datetime


# ============================================================
# STEP 76
# FINAL KJU MASTER FREEZE
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "final_integrated_master"
)

MASTER = ROOT / "FINAL_204_KJU_ADJUDICATED_MASTER.csv"
BIO = ROOT / "FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER.csv"
QC = ROOT / "STEP75_FINAL_QC_REPORT.txt"

FREEZE = ROOT / "FINAL_KJU_MASTER_FREEZE_V1"
FREEZE.mkdir(parents=True, exist_ok=True)


FINAL_MASTER = FREEZE / "FINAL_204_KJU_ADJUDICATED_MASTER_V1.csv"
FINAL_BIO = FREEZE / "FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER_V1.csv"
FINAL_QC = FREEZE / "STEP75_FINAL_QC_REPORT_V1.txt"
MANIFEST = FREEZE / "FINAL_KJU_MASTER_FREEZE_MANIFEST.txt"


print("=" * 110)
print("STEP 76: FINAL KJU MASTER FREEZE V1")
print("=" * 110)


# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(MASTER, dtype=str).fillna("")
bio = pd.read_csv(BIO, dtype=str).fillna("")


# ============================================================
# BASIC QC
# ============================================================

assert len(df) == 204
assert df["KJU_ID"].nunique() == 204

assert len(bio) == 203
assert bio["Biological_identity"].nunique() == 203


# ============================================================
# TAXONOMY QC
# ============================================================

resolved = df[
    df["Adjudication_type"] != "MASTER_ONLY_UNRESOLVED"
]

unresolved = df[
    df["Adjudication_type"] == "MASTER_ONLY_UNRESOLVED"
]


assert len(resolved) == 112
assert len(unresolved) == 92


# No resolved record may have blank taxonomy
bad_resolved_taxonomy = resolved[
    resolved["Accepted_taxonomy"].str.strip() == ""
]

assert len(bad_resolved_taxonomy) == 0


# Every unresolved record must have blank taxonomy
bad_unresolved_taxonomy = unresolved[
    unresolved["Accepted_taxonomy"].str.strip() != ""
]

assert len(bad_unresolved_taxonomy) == 0


# Every unresolved record must have unique biological identity
bad_unresolved_identity = unresolved[
    ~unresolved["Biological_identity"].str.startswith("UNRESOLVED:")
]

assert len(bad_unresolved_identity) == 0


# ============================================================
# KNOWN DUPLICATE QC
# ============================================================

dup = df[
    df["Biological_identity"] == "PZ666011"
]

assert set(dup["KJU_ID"]) == {"KJU0038", "KJU0039"}


# ============================================================
# ACCESSION QC
# ============================================================

resolved_accessions = (
    resolved["Accepted_accession"]
    .str.strip()
)

assert not resolved_accessions.eq("").any()


# ============================================================
# COPY FINAL FILES
# ============================================================

shutil.copy2(MASTER, FINAL_MASTER)
shutil.copy2(BIO, FINAL_BIO)
shutil.copy2(QC, FINAL_QC)


# ============================================================
# MANIFEST
# ============================================================

timestamp = datetime.now().isoformat(timespec="seconds")

lines = []

lines.append("FINAL KJU MASTER FREEZE V1")
lines.append("=" * 100)
lines.append("")

lines.append(f"Freeze timestamp : {timestamp}")
lines.append("")

lines.append("INPUT MASTER")
lines.append(str(MASTER))
lines.append("")

lines.append("FINAL MASTER")
lines.append(str(FINAL_MASTER))
lines.append("")

lines.append("BIOLOGICALLY UNIQUE MASTER")
lines.append(str(FINAL_BIO))
lines.append("")

lines.append("QC REPORT")
lines.append(str(FINAL_QC))
lines.append("")

lines.append("DATASET COUNTS")
lines.append("-" * 100)
lines.append(f"Total KJU records          : {len(df)}")
lines.append(f"Unique KJU IDs             : {df['KJU_ID'].nunique()}")
lines.append(f"Resolved records           : {len(resolved)}")
lines.append(f"Unresolved MASTER_ONLY     : {len(unresolved)}")
lines.append(f"Biological identities      : {bio['Biological_identity'].nunique()}")
lines.append("")

lines.append("TAXONOMY POLICY")
lines.append("-" * 100)
lines.append("Resolved records have accepted taxonomy.")
lines.append("MASTER_ONLY unresolved records retain blank taxonomy.")
lines.append("No taxonomy was fabricated.")
lines.append("")

lines.append("KNOWN BIOLOGICAL DUPLICATE")
lines.append("-" * 100)
lines.append("KJU0038 + KJU0039 -> PZ666011")
lines.append("")

lines.append("QC STATUS")
lines.append("-" * 100)
lines.append("STEP 75 QC PASSED")
lines.append("FINAL MASTER FREEZE V1 CREATED")
lines.append("")


MANIFEST.write_text("\n".join(lines) + "\n")


# ============================================================
# FINAL REPORT
# ============================================================

print()
print("=" * 110)
print("STEP 76 COMPLETE")
print("=" * 110)

print()
print("FINAL MASTER:")
print(FINAL_MASTER)

print()
print("BIOLOGICALLY UNIQUE MASTER:")
print(FINAL_BIO)

print()
print("QC REPORT:")
print(FINAL_QC)

print()
print("MANIFEST:")
print(MANIFEST)

print()
print("COUNTS")
print("-" * 100)
print(f"Total KJU records      : {len(df)}")
print(f"Resolved records       : {len(resolved)}")
print(f"Unresolved records     : {len(unresolved)}")
print(f"Biological identities  : {bio['Biological_identity'].nunique()}")

print()
print("STATUS: FINAL MASTER FROZEN — V1")
print("=" * 110)
