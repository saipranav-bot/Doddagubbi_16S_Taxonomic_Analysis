from pathlib import Path
import hashlib
import shutil
import pandas as pd
from datetime import datetime

# ============================================================
# DODDAGGUBBI
# FREEZE CANONICAL KJU METADATA
#
# SAFE / NON-DESTRUCTIVE
# ============================================================

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

CANONICAL = META / "CANONICAL_KJU_METADATA_V1.csv"
SUMMARY = META / "CANONICAL_KJU_METADATA_SUMMARY_V1.txt"

FREEZE = META / "CANONICAL_KJU_METADATA_FREEZE_V1"

FREEZE.mkdir(exist_ok=True)

print("=" * 90)
print(" DODDAGGUBBI - FREEZE CANONICAL KJU METADATA")
print("=" * 90)

# ------------------------------------------------------------
# Validate canonical file
# ------------------------------------------------------------

if not CANONICAL.exists():
    raise SystemExit(
        f"ERROR: Missing canonical file:\n{CANONICAL}"
    )

df = pd.read_csv(
    CANONICAL,
    dtype=str
).fillna("")

# ------------------------------------------------------------
# Structural validation
# ------------------------------------------------------------

if len(df) != 183:
    raise SystemExit(
        f"ERROR: Expected 183 records, found {len(df)}"
    )

if df["kju_id"].duplicated().any():
    raise SystemExit(
        "ERROR: Duplicate KJU IDs detected."
    )

if df["kju_id"].nunique() != 183:
    raise SystemExit(
        "ERROR: Unique KJU ID count is not 183."
    )

# ------------------------------------------------------------
# KJU0038 safety
# ------------------------------------------------------------

r38 = df[df["kju_id"] == "KJU0038"]

if len(r38) != 1:
    raise SystemExit(
        "ERROR: KJU0038 missing."
    )

r38 = r38.iloc[0]

if r38["canonical_accession"] != "PZ666011.1":
    raise SystemExit(
        "ERROR: KJU0038 accession safety check failed."
    )

if r38["canonical_organism"] != "Shewanella algae":
    raise SystemExit(
        "ERROR: KJU0038 organism safety check failed."
    )

# ------------------------------------------------------------
# KJU0039 safety
# ------------------------------------------------------------

r39 = df[df["kju_id"] == "KJU0039"]

if len(r39) != 1:
    raise SystemExit(
        "ERROR: KJU0039 missing."
    )

r39 = r39.iloc[0]

if r39["canonical_accession"] != "":
    raise SystemExit(
        "SAFETY ERROR: KJU0039 has canonical accession."
    )

if r39["canonical_organism"] != "":
    raise SystemExit(
        "SAFETY ERROR: KJU0039 has canonical organism."
    )

if r39["resolution_status"] != \
        "UNRESOLVED_ACCESSION_COLLISION":

    raise SystemExit(
        "SAFETY ERROR: KJU0039 resolution changed."
    )

# ------------------------------------------------------------
# Freeze canonical CSV
# ------------------------------------------------------------

frozen_csv = FREEZE / CANONICAL.name

shutil.copy2(
    CANONICAL,
    frozen_csv
)

print(
    f"FROZEN: {CANONICAL.name}"
)

# ------------------------------------------------------------
# Freeze summary
# ------------------------------------------------------------

if SUMMARY.exists():

    frozen_summary = FREEZE / SUMMARY.name

    shutil.copy2(
        SUMMARY,
        frozen_summary
    )

    print(
        f"FROZEN: {SUMMARY.name}"
    )

# ------------------------------------------------------------
# SHA256 manifest
# ------------------------------------------------------------

manifest_rows = []

for file in sorted(FREEZE.iterdir()):

    if not file.is_file():
        continue

    if file.name == "SHA256_MANIFEST_V1.csv":
        continue

    sha256 = hashlib.sha256()

    with open(file, "rb") as fh:

        for chunk in iter(
            lambda: fh.read(1024 * 1024),
            b""
        ):
            sha256.update(chunk)

    manifest_rows.append({
        "filename": file.name,
        "bytes": file.stat().st_size,
        "sha256": sha256.hexdigest(),
    })

manifest = pd.DataFrame(manifest_rows)

manifest_file = (
    FREEZE / "SHA256_MANIFEST_V1.csv"
)

manifest.to_csv(
    manifest_file,
    index=False
)

print("")
print("SHA256 manifest created:")
print(manifest_file)

# ------------------------------------------------------------
# Integrity report
# ------------------------------------------------------------

integrity = []

integrity.append(
    "DODDAGGUBBI - CANONICAL METADATA FREEZE"
)

integrity.append("=" * 80)

integrity.append("")

integrity.append(
    f"Freeze timestamp: "
    f"{datetime.now().isoformat(timespec='seconds')}"
)

integrity.append(
    "Canonical records: 183"
)

integrity.append(
    "Unique KJU IDs: 183"
)

integrity.append(
    "Duplicate KJU IDs: 0"
)

integrity.append("")

integrity.append("KJU0038")

integrity.append("-" * 80)

integrity.append(
    "canonical_accession = PZ666011.1"
)

integrity.append(
    "canonical_organism = Shewanella algae"
)

integrity.append(
    "status = NCBI_MASTER_SUPPORTED"
)

integrity.append("")

integrity.append("KJU0039")

integrity.append("-" * 80)

integrity.append(
    "canonical_accession = [BLANK]"
)

integrity.append(
    "canonical_organism = [BLANK]"
)

integrity.append(
    "status = UNRESOLVED_ACCESSION_COLLISION"
)

integrity.append(
    "action = DO_NOT_MERGE"
)

integrity.append("")

integrity.append("IMPORTANT")

integrity.append("-" * 80)

integrity.append(
    "MASTER_KJU_RECONCILIATION_001_183.csv "
    "was NOT modified."
)

integrity.append(
    "Canonical metadata is a derived downstream layer."
)

integrity.append(
    "KJU0039 was NOT assigned PZ666011.1."
)

integrity.append(
    "This directory represents the frozen canonical state."
)

integrity_file = (
    FREEZE / "CANONICAL_INTEGRITY_CHECK_V1.txt"
)

integrity_file.write_text(
    "\n".join(integrity) + "\n",
    encoding="utf-8"
)

print(
    f"Integrity report: {integrity_file}"
)

# ------------------------------------------------------------
# Final output
# ------------------------------------------------------------

print("")
print("=" * 90)
print(" CANONICAL FREEZE COMPLETE")
print("=" * 90)

print("")
print("Frozen directory:")
print(FREEZE)

print("")
print("Files frozen:")

for row in manifest_rows:
    print(
        f"  {row['filename']} "
        f"({row['bytes']} bytes)"
    )

print("")
print("KJU0038: PASS")
print("KJU0039: PASS")
print("183 KJU records: PASS")

print("")
print("MASTER DATA WAS NOT MODIFIED.")
