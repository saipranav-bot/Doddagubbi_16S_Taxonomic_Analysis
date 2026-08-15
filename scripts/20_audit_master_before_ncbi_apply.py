from pathlib import Path
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

MASTER_FILE = META / "MASTER_KJU_RECONCILIATION_001_183.csv"
NCBI_RES = META / "NCBI_RESOLUTION_TABLE_V1.csv"
KJU0039_RES = META / "KJU0039_FINAL_RESOLUTION_V1.csv"

print("=" * 90)
print(" DODDAGGUBBI - PRE-APPLICATION MASTER AUDIT")
print("=" * 90)

for f in [MASTER_FILE, NCBI_RES, KJU0039_RES]:
    print(f"\nFILE: {f}")
    print("EXISTS:", f.exists())

    if not f.exists():
        continue

    df = pd.read_csv(f, dtype=str).fillna("")

    print("Rows:", len(df))
    print("Columns:", len(df.columns))
    print("Column names:")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i:02d}. {col}")

# ------------------------------------------------------------
# Master KJU structure
# ------------------------------------------------------------

master = pd.read_csv(MASTER_FILE, dtype=str).fillna("")
resolution = pd.read_csv(NCBI_RES, dtype=str).fillna("")

print("\n" + "=" * 90)
print(" MASTER DATASET CHECK")
print("=" * 90)

if "kju_id" in master.columns:
    print("KJU column: kju_id")
elif "KJU_ID" in master.columns:
    print("KJU column: KJU_ID")
else:
    print("WARNING: Could not identify KJU column")

# ------------------------------------------------------------
# Duplicate KJU IDs
# ------------------------------------------------------------

kju_col = None

for candidate in ["kju_id", "KJU_ID", "KJU"]:
    if candidate in master.columns:
        kju_col = candidate
        break

if kju_col:
    duplicates = master[
        master[kju_col].duplicated(keep=False)
    ]

    print("\nKJU duplicates:", len(duplicates))

    if not duplicates.empty:
        print(duplicates[[kju_col]].value_counts())

# ------------------------------------------------------------
# Expected range
# ------------------------------------------------------------

if kju_col:
    ids = master[kju_col].tolist()

    expected = {
        f"KJU{i:04d}"
        for i in range(1, 184)
    }

    observed = set(x for x in ids if x)

    missing = sorted(expected - observed)
    extra = sorted(observed - expected)

    print("\nExpected KJU IDs: 183")
    print("Observed KJU IDs:", len(observed))
    print("Missing:", len(missing))
    print("Extra:", len(extra))

    if missing:
        print("\nMissing KJU IDs:")
        print(", ".join(missing))

    if extra:
        print("\nExtra KJU IDs:")
        print(", ".join(extra))

# ------------------------------------------------------------
# NCBI resolution coverage
# ------------------------------------------------------------

print("\n" + "=" * 90)
print(" NCBI RESOLUTION COVERAGE")
print("=" * 90)

if "kju_id" in resolution:

    print("Resolution records:", len(resolution))

    print("\nResolution classifications:")

    if "classification" in resolution.columns:
        print(
            resolution["classification"]
            .value_counts()
            .to_string()
        )

    if "resolution_action" in resolution.columns:
        print(
            resolution["resolution_action"]
            .value_counts()
            .to_string()
        )

# ------------------------------------------------------------
# Critical KJU0039 check
# ------------------------------------------------------------

print("\n" + "=" * 90)
print(" KJU0039 SAFETY CHECK")
print("=" * 90)

if kju_col:

    k39 = master[
        master[kju_col] == "KJU0039"
    ]

    if k39.empty:
        print("KJU0039 is NOT present in master.")

    else:
        print(k39.to_string(index=False))

print("\n" + "=" * 90)
print(" AUDIT COMPLETE")
print("=" * 90)

print("\nNO FILES WERE MODIFIED.")

