from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

# ------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------

OPENCITY = Path("/mnt/c/Users/saipr/Desktop/731b3fb2-c763-4491-9621-e04f9e1a9740.csv")

RECON_FILE = META / "MASTER_KJU_RECONCILIATION_001_183.csv"

OUT = META

print("=" * 80)
print(" DODDAGGUBBI - OPENCITY ↔ MASTER KJU RECONCILIATION")
print("=" * 80)

# ------------------------------------------------------------
# 1. READ OPENCITY
# ------------------------------------------------------------

print("\n[1/7] Reading OpenCity dataset")
print("-" * 80)

if not OPENCITY.exists():
    raise SystemExit(
        f"ERROR: OpenCity file not found:\n{OPENCITY}"
    )

if not RECON_FILE.exists():
    raise SystemExit(
        f"ERROR: Master reconciliation file not found:\n{RECON_FILE}"
    )

oc = pd.read_csv(OPENCITY)

print(f"OpenCity records : {len(oc)}")
print(f"Columns          : {list(oc.columns)}")

# ------------------------------------------------------------
# 2. STANDARDIZE OPENCITY COLUMN NAMES
# ------------------------------------------------------------

print("\n[2/7] Standardizing OpenCity columns")
print("-" * 80)

rename_map = {
    "Strain ID": "kju_id",
    "GenBank Acc. No.": "accession",
    "Species": "organism",
    "Hyperlink": "ncbi_url",
}

oc = oc.rename(columns=rename_map)

print("Standardized columns:")
print(list(oc.columns))

# ------------------------------------------------------------
# 3. CLEAN IDENTIFIERS
# ------------------------------------------------------------

print("\n[3/7] Cleaning identifiers")
print("-" * 80)

if "kju_id" not in oc.columns:
    raise SystemExit(
        "ERROR: OpenCity dataset does not contain 'Strain ID' column."
    )

if "accession" not in oc.columns:
    raise SystemExit(
        "ERROR: OpenCity dataset does not contain 'GenBank Acc. No.' column."
    )

oc["kju_id"] = (
    oc["kju_id"]
    .astype(str)
    .str.strip()
    .str.upper()
)

oc["accession"] = (
    oc["accession"]
    .astype(str)
    .str.strip()
)

# Normalize accession where possible
oc["accession_base"] = (
    oc["accession"]
    .str.replace(r"\s+", "", regex=True)
    .str.replace(r"\.$", "", regex=True)
)

# Extract valid KJU IDs
oc["kju_valid"] = oc["kju_id"].str.fullmatch(
    r"KJU\d{4}",
    na=False
)

print(f"Valid KJU IDs: {oc['kju_valid'].sum()}")

# ------------------------------------------------------------
# 4. DUPLICATE CHECKS
# ------------------------------------------------------------

print("\n[4/7] Checking OpenCity duplicates")
print("-" * 80)

dup_kju = oc[oc["kju_id"].duplicated(keep=False)]
dup_acc = oc[oc["accession_base"].duplicated(keep=False)]

print(f"Duplicate KJU rows       : {len(dup_kju)}")
print(f"Duplicate accession rows : {len(dup_acc)}")

if len(dup_kju) > 0:
    print("\nDuplicate KJU records:")
    print(dup_kju.to_string(index=False))

if len(dup_acc) > 0:
    print("\nDuplicate accession records:")
    print(dup_acc.to_string(index=False))

# ------------------------------------------------------------
# 5. READ MASTER RECONCILIATION
# ------------------------------------------------------------

print("\n[5/7] Comparing with MASTER 001-183 reconciliation")
print("-" * 80)

master = pd.read_csv(RECON_FILE)

print(f"Master reconciliation records: {len(master)}")

master["kju_id"] = (
    master["kju_id"]
    .astype(str)
    .str.strip()
    .str.upper()
)

if "accession" in master.columns:
    master["accession"] = (
        master["accession"]
        .astype(str)
        .str.strip()
    )

# ------------------------------------------------------------
# 6. KJU-BY-KJU RECONCILIATION
# ------------------------------------------------------------

print("\n[6/7] Building OpenCity ↔ MASTER reconciliation")
print("-" * 80)

# Select useful OpenCity columns
oc_cols = [
    c for c in [
        "kju_id",
        "accession",
        "accession_base",
        "organism",
        "ncbi_url"
    ]
    if c in oc.columns
]

oc_small = oc[oc_cols].copy()

# Prevent accidental many-to-many merge
oc_small = oc_small.drop_duplicates(
    subset=["kju_id"],
    keep="first"
)

# Prefix OpenCity columns
oc_small = oc_small.rename(
    columns={
        c: f"opencity_{c}"
        for c in oc_small.columns
        if c != "kju_id"
    }
)

merged = master.merge(
    oc_small,
    left_on="kju_id",
    right_on="kju_id",
    how="outer",
    indicator=True
)

# ------------------------------------------------------------
# Determine reconciliation status
# ------------------------------------------------------------

def classify(row):

    master_present = row["_merge"] in ["both", "left_only"]
    oc_present = row["_merge"] in ["both", "right_only"]

    if row["_merge"] == "both":
        master_acc = str(row.get("accession", "")).strip()
        oc_acc = str(row.get("opencity_accession", "")).strip()

        master_base = re.sub(r"\.$", "", master_acc)
        oc_base = re.sub(r"\.$", "", oc_acc)

        if master_base and oc_base and master_base == oc_base:
            return "MATCH_KJU_ACCESSION"

        return "MATCH_KJU_ACCESSION_MISMATCH"

    if row["_merge"] == "left_only":
        return "MASTER_ONLY"

    if row["_merge"] == "right_only":
        return "OPENCITY_ONLY"

    return "UNKNOWN"


merged["reconciliation_status"] = merged.apply(
    classify,
    axis=1
)

# ------------------------------------------------------------
# Additional comparison fields
# ------------------------------------------------------------

if "organism" in merged.columns and "opencity_organism" in merged.columns:

    merged["organism_match"] = (
        merged["organism"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
        ==
        merged["opencity_organism"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )

else:
    merged["organism_match"] = pd.NA

# ------------------------------------------------------------
# Save complete reconciliation
# ------------------------------------------------------------

complete_file = OUT / "OPENCITY_MASTER_KJU_RECONCILIATION.csv"

merged.to_csv(
    complete_file,
    index=False
)

# ------------------------------------------------------------
# Separate outputs
# ------------------------------------------------------------

matches = merged[
    merged["reconciliation_status"] == "MATCH_KJU_ACCESSION"
].copy()

acc_mismatch = merged[
    merged["reconciliation_status"] == "MATCH_KJU_ACCESSION_MISMATCH"
].copy()

master_only = merged[
    merged["reconciliation_status"] == "MASTER_ONLY"
].copy()

opencity_only = merged[
    merged["reconciliation_status"] == "OPENCITY_ONLY"
].copy()

matches.to_csv(
    OUT / "OPENCITY_MASTER_MATCHES.csv",
    index=False
)

acc_mismatch.to_csv(
    OUT / "OPENCITY_ACCESSION_MISMATCH.csv",
    index=False
)

master_only.to_csv(
    OUT / "OPENCITY_MASTER_ONLY.csv",
    index=False
)

opencity_only.to_csv(
    OUT / "OPENCITY_ONLY.csv",
    index=False
)

# ------------------------------------------------------------
# 7. FINAL SUMMARY
# ------------------------------------------------------------

print("\n[7/7] FINAL RECONCILIATION SUMMARY")
print("=" * 80)

print(f"OpenCity records                  : {len(oc)}")
print(f"Master reconciliation records     : {len(master)}")
print(f"Combined reconciliation records   : {len(merged)}")

print("\nStatus breakdown:")
print(
    merged["reconciliation_status"]
    .value_counts()
    .to_string()
)

print("\nOrganism comparison:")
print(
    merged["organism_match"]
    .value_counts(dropna=False)
    .to_string()
)

print("\n" + "=" * 80)
print(" IMPORTANT RESULTS")
print("=" * 80)

print(
    f"\nExact KJU + accession matches : {len(matches)}"
)

print(
    f"KJU matched but accession differs : {len(acc_mismatch)}"
)

print(
    f"Master-only KJU IDs           : {len(master_only)}"
)

print(
    f"OpenCity-only KJU IDs         : {len(opencity_only)}"
)

if len(acc_mismatch) > 0:
    print("\nACCESSION MISMATCHES:")
    print("-" * 80)

    cols = [
        c for c in [
            "kju_id",
            "accession",
            "opencity_accession",
            "organism",
            "opencity_organism"
        ]
        if c in acc_mismatch.columns
    ]

    print(
        acc_mismatch[cols]
        .to_string(index=False)
    )

if len(opencity_only) > 0:
    print("\nKJU IDs PRESENT IN OPENCITY BUT NOT MASTER:")
    print("-" * 80)

    print(
        opencity_only[
            [
                c for c in [
                    "kju_id",
                    "opencity_accession",
                    "opencity_organism"
                ]
                if c in opencity_only.columns
            ]
        ].to_string(index=False)
    )

print("\n" + "=" * 80)
print(" RECONCILIATION COMPLETE")
print("=" * 80)

print("\nOutput files:")
print("-" * 80)

print(complete_file)
print(OUT / "OPENCITY_MASTER_MATCHES.csv")
print(OUT / "OPENCITY_ACCESSION_MISMATCH.csv")
print(OUT / "OPENCITY_MASTER_ONLY.csv")
print(OUT / "OPENCITY_ONLY.csv")

print("\n" + "=" * 80)
