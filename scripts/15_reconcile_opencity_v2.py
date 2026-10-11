from pathlib import Path
import pandas as pd
import re

# ============================================================
# DODDAGGUBBI
# OpenCity ↔ MASTER KJU Reconciliation v2
#
# Purpose:
#   1. Compare OpenCity against MASTER 001-183
#   2. Normalize accession version suffixes (.1, .2, ...)
#   3. Preserve original accession values
#   4. Separate accession-version differences from real conflicts
#   5. Detect organism mismatches
#   6. Detect obvious organism spelling/case differences
#   7. Identify MASTER_ONLY and OPENCITY_ONLY records
#   8. Produce audit-ready output files
# ============================================================

from pathlib import Path
import pandas as pd
import re
import sys
import os

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

MASTER = Path(os.environ.get("DODDAGUBBI_MASTER", Path.home() / "Doddagubbi_Master"))
META = MASTER / "07_Master_Metadata"

# External OpenCity export. Override with OPENCITY_CSV if stored elsewhere.
OPENCITY = Path(
    os.environ.get("OPENCITY_CSV", MASTER / "opencity_731b3fb2.csv")
)

RECON_FILE = META / "MASTER_KJU_RECONCILIATION_001_183.csv"

OUT = META

# ------------------------------------------------------------
# OUTPUT FILES
# ------------------------------------------------------------

FULL_OUT = OUT / "OPENCITY_MASTER_RECONCILIATION_V2.csv"
MATCH_OUT = OUT / "OPENCITY_MASTER_MATCHES_V2.csv"
VERSION_OUT = OUT / "OPENCITY_ACCESSION_VERSION_ONLY_V2.csv"
ORGANISM_OUT = OUT / "OPENCITY_ORGANISM_CONFLICTS_V2.csv"
CONFLICT_OUT = OUT / "OPENCITY_TRUE_CONFLICTS_V2.csv"
MASTER_ONLY_OUT = OUT / "OPENCITY_MASTER_ONLY_V2.csv"
OPENCITY_ONLY_OUT = OUT / "OPENCITY_ONLY_V2.csv"
SUMMARY_OUT = OUT / "OPENCITY_RECONCILIATION_SUMMARY_V2.txt"

# ============================================================
# HELPERS
# ============================================================

def clean_text(value):
    """
    Convert values to clean strings.
    Preserve missing values as None.
    """
    if pd.isna(value):
        return None

    value = str(value).strip()

    if value == "":
        return None

    if value.lower() in {
        "nan",
        "none",
        "null",
        "na",
        "n/a",
    }:
        return None

    return value


def normalize_kju(value):
    """
    Normalize KJU identifiers.

    Examples:
        KJU0001 -> KJU0001
        kju0001 -> KJU0001
        KJU 0001 -> KJU0001
    """
    value = clean_text(value)

    if value is None:
        return None

    value = value.upper()
    value = re.sub(r"\s+", "", value)

    match = re.search(r"KJU0*(\d+)", value)

    if match:
        return f"KJU{int(match.group(1)):04d}"

    return value


def normalize_accession(value):
    """
    Normalize accession ONLY for comparison.

    Examples:
        PZ638859.1 -> PZ638859
        PZ638859.2 -> PZ638859
        PZ638859   -> PZ638859

    Original accession is preserved separately.
    """
    value = clean_text(value)

    if value is None:
        return None

    value = value.strip().upper()

    # Remove whitespace
    value = re.sub(r"\s+", "", value)

    # Remove NCBI accession version suffix
    value = re.sub(r"\.\d+$", "", value)

    return value


def normalize_organism(value):
    """
    Conservative organism normalization.

    This DOES NOT automatically correct taxonomy.

    It only normalizes:
        - case
        - whitespace
        - punctuation spacing

    Scientific names are otherwise preserved.
    """
    value = clean_text(value)

    if value is None:
        return None

    value = value.lower().strip()

    value = re.sub(r"\s+", " ", value)

    value = value.replace(",", "")

    return value


def organism_first_two_words(value):
    """
    Return first two tokens of organism name.

    Example:
        'Aeromonas caviae strain X'
        -> 'aeromonas caviae'

    Used only as supporting information.
    """
    value = normalize_organism(value)

    if value is None:
        return None

    tokens = value.split()

    if len(tokens) >= 2:
        return " ".join(tokens[:2])

    return value


def safe_bool(value):
    """
    Convert nullable comparison to explicit boolean.
    """
    return bool(value) if pd.notna(value) else False


# ============================================================
# HEADER
# ============================================================

print("=" * 90)
print(" DODDAGGUBBI - OPENCITY ↔ MASTER KJU RECONCILIATION v2")
print("=" * 90)

print(f"\nOpenCity file : {OPENCITY}")
print(f"Master file   : {RECON_FILE}")
print(f"Output folder : {OUT}")

# ============================================================
# 1. CHECK INPUTS
# ============================================================

print("\n[1/10] Checking input files")
print("-" * 90)

if not OPENCITY.exists():
    raise SystemExit(
        f"\nERROR: OpenCity file not found:\n{OPENCITY}\n"
    )

if not RECON_FILE.exists():
    raise SystemExit(
        f"\nERROR: Master reconciliation file not found:\n{RECON_FILE}\n"
    )

print("OpenCity file : FOUND")
print("Master file   : FOUND")

# ============================================================
# 2. READ OPENCITY
# ============================================================

print("\n[2/10] Reading OpenCity dataset")
print("-" * 90)

try:
    oc = pd.read_csv(OPENCITY, dtype=str)
except Exception as e:
    raise SystemExit(
        f"\nERROR reading OpenCity CSV:\n{e}"
    )

print(f"OpenCity records : {len(oc)}")
print(f"Columns          : {list(oc.columns)}")

# ============================================================
# 3. STANDARDIZE OPENCITY COLUMNS
# ============================================================

print("\n[3/10] Standardizing OpenCity columns")
print("-" * 90)

oc_column_map = {}

for col in oc.columns:
    c = col.strip().lower()

    if c in {"strain id", "strain_id", "kju id", "kju_id", "kju"}:
        oc_column_map[col] = "kju_id"

    elif c in {
        "genbank acc. no.",
        "genbank acc no.",
        "genbank accession",
        "accession",
        "genbank acc",
        "genbank accession no.",
    }:
        oc_column_map[col] = "accession"

    elif c in {
        "species",
        "organism",
        "species name",
        "organism name",
    }:
        oc_column_map[col] = "organism"

    elif c in {
        "hyperlink",
        "ncbi url",
        "ncbi link",
        "url",
    }:
        oc_column_map[col] = "ncbi_url"

    elif c == "_id":
        oc_column_map[col] = "_id"

oc = oc.rename(columns=oc_column_map)

required_oc = {
    "kju_id",
    "accession",
    "organism",
}

missing_oc = required_oc - set(oc.columns)

if missing_oc:
    raise SystemExit(
        "\nERROR: Missing required OpenCity columns:\n"
        + "\n".join(sorted(missing_oc))
    )

for col in ["_id", "kju_id", "accession", "organism", "ncbi_url"]:
    if col not in oc.columns:
        oc[col] = None

oc = oc[
    [
        "_id",
        "kju_id",
        "accession",
        "organism",
        "ncbi_url",
    ]
].copy()

print("Standardized columns:")
print(list(oc.columns))

# ============================================================
# 4. CLEAN OPENCITY
# ============================================================

print("\n[4/10] Cleaning and normalizing OpenCity identifiers")
print("-" * 90)

oc["kju_id_original"] = oc["kju_id"]
oc["accession_original"] = oc["accession"]
oc["organism_original"] = oc["organism"]

oc["kju_id"] = oc["kju_id"].apply(normalize_kju)

oc["accession_normalized"] = oc["accession"].apply(
    normalize_accession
)

oc["organism_normalized"] = oc["organism"].apply(
    normalize_organism
)

oc["organism_genus_species"] = oc["organism"].apply(
    organism_first_two_words
)

valid_kju = oc["kju_id"].notna().sum()

print(f"Valid KJU IDs : {valid_kju}")

invalid_kju = oc[oc["kju_id"].isna()].copy()

if len(invalid_kju) > 0:
    print(f"WARNING: Invalid/missing KJU rows : {len(invalid_kju)}")

# ------------------------------------------------------------
# OpenCity duplicate checks
# ------------------------------------------------------------

print("\nChecking OpenCity duplicates")

kju_dup = oc[
    oc["kju_id"].duplicated(
        keep=False
    )
].sort_values("kju_id")

acc_dup = oc[
    oc["accession_normalized"].notna()
    &
    oc["accession_normalized"].duplicated(
        keep=False
    )
].sort_values("accession_normalized")

print(
    f"Duplicate KJU rows              : "
    f"{len(kju_dup)}"
)

print(
    f"Duplicate normalized accession rows : "
    f"{len(acc_dup)}"
)

if len(kju_dup) > 0:
    print("\nWARNING: Duplicate KJU IDs detected:")
    print(
        kju_dup[
            [
                "kju_id",
                "accession",
                "organism",
            ]
        ].to_string(index=False)
    )

# ============================================================
# 5. READ MASTER
# ============================================================

print("\n[5/10] Reading MASTER 001-183 reconciliation")
print("-" * 90)

try:
    master = pd.read_csv(RECON_FILE, dtype=str)
except Exception as e:
    raise SystemExit(
        f"\nERROR reading master reconciliation CSV:\n{e}"
    )

print(f"Master records : {len(master)}")
print(f"Columns        : {list(master.columns)}")

# ============================================================
# 6. IDENTIFY MASTER COLUMNS
# ============================================================

print("\n[6/10] Standardizing MASTER columns")
print("-" * 90)

master_column_map = {}

for col in master.columns:
    c = col.strip().lower()

    if c in {
        "kju id",
        "kju_id",
        "kju",
        "strain id",
        "strain_id",
    }:
        master_column_map[col] = "kju_id"

    elif c in {
        "accession",
        "genbank accession",
        "genbank acc. no.",
        "genbank acc no.",
        "genbank accession no.",
    }:
        master_column_map[col] = "accession"

    elif c in {
        "organism",
        "species",
        "species name",
    }:
        master_column_map[col] = "organism"

master = master.rename(columns=master_column_map)

if "kju_id" not in master.columns:
    raise SystemExit(
        "\nERROR: Could not identify KJU ID column in MASTER.\n"
        f"Available columns:\n{list(master.columns)}"
    )

# If accession is missing, create it.
if "accession" not in master.columns:
    master["accession"] = None

if "organism" not in master.columns:
    master["organism"] = None

master["kju_id_original"] = master["kju_id"]
master["accession_original"] = master["accession"]
master["organism_original"] = master["organism"]

master["kju_id"] = master["kju_id"].apply(
    normalize_kju
)

master["accession_normalized"] = master["accession"].apply(
    normalize_accession
)

master["organism_normalized"] = master["organism"].apply(
    normalize_organism
)

master["organism_genus_species"] = master["organism"].apply(
    organism_first_two_words
)

print("Master columns standardized.")

# ============================================================
# 7. MASTER DUPLICATE CHECK
# ============================================================

print("\n[7/10] Checking MASTER duplicates")
print("-" * 90)

master_kju_dup = master[
    master["kju_id"].duplicated(
        keep=False
    )
].sort_values("kju_id")

print(
    f"Duplicate MASTER KJU rows : "
    f"{len(master_kju_dup)}"
)

if len(master_kju_dup) > 0:
    print("\nWARNING: Duplicate MASTER KJU IDs detected:")
    print(
        master_kju_dup[
            [
                "kju_id",
                "accession",
                "organism",
            ]
        ].to_string(index=False)
    )

# ============================================================
# 8. BUILD RECONCILIATION
# ============================================================

print("\n[8/10] Building detailed reconciliation")
print("-" * 90)

# ------------------------------------------------------------
# Keep only one master row per KJU for the primary comparison.
#
# If duplicates exist, the first row is retained for the main
# comparison while duplicate information remains in audit data.
# ------------------------------------------------------------

master_primary = master.drop_duplicates(
    subset=["kju_id"],
    keep="first"
).copy()

oc_primary = oc.drop_duplicates(
    subset=["kju_id"],
    keep="first"
).copy()

# ------------------------------------------------------------
# Outer join on normalized KJU ID
# ------------------------------------------------------------

merged = pd.merge(
    master_primary,
    oc_primary,
    on="kju_id",
    how="outer",
    suffixes=("_master", "_opencity"),
    indicator=True,
)

# ------------------------------------------------------------
# Presence
# ------------------------------------------------------------

merged["in_master"] = (
    merged["_merge"]
    .isin(["both", "left_only"])
)

merged["in_opencity"] = (
    merged["_merge"]
    .isin(["both", "right_only"])
)

# ------------------------------------------------------------
# Accession comparison
# ------------------------------------------------------------

merged["accession_exact_match"] = (
    merged["accession_original_master"].notna()
    &
    merged["accession_original_opencity"].notna()
    &
    (
        merged["accession_original_master"]
        ==
        merged["accession_original_opencity"]
    )
)

merged["accession_normalized_match"] = (
    merged["accession_normalized_master"].notna()
    &
    merged["accession_normalized_opencity"].notna()
    &
    (
        merged["accession_normalized_master"]
        ==
        merged["accession_normalized_opencity"]
    )
)

merged["accession_version_only_difference"] = (
    merged["accession_normalized_match"]
    &
    ~merged["accession_exact_match"]
)

# ------------------------------------------------------------
# Organism comparison
# ------------------------------------------------------------

merged["organism_exact_match"] = (
    merged["organism_original_master"].notna()
    &
    merged["organism_original_opencity"].notna()
    &
    (
        merged["organism_original_master"]
        ==
        merged["organism_original_opencity"]
    )
)

merged["organism_normalized_match"] = (
    merged["organism_normalized_master"].notna()
    &
    merged["organism_normalized_opencity"].notna()
    &
    (
        merged["organism_normalized_master"]
        ==
        merged["organism_normalized_opencity"]
    )
)

merged["organism_genus_species_match"] = (
    merged["organism_genus_species_master"].notna()
    &
    merged["organism_genus_species_opencity"].notna()
    &
    (
        merged["organism_genus_species_master"]
        ==
        merged["organism_genus_species_opencity"]
    )
)

# ------------------------------------------------------------
# Record classification
# ------------------------------------------------------------

def classify(row):

    in_master = row["in_master"]
    in_opencity = row["in_opencity"]

    if in_master and not in_opencity:
        return "MASTER_ONLY"

    if in_opencity and not in_master:
        return "OPENCITY_ONLY"

    # Both exist
    acc_master = row["accession_normalized_master"]
    acc_oc = row["accession_normalized_opencity"]

    org_master = row["organism_normalized_master"]
    org_oc = row["organism_normalized_opencity"]

    acc_available = (
        pd.notna(acc_master)
        and pd.notna(acc_oc)
    )

    org_available = (
        pd.notna(org_master)
        and pd.notna(org_oc)
    )

    acc_match = (
        acc_available
        and acc_master == acc_oc
    )

    org_match = (
        org_available
        and org_master == org_oc
    )

    if acc_match and org_match:
        if row["accession_exact_match"]:
            return "EXACT_MATCH"
        else:
            return "ACCESSION_VERSION_ONLY"

    if acc_match and not org_match:
        return "ORGANISM_CONFLICT"

    if not acc_match and org_match:
        return "ACCESSION_CONFLICT"

    if not acc_available and org_available:
        return "MASTER_OR_OPENCITY_MISSING_ACCESSION"

    if acc_available and not org_available:
        return "MASTER_OR_OPENCITY_MISSING_ORGANISM"

    if not acc_available and not org_available:
        return "BOTH_ACCESSION_AND_ORGANISM_UNRESOLVED"

    return "TRUE_CONFLICT"


merged["reconciliation_status"] = merged.apply(
    classify,
    axis=1
)

# ------------------------------------------------------------
# Additional conflict labels
# ------------------------------------------------------------

merged["accession_comparison"] = "NOT_COMPARABLE"

merged.loc[
    merged["accession_normalized_match"],
    "accession_comparison"
] = "NORMALIZED_MATCH"

merged.loc[
    merged["accession_exact_match"],
    "accession_comparison"
] = "EXACT_MATCH"

merged.loc[
    merged["in_master"]
    &
    merged["in_opencity"]
    &
    ~merged["accession_normalized_match"]
    &
    merged["accession_normalized_master"].notna()
    &
    merged["accession_normalized_opencity"].notna(),
    "accession_comparison"
] = "REAL_ACCESSION_DIFFERENCE"

# ------------------------------------------------------------
# Organism comparison label
# ------------------------------------------------------------

merged["organism_comparison"] = "NOT_COMPARABLE"

merged.loc[
    merged["organism_normalized_match"],
    "organism_comparison"
] = "NORMALIZED_MATCH"

merged.loc[
    merged["organism_exact_match"],
    "organism_comparison"
] = "EXACT_MATCH"

merged.loc[
    merged["organism_genus_species_match"]
    &
    ~merged["organism_normalized_match"],
    "organism_comparison"
] = "POSSIBLE_SPELLING_OR_FORMAT_DIFFERENCE"

merged.loc[
    merged["in_master"]
    &
    merged["in_opencity"]
    &
    merged["organism_normalized_master"].notna()
    &
    merged["organism_normalized_opencity"].notna()
    &
    ~merged["organism_genus_species_match"],
    "organism_comparison"
] = "POSSIBLE_TRUE_TAXONOMIC_CONFLICT"

# ------------------------------------------------------------
# KJU number
# ------------------------------------------------------------

merged["kju_number"] = pd.to_numeric(
    merged["kju_id"].str.extract(
        r"KJU(\d+)",
        expand=False
    ),
    errors="coerce"
)

# ------------------------------------------------------------
# Sort
# ------------------------------------------------------------

merged = merged.sort_values(
    by=["kju_number", "kju_id"],
    na_position="last"
)

# ============================================================
# 9. WRITE OUTPUTS
# ============================================================

print("\n[9/10] Writing reconciliation outputs")
print("-" * 90)

# ------------------------------------------------------------
# Select useful audit columns
# ------------------------------------------------------------

preferred_columns = [
    "kju_id",

    "reconciliation_status",

    "in_master",
    "in_opencity",

    "accession_original_master",
    "accession_original_opencity",

    "accession_normalized_master",
    "accession_normalized_opencity",

    "accession_exact_match",
    "accession_normalized_match",
    "accession_version_only_difference",

    "accession_comparison",

    "organism_original_master",
    "organism_original_opencity",

    "organism_normalized_master",
    "organism_normalized_opencity",

    "organism_exact_match",
    "organism_normalized_match",
    "organism_genus_species_match",

    "organism_comparison",

    "organism_genus_species_master",
    "organism_genus_species_opencity",

    "_id_opencity",
    "ncbi_url_opencity",

    "kju_id_original_master",
    "kju_id_original_opencity",

    "kju_number",
]

final_columns = [
    c for c in preferred_columns
    if c in merged.columns
]

reconciliation = merged[final_columns].copy()

# ------------------------------------------------------------
# Full reconciliation
# ------------------------------------------------------------

reconciliation.to_csv(
    FULL_OUT,
    index=False
)

# ------------------------------------------------------------
# Matches
# ------------------------------------------------------------

matches = reconciliation[
    reconciliation["reconciliation_status"].isin(
        [
            "EXACT_MATCH",
            "ACCESSION_VERSION_ONLY",
        ]
    )
].copy()

matches.to_csv(
    MATCH_OUT,
    index=False
)

# ------------------------------------------------------------
# Accession version-only
# ------------------------------------------------------------

version_only = reconciliation[
    reconciliation["reconciliation_status"]
    == "ACCESSION_VERSION_ONLY"
].copy()

version_only.to_csv(
    VERSION_OUT,
    index=False
)

# ------------------------------------------------------------
# Organism conflicts
# ------------------------------------------------------------

organism_conflicts = reconciliation[
    reconciliation["reconciliation_status"]
    == "ORGANISM_CONFLICT"
].copy()

organism_conflicts.to_csv(
    ORGANISM_OUT,
    index=False
)

# ------------------------------------------------------------
# True conflicts
# ------------------------------------------------------------

true_conflicts = reconciliation[
    reconciliation["reconciliation_status"]
    == "TRUE_CONFLICT"
].copy()

true_conflicts.to_csv(
    CONFLICT_OUT,
    index=False
)

# ------------------------------------------------------------
# Master only
# ------------------------------------------------------------

master_only = reconciliation[
    reconciliation["reconciliation_status"]
    == "MASTER_ONLY"
].copy()

master_only.to_csv(
    MASTER_ONLY_OUT,
    index=False
)

# ------------------------------------------------------------
# OpenCity only
# ------------------------------------------------------------

opencity_only = reconciliation[
    reconciliation["reconciliation_status"]
    == "OPENCITY_ONLY"
].copy()

opencity_only.to_csv(
    OPENCITY_ONLY_OUT,
    index=False
)

# ============================================================
# 10. SUMMARY
# ============================================================

print("\n[10/10] Generating final summary")
print("-" * 90)

status_counts = (
    reconciliation["reconciliation_status"]
    .value_counts()
    .sort_index()
)

# ------------------------------------------------------------
# Basic statistics
# ------------------------------------------------------------

total_records = len(reconciliation)

master_records = reconciliation[
    reconciliation["in_master"]
].shape[0]

opencity_records = reconciliation[
    reconciliation["in_opencity"]
].shape[0]

both_records = reconciliation[
    reconciliation["in_master"]
    &
    reconciliation["in_opencity"]
].shape[0]

exact_matches = len(
    reconciliation[
        reconciliation["reconciliation_status"]
        == "EXACT_MATCH"
    ]
)

version_matches = len(version_only)

organism_conflict_count = len(
    organism_conflicts
)

true_conflict_count = len(
    true_conflicts
)

master_only_count = len(
    master_only
)

opencity_only_count = len(
    opencity_only
)

real_accession_difference = len(
    reconciliation[
        reconciliation["accession_comparison"]
        == "REAL_ACCESSION_DIFFERENCE"
    ]
)

possible_spelling = len(
    reconciliation[
        reconciliation["organism_comparison"]
        == "POSSIBLE_SPELLING_OR_FORMAT_DIFFERENCE"
    ]
)

possible_taxonomic = len(
    reconciliation[
        reconciliation["organism_comparison"]
        == "POSSIBLE_TRUE_TAXONOMIC_CONFLICT"
    ]
)

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

print("\n")
print("=" * 90)
print(" FINAL RECONCILIATION SUMMARY")
print("=" * 90)

print(
    f"\nOpenCity records                  : "
    f"{opencity_records}"
)

print(
    f"Master records                    : "
    f"{master_records}"
)

print(
    f"Combined unique KJU records       : "
    f"{total_records}"
)

print(
    f"KJU IDs present in both           : "
    f"{both_records}"
)

print("\nStatus breakdown:")
print("-" * 90)

for status, count in status_counts.items():
    print(
        f"{status:<50} {count:>6}"
    )

print("\n")
print("=" * 90)
print(" IMPORTANT RESULTS")
print("=" * 90)

print(
    f"\nExact matches                     : "
    f"{exact_matches}"
)

print(
    f"Accession version-only matches    : "
    f"{version_matches}"
)

print(
    f"Real accession differences        : "
    f"{real_accession_difference}"
)

print(
    f"Organism conflicts                : "
    f"{organism_conflict_count}"
)

print(
    f"Possible organism spelling issues : "
    f"{possible_spelling}"
)

print(
    f"Possible taxonomic conflicts       : "
    f"{possible_taxonomic}"
)

print(
    f"MASTER-only KJU IDs               : "
    f"{master_only_count}"
)

print(
    f"OpenCity-only KJU IDs             : "
    f"{opencity_only_count}"
)

# ============================================================
# PRINT ACCESSION VERSION MATCHES
# ============================================================

if len(version_only) > 0:

    print("\n")
    print("=" * 90)
    print(" ACCESSION VERSION-ONLY DIFFERENCES")
    print("=" * 90)

    print(
        version_only[
            [
                "kju_id",
                "accession_original_master",
                "accession_original_opencity",
                "organism_original_master",
                "organism_original_opencity",
            ]
        ].to_string(index=False)
    )

# ============================================================
# PRINT ORGANISM CONFLICTS
# ============================================================

if len(organism_conflicts) > 0:

    print("\n")
    print("=" * 90)
    print(" ORGANISM CONFLICTS")
    print("=" * 90)

    print(
        organism_conflicts[
            [
                "kju_id",
                "accession_normalized_master",
                "accession_normalized_opencity",
                "organism_original_master",
                "organism_original_opencity",
            ]
        ].to_string(index=False)
    )

# ============================================================
# PRINT TRUE ACCESSION CONFLICTS
# ============================================================

real_acc_conflicts = reconciliation[
    reconciliation["accession_comparison"]
    == "REAL_ACCESSION_DIFFERENCE"
].copy()

if len(real_acc_conflicts) > 0:

    print("\n")
    print("=" * 90)
    print(" REAL ACCESSION DIFFERENCES")
    print("=" * 90)

    print(
        real_acc_conflicts[
            [
                "kju_id",
                "accession_original_master",
                "accession_original_opencity",
                "organism_original_master",
                "organism_original_opencity",
            ]
        ].to_string(index=False)
    )

# ============================================================
# PRINT MASTER ONLY
# ============================================================

if len(master_only) > 0:

    print("\n")
    print("=" * 90)
    print(" KJU IDs PRESENT IN MASTER BUT NOT OPENCITY")
    print("=" * 90)

    print(
        master_only[
            [
                "kju_id",
                "accession_original_master",
                "organism_original_master",
            ]
        ].to_string(index=False)
    )

# ============================================================
# PRINT OPENCITY ONLY
# ============================================================

if len(opencity_only) > 0:

    print("\n")
    print("=" * 90)
    print(" KJU IDs PRESENT IN OPENCITY BUT NOT MASTER")
    print("=" * 90)

    print(
        opencity_only[
            [
                "kju_id",
                "accession_original_opencity",
                "organism_original_opencity",
            ]
        ].to_string(index=False)
    )

# ============================================================
# WRITE SUMMARY TXT
# ============================================================

summary_lines = []

summary_lines.append(
    "DODDAGGUBBI - OPENCITY ↔ MASTER KJU RECONCILIATION v2"
)

summary_lines.append("=" * 80)

summary_lines.append(
    f"OpenCity records: {opencity_records}"
)

summary_lines.append(
    f"Master records: {master_records}"
)

summary_lines.append(
    f"Combined unique KJU records: {total_records}"
)

summary_lines.append(
    f"KJU IDs present in both: {both_records}"
)

summary_lines.append("")
summary_lines.append("STATUS BREAKDOWN")
summary_lines.append("-" * 80)

for status, count in status_counts.items():
    summary_lines.append(
        f"{status}: {count}"
    )

summary_lines.append("")
summary_lines.append("KEY RESULTS")
summary_lines.append("-" * 80)

summary_lines.append(
    f"Exact matches: {exact_matches}"
)

summary_lines.append(
    f"Accession version-only matches: {version_matches}"
)

summary_lines.append(
    f"Real accession differences: {real_accession_difference}"
)

summary_lines.append(
    f"Organism conflicts: {organism_conflict_count}"
)

summary_lines.append(
    f"Possible spelling/format differences: {possible_spelling}"
)

summary_lines.append(
    f"Possible taxonomic conflicts: {possible_taxonomic}"
)

summary_lines.append(
    f"MASTER-only: {master_only_count}"
)

summary_lines.append(
    f"OpenCity-only: {opencity_only_count}"
)

summary_lines.append("")
summary_lines.append("OUTPUT FILES")
summary_lines.append("-" * 80)

for path in [
    FULL_OUT,
    MATCH_OUT,
    VERSION_OUT,
    ORGANISM_OUT,
    CONFLICT_OUT,
    MASTER_ONLY_OUT,
    OPENCITY_ONLY_OUT,
]:
    summary_lines.append(str(path))

SUMMARY_OUT.write_text(
    "\n".join(summary_lines),
    encoding="utf-8"
)

# ============================================================
# FINAL
# ============================================================

print("\n")
print("=" * 90)
print(" RECONCILIATION v2 COMPLETE")
print("=" * 90)

print("\nOutput files:")
print("-" * 90)

for path in [
    FULL_OUT,
    MATCH_OUT,
    VERSION_OUT,
    ORGANISM_OUT,
    CONFLICT_OUT,
    MASTER_ONLY_OUT,
    OPENCITY_ONLY_OUT,
    SUMMARY_OUT,
]:
    print(path)

print("\n")
print("=" * 90)
print(" IMPORTANT")
print("=" * 90)

print(
    """
Accession version suffixes (.1, .2, etc.) were normalized ONLY
for comparison. Original accession values were preserved.

No organism names were automatically corrected.

Possible spelling differences and possible true taxonomic
conflicts are separated for manual review.

MASTER and OpenCity records were not overwritten.
"""
)

print("=" * 90)

