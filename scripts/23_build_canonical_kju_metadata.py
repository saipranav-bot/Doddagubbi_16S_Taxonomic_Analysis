from pathlib import Path
import pandas as pd

# ============================================================
# DODDAGGUBBI
# BUILD CANONICAL KJU METADATA
#
# SAFE / NON-DESTRUCTIVE
#
# MASTER FILE IS NEVER MODIFIED
# ============================================================

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

MASTER_FILE = META / "MASTER_KJU_RECONCILIATION_001_183.csv"
NCBI_RESOLUTION = META / "NCBI_RESOLUTION_TABLE_V1.csv"
KJU0039_RESOLUTION = META / "KJU0039_FINAL_RESOLUTION_V1.csv"

OUT = META / "CANONICAL_KJU_METADATA_V1.csv"
SUMMARY = META / "CANONICAL_KJU_METADATA_SUMMARY_V1.txt"

print("=" * 90)
print(" DODDAGGUBBI - BUILD CANONICAL KJU METADATA")
print("=" * 90)

# ------------------------------------------------------------
# Validate input files
# ------------------------------------------------------------

for f in [MASTER_FILE, NCBI_RESOLUTION, KJU0039_RESOLUTION]:
    if not f.exists():
        raise SystemExit(f"ERROR: Missing required file:\n{f}")

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------

master = pd.read_csv(
    MASTER_FILE,
    dtype=str
).fillna("")

ncbi = pd.read_csv(
    NCBI_RESOLUTION,
    dtype=str
).fillna("")

kju0039 = pd.read_csv(
    KJU0039_RESOLUTION,
    dtype=str
).fillna("")

# ------------------------------------------------------------
# Validate master
# ------------------------------------------------------------

required_master = [
    "kju_id",
    "status",
    "accession",
    "organism",
    "sequence_type",
    "length",
]

missing = [c for c in required_master if c not in master.columns]

if missing:
    raise SystemExit(
        "ERROR: MASTER is missing columns:\n"
        + "\n".join(missing)
    )

if master["kju_id"].duplicated().any():
    raise SystemExit("ERROR: Duplicate KJU IDs found in master.")

if len(master) != 183:
    raise SystemExit(
        f"ERROR: Expected 183 master rows, found {len(master)}"
    )

# ------------------------------------------------------------
# Validate NCBI resolution
# ------------------------------------------------------------

required_ncbi = [
    "kju_id",
    "ncbi_accession",
    "ncbi_organism",
    "master_organism",
    "opencity_organism",
    "master_supported_by_ncbi",
    "opencity_supported_by_ncbi",
    "resolution",
    "classification",
    "recommended_organism",
    "recommended_action",
]

missing = [c for c in required_ncbi if c not in ncbi.columns]

if missing:
    raise SystemExit(
        "ERROR: NCBI resolution is missing columns:\n"
        + "\n".join(missing)
    )

# ------------------------------------------------------------
# Validate KJU0039 resolution
# ------------------------------------------------------------

if len(kju0039) != 1:
    raise SystemExit(
        f"ERROR: Expected exactly one KJU0039 resolution row, "
        f"found {len(kju0039)}"
    )

k39 = kju0039.iloc[0]

if k39["kju_id"] != "KJU0039":
    raise SystemExit("ERROR: KJU0039 resolution file has wrong KJU ID.")

# ------------------------------------------------------------
# Prepare resolution lookup
# ------------------------------------------------------------

ncbi_lookup = {}

for _, row in ncbi.iterrows():

    kju = row["kju_id"].strip()

    if not kju:
        continue

    if kju in ncbi_lookup:
        raise SystemExit(
            f"ERROR: Multiple NCBI resolution rows for {kju}"
        )

    ncbi_lookup[kju] = row.to_dict()

# ------------------------------------------------------------
# Build canonical records
# ------------------------------------------------------------

records = []

for _, m in master.iterrows():

    kju = m["kju_id"].strip()

    n = ncbi_lookup.get(kju, {})

    master_accession = m["accession"].strip()
    master_organism = m["organism"].strip()
    master_sequence_type = m["sequence_type"].strip()
    master_length = m["length"].strip()

    ncbi_accession = n.get("ncbi_accession", "").strip()
    ncbi_organism = n.get("ncbi_organism", "").strip()

    ncbi_resolution = n.get("resolution", "").strip()
    ncbi_classification = n.get("classification", "").strip()

    recommended_organism = n.get(
        "recommended_organism", ""
    ).strip()

    recommended_action = n.get(
        "recommended_action", ""
    ).strip()

    master_supported = n.get(
        "master_supported_by_ncbi", ""
    ).strip()

    opencity_supported = n.get(
        "opencity_supported_by_ncbi", ""
    ).strip()

    opencity_organism = n.get(
        "opencity_organism", ""
    ).strip()

    # --------------------------------------------------------
    # Default canonical values
    # --------------------------------------------------------

    canonical_accession = master_accession
    canonical_organism = master_organism
    canonical_sequence_type = master_sequence_type
    canonical_length = master_length

    provenance_source = "MASTER"

    evidence_level = "MASTER_BASELINE"

    resolution_status = "NO_NCBI_DISPUTE"

    # --------------------------------------------------------
    # Apply only defensible NCBI-supported organism correction
    # --------------------------------------------------------

    if (
        recommended_organism
        and recommended_organism != master_organism
        and ncbi_classification == "MASTER_SUPPORTED"
    ):

        canonical_organism = recommended_organism

        provenance_source = "MASTER+NCBI"

        evidence_level = "NCBI_VALIDATED"

        resolution_status = "NCBI_ORGANISM_VALIDATED"

    elif ncbi_classification == "MASTER_SUPPORTED":

        provenance_source = "MASTER+NCBI"

        evidence_level = "NCBI_VALIDATED"

        resolution_status = "NCBI_MASTER_SUPPORTED"

    elif ncbi_classification == "OPENCITY_SUPPORTED":

        # ----------------------------------------------------
        # OpenCity-supported records without a master row
        # are NOT allowed to overwrite the master.
        # ----------------------------------------------------

        provenance_source = "MASTER"

        evidence_level = "MASTER_BASELINE"

        resolution_status = "OPENCITY_ONLY_NOT_MASTER"

    elif ncbi_classification == "SPELLING_OR_TAXONOMIC_NAME_STANDARDIZATION":

        # ----------------------------------------------------
        # Taxonomic spelling correction is represented in the
        # canonical layer only when a master organism exists.
        # ----------------------------------------------------

        if master_organism and recommended_organism:

            canonical_organism = recommended_organism

            provenance_source = "MASTER+NCBI_TAXONOMY"

            evidence_level = "NCBI_TAXONOMIC_STANDARDIZATION"

            resolution_status = "TAXONOMIC_NAME_STANDARDIZED"

        else:

            provenance_source = "MASTER"

            evidence_level = "MASTER_BASELINE"

            resolution_status = "TAXONOMIC_REVIEW_NO_MASTER"

    elif ncbi_classification == "ACCESSION_KJU_COLLISION":

        resolution_status = "UNRESOLVED_ACCESSION_COLLISION"

        evidence_level = "UNRESOLVED"

        provenance_source = "MASTER+NCBI_REVIEW"

    # --------------------------------------------------------
    # SPECIAL KJU0039 SAFETY OVERRIDE
    # --------------------------------------------------------

    if kju == "KJU0039":

        canonical_accession = ""

        canonical_organism = ""

        canonical_sequence_type = ""

        canonical_length = ""

        provenance_source = "UNRESOLVED"

        evidence_level = "UNRESOLVED"

        resolution_status = k39["resolution_status"].strip()

        recommended_action = k39["action"].strip()

    # --------------------------------------------------------
    # Record
    # --------------------------------------------------------

    records.append({

        "kju_id": kju,

        "master_status": m["status"].strip(),

        "canonical_accession": canonical_accession,
        "canonical_organism": canonical_organism,
        "canonical_sequence_type": canonical_sequence_type,
        "canonical_length": canonical_length,

        "master_accession": master_accession,
        "master_organism": master_organism,
        "master_sequence_type": master_sequence_type,
        "master_length": master_length,

        "ncbi_accession": ncbi_accession,
        "ncbi_organism": ncbi_organism,

        "opencity_organism": opencity_organism,

        "ncbi_resolution": ncbi_resolution,
        "ncbi_classification": ncbi_classification,

        "master_supported_by_ncbi": master_supported,
        "opencity_supported_by_ncbi": opencity_supported,

        "provenance_source": provenance_source,
        "evidence_level": evidence_level,
        "resolution_status": resolution_status,

        "recommended_action": recommended_action,
    })

# ------------------------------------------------------------
# Create dataframe
# ------------------------------------------------------------

canonical = pd.DataFrame(records)

# ------------------------------------------------------------
# Final structural checks
# ------------------------------------------------------------

if len(canonical) != 183:
    raise SystemExit(
        f"ERROR: Canonical table has {len(canonical)} rows; "
        "expected 183."
    )

if canonical["kju_id"].duplicated().any():
    raise SystemExit(
        "ERROR: Duplicate KJU IDs in canonical table."
    )

if canonical["kju_id"].nunique() != 183:
    raise SystemExit(
        "ERROR: Canonical KJU ID count is not 183."
    )

# ------------------------------------------------------------
# KJU0038 safety
# ------------------------------------------------------------

r38 = canonical[
    canonical["kju_id"] == "KJU0038"
]

if len(r38) != 1:
    raise SystemExit("ERROR: KJU0038 missing from canonical table.")

r38 = r38.iloc[0]

if r38["canonical_accession"] != "PZ666011.1":
    raise SystemExit(
        "ERROR: KJU0038 canonical accession is incorrect."
    )

if r38["canonical_organism"] != "Shewanella algae":
    raise SystemExit(
        "ERROR: KJU0038 canonical organism is incorrect."
    )

# ------------------------------------------------------------
# KJU0039 safety
# ------------------------------------------------------------

r39 = canonical[
    canonical["kju_id"] == "KJU0039"
]

if len(r39) != 1:
    raise SystemExit("ERROR: KJU0039 missing from canonical table.")

r39 = r39.iloc[0]

if r39["canonical_accession"] != "":
    raise SystemExit(
        "SAFETY ERROR: KJU0039 received an accession."
    )

if r39["canonical_organism"] != "":
    raise SystemExit(
        "SAFETY ERROR: KJU0039 received an organism."
    )

if r39["resolution_status"] != "UNRESOLVED_ACCESSION_COLLISION":
    raise SystemExit(
        "SAFETY ERROR: KJU0039 resolution status changed."
    )

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

canonical.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary_lines = []

summary_lines.append(
    "DODDAGGUBBI - CANONICAL KJU METADATA SUMMARY"
)

summary_lines.append("=" * 80)

summary_lines.append("")

summary_lines.append(
    f"Canonical records: {len(canonical)}"
)

summary_lines.append(
    f"Unique KJU IDs: {canonical['kju_id'].nunique()}"
)

summary_lines.append("")

summary_lines.append("EVIDENCE LEVEL")

summary_lines.append("-" * 80)

for key, value in canonical["evidence_level"].value_counts().items():
    summary_lines.append(
        f"{key}: {value}"
    )

summary_lines.append("")

summary_lines.append("RESOLUTION STATUS")

summary_lines.append("-" * 80)

for key, value in canonical["resolution_status"].value_counts().items():
    summary_lines.append(
        f"{key}: {value}"
    )

summary_lines.append("")

summary_lines.append("CANONICAL ACCESSION COVERAGE")

summary_lines.append("-" * 80)

with_accession = (
    canonical["canonical_accession"].str.strip() != ""
).sum()

without_accession = len(canonical) - with_accession

summary_lines.append(
    f"Records with canonical accession: {with_accession}"
)

summary_lines.append(
    f"Records without canonical accession: {without_accession}"
)

summary_lines.append("")

summary_lines.append("KJU0038")

summary_lines.append("-" * 80)

summary_lines.append(
    "KJU0038 canonical accession: "
    + r38["canonical_accession"]
)

summary_lines.append(
    "KJU0038 canonical organism: "
    + r38["canonical_organism"]
)

summary_lines.append("")

summary_lines.append("KJU0039")

summary_lines.append("-" * 80)

summary_lines.append(
    "KJU0039 canonical accession: [BLANK]"
)

summary_lines.append(
    "KJU0039 canonical organism: [BLANK]"
)

summary_lines.append(
    "KJU0039 status: UNRESOLVED_ACCESSION_COLLISION"
)

summary_lines.append(
    "KJU0039 action: DO NOT MERGE"
)

summary_lines.append("")

summary_lines.append("IMPORTANT")

summary_lines.append("-" * 80)

summary_lines.append(
    "MASTER_KJU_RECONCILIATION_001_183.csv was NOT modified."
)

summary_lines.append(
    "This canonical table is a downstream derived layer."
)

summary_lines.append(
    "NCBI and OpenCity evidence are retained separately."
)

summary_lines.append(
    "KJU0039 was NOT assigned PZ666011.1."
)

SUMMARY.write_text(
    "\n".join(summary_lines) + "\n",
    encoding="utf-8"
)

print("")
print("=" * 90)
print(" CANONICAL METADATA BUILD COMPLETE")
print("=" * 90)

print(f"Output : {OUT}")
print(f"Summary: {SUMMARY}")

print("")
print("Records:", len(canonical))

print("")
print("Evidence levels:")
print(
    canonical["evidence_level"]
    .value_counts()
    .to_string()
)

print("")
print("KJU0038:")
print(
    canonical[
        canonical["kju_id"] == "KJU0038"
    ].to_string(index=False)
)

print("")
print("KJU0039:")
print(
    canonical[
        canonical["kju_id"] == "KJU0039"
    ].to_string(index=False)
)

print("")
print("MASTER DATA WAS NOT MODIFIED.")
