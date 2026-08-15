from pathlib import Path
import pandas as pd

# ============================================================
# DODDAGGUBBI
# BUILD NCBI-BASED RESOLUTION TABLE
# ============================================================

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

NCBI_FILE = META / "NCBI_VALIDATION_DISPUTED_RECORDS.csv"

OUT = META / "NCBI_RESOLUTION_TABLE_V1.csv"
SUMMARY = META / "NCBI_RESOLUTION_SUMMARY_V1.txt"

print("=" * 90)
print(" DODDAGGUBBI - BUILD NCBI RESOLUTION TABLE")
print("=" * 90)

if not NCBI_FILE.exists():
    raise SystemExit(f"ERROR: Missing file:\n{NCBI_FILE}")

df = pd.read_csv(NCBI_FILE, dtype=str).fillna("")

# ------------------------------------------------------------
# Collapse duplicate accession checks
# ------------------------------------------------------------

records = []

for kju, g in df.groupby("kju_id", sort=True):

    # NCBI should normally return one biological identity
    ncbi_orgs = sorted(set(
        x.strip() for x in g["ncbi_organism"]
        if x.strip()
    ))

    ncbi_accessions = sorted(set(
        x.strip() for x in g["ncbi_accession"]
        if x.strip()
    ))

    master_orgs = sorted(set(
        x.strip() for x in g["organism_original_master"]
        if x.strip()
    ))

    opencity_orgs = sorted(set(
        x.strip() for x in g["organism_original_opencity"]
        if x.strip()
    ))

    master_supported = (
        g["master_matches_ncbi"].astype(str).str.lower() == "true"
    ).any()

    opencity_supported = (
        g["opencity_matches_ncbi"].astype(str).str.lower() == "true"
    ).any()

    # --------------------------------------------------------
    # Determine resolution
    # --------------------------------------------------------

    if kju == "KJU0039":
        resolution = "ACCESSION_IDENTITY_CONFLICT_REVIEW"

    elif master_supported and not opencity_supported:
        resolution = "USE_MASTER_NCBI_VALIDATED"

    elif opencity_supported and not master_supported:
        resolution = "USE_OPENCITY_NCBI_VALIDATED"

    elif master_supported and opencity_supported:
        resolution = "BOTH_SOURCES_NCBI_SUPPORTED"

    else:
        resolution = "USE_NCBI_TAXONOMY_SPELLING_REVIEW"

    # --------------------------------------------------------
    # Classification
    # --------------------------------------------------------

    if kju == "KJU0039":
        classification = "ACCESSION_KJU_COLLISION"

    elif resolution == "USE_MASTER_NCBI_VALIDATED":
        classification = "MASTER_SUPPORTED"

    elif resolution == "USE_OPENCITY_NCBI_VALIDATED":
        classification = "OPENCITY_SUPPORTED"

    elif resolution == "BOTH_SOURCES_NCBI_SUPPORTED":
        classification = "NO_BIOLOGICAL_CONFLICT"

    else:
        classification = "SPELLING_OR_TAXONOMIC_NAME_STANDARDIZATION"

    records.append({
        "kju_id": kju,
        "ncbi_accession": ";".join(ncbi_accessions),
        "ncbi_organism": ";".join(ncbi_orgs),
        "master_organism": ";".join(master_orgs),
        "opencity_organism": ";".join(opencity_orgs),
        "master_supported_by_ncbi": master_supported,
        "opencity_supported_by_ncbi": opencity_supported,
        "resolution": resolution,
        "classification": classification,
    })

out = pd.DataFrame(records)

# ------------------------------------------------------------
# Add recommended organism
# ------------------------------------------------------------

def choose_organism(row):

    if row["classification"] == "ACCESSION_KJU_COLLISION":
        return ""

    if row["resolution"] == "USE_MASTER_NCBI_VALIDATED":
        return row["ncbi_organism"]

    if row["resolution"] == "USE_OPENCITY_NCBI_VALIDATED":
        return row["ncbi_organism"]

    if row["resolution"] == "USE_NCBI_TAXONOMY_SPELLING_REVIEW":
        return row["ncbi_organism"]

    if row["resolution"] == "BOTH_SOURCES_NCBI_SUPPORTED":
        return row["ncbi_organism"]

    return ""

out["recommended_organism"] = out.apply(choose_organism, axis=1)

# ------------------------------------------------------------
# Recommended action
# ------------------------------------------------------------

def action(row):

    if row["classification"] == "ACCESSION_KJU_COLLISION":
        return (
            "DO NOT MERGE; investigate KJU0039 because "
            "PZ666011 belongs to KJU0038 in NCBI"
        )

    if row["resolution"] == "USE_MASTER_NCBI_VALIDATED":
        return (
            "RETAIN MASTER ORGANISM; mark OpenCity value as "
            "incorrect/variant"
        )

    if row["resolution"] == "USE_OPENCITY_NCBI_VALIDATED":
        return (
            "RETAIN OPENCITY ORGANISM; NCBI validates OpenCity"
        )

    if row["resolution"] == "USE_NCBI_TAXONOMY_SPELLING_REVIEW":
        return (
            "STANDARDIZE TO NCBI ORGANISM NAME; retain original "
            "source values for provenance"
        )

    return "NO CHANGE REQUIRED"

out["recommended_action"] = out.apply(action, axis=1)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

out.to_csv(OUT, index=False)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

lines = []

lines.append("DODDAGGUBBI - NCBI RESOLUTION SUMMARY")
lines.append("=" * 80)
lines.append("")
lines.append(f"Records reviewed: {len(out)}")
lines.append("")

lines.append("CLASSIFICATION")
lines.append("-" * 80)

for key, value in out["classification"].value_counts().items():
    lines.append(f"{key}: {value}")

lines.append("")
lines.append("RESOLUTION")
lines.append("-" * 80)

for key, value in out["resolution"].value_counts().items():
    lines.append(f"{key}: {value}")

lines.append("")
lines.append("ACCESSION/KJU COLLISIONS")
lines.append("-" * 80)

collisions = out[
    out["classification"] == "ACCESSION_KJU_COLLISION"
]

if len(collisions) == 0:
    lines.append("None")

else:
    for _, r in collisions.iterrows():
        lines.append(
            f"{r['kju_id']} -> {r['ncbi_accession']} -> "
            f"{r['ncbi_organism']}"
        )

lines.append("")
lines.append("RECOMMENDED ORGANISM CHANGES")
lines.append("-" * 80)

for _, r in out.iterrows():

    if r["recommended_organism"]:
        lines.append(
            f"{r['kju_id']}: "
            f"{r['master_organism']} | "
            f"{r['opencity_organism']} -> "
            f"{r['recommended_organism']}"
        )

lines.append("")
lines.append("IMPORTANT")
lines.append("-" * 80)
lines.append(
    "This table is a resolution/review layer only."
)
lines.append(
    "MASTER and OpenCity source files were NOT modified."
)

SUMMARY.write_text("\n".join(lines), encoding="utf-8")

print("")
print("NCBI RESOLUTION COMPLETE")
print("")
print(f"Output table:")
print(OUT)
print("")
print(f"Summary:")
print(SUMMARY)
print("")
print("MASTER DATA WAS NOT MODIFIED.")
