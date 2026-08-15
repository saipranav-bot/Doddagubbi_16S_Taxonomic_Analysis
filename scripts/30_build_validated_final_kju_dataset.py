from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

EXTENDED = META / "MASTER_KJU_16S_DATASET.csv"
DIRECT = META / "direct_genbank_kju_metadata.csv"
ORPHAN = META / "ORPHAN_NCBI_RECONCILIATION_V1.csv"
CANONICAL = META / "MASTER_KJU_RECONCILIATION_001_183.csv"

OUT = META / "FINAL_KJU_16S_DATASET_VALIDATED_V1.csv"
UNRESOLVED = META / "KJU_UNRESOLVED_METADATA_V1.csv"
AUDIT = META / "FINAL_KJU_DATASET_VALIDATION_AUDIT_V1.csv"
SUMMARY = META / "FINAL_KJU_DATASET_VALIDATION_SUMMARY_V1.txt"

print("=" * 110)
print(" DODDAGGUBBI - VALIDATED FINAL KJU 16S DATASET")
print("=" * 110)


def load(path):
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str).fillna("")


def norm_kju(x):
    m = re.search(r"KJU\d{4}", str(x).upper())
    return m.group(0) if m else ""


def clean(x):
    return str(x).strip()


extended = load(EXTENDED)
direct = load(DIRECT)
orphan = load(ORPHAN)
canonical = load(CANONICAL)

print("\nSOURCE COUNTS")
print("-" * 110)
print("Extended dataset :", len(extended))
print("Direct GenBank   :", len(direct))
print("Orphan table     :", len(orphan))
print("Canonical IDs    :", len(canonical))


# ------------------------------------------------------------------
# Normalize KJU IDs
# ------------------------------------------------------------------

for df in [extended, direct, orphan, canonical]:
    if not df.empty and "kju_id" in df.columns:
        df["_kju"] = df["kju_id"].map(norm_kju)


# ------------------------------------------------------------------
# Build evidence priority
#
# Priority:
#   1. Direct GenBank metadata
#   2. Extended 16S dataset
#   3. Orphan NCBI metadata
#
# We NEVER use the incomplete canonical table to invent metadata.
# ------------------------------------------------------------------

records = {}


def add_record(kju, accession, organism, source, row):
    kju = norm_kju(kju)
    accession = clean(accession)
    organism = clean(organism)

    if not kju:
        return

    # Critical metadata must exist.
    if not accession or not organism:
        return

    # Existing record: retain higher-priority source.
    if kju not in records:
        records[kju] = {
            "kju_id": kju,
            "accession": accession,
            "organism": organism,
            "sequence_type": clean(row.get("sequence_type", "")),
            "length": clean(row.get("length", "")),
            "molecule_type": clean(row.get("molecule_type", "")),
            "description": clean(row.get("description", "")),
            "strain": clean(row.get("strain", "")),
            "isolation_source": clean(row.get("isolation_source", "")),
            "geo_loc_name": clean(row.get("geo_loc_name", "")),
            "collection_date": clean(row.get("collection_date", "")),
            "collected_by": clean(row.get("collected_by", "")),
            "genbank_file": clean(row.get("genbank_file", "")),
            "metadata_source": source,
        }


# ------------------------------------------------------------------
# 1. DIRECT GENBANK
# ------------------------------------------------------------------

if not direct.empty:

    for _, r in direct.iterrows():

        add_record(
            r.get("_kju", ""),
            r.get("accession", ""),
            r.get("organism", ""),
            "direct_genbank_kju_metadata.csv",
            r
        )


# ------------------------------------------------------------------
# 2. EXTENDED DATASET
# ------------------------------------------------------------------

if not extended.empty:

    for _, r in extended.iterrows():

        kju = r.get("_kju", "")

        accession = (
            r.get("accession", "")
            or r.get("ncbi_accession", "")
        )

        organism = (
            r.get("organism", "")
            or r.get("organism_name", "")
        )

        # If not already present, add it.
        if norm_kju(kju) not in records:

            add_record(
                kju,
                accession,
                organism,
                "MASTER_KJU_16S_DATASET.csv",
                r
            )


# ------------------------------------------------------------------
# 3. ORPHAN NCBI TABLE
# ------------------------------------------------------------------

if not orphan.empty:

    for _, r in orphan.iterrows():

        kju = r.get("_kju", "")

        accession = (
            r.get("ncbi_accession", "")
            or r.get("accession", "")
        )

        organism = (
            r.get("organism_in_file", "")
            or r.get("organism", "")
        )

        if norm_kju(kju) not in records:

            add_record(
                kju,
                accession,
                organism,
                "ORPHAN_NCBI_RECONCILIATION_V1.csv",
                r
            )


# ------------------------------------------------------------------
# Create final dataframe
# ------------------------------------------------------------------

final = pd.DataFrame(list(records.values()))

if final.empty:
    raise RuntimeError("No valid records were found.")


# ------------------------------------------------------------------
# Normalize columns
# ------------------------------------------------------------------

required_columns = [
    "kju_id",
    "accession",
    "organism",
    "sequence_type",
    "length",
    "molecule_type",
    "description",
    "strain",
    "isolation_source",
    "geo_loc_name",
    "collection_date",
    "collected_by",
    "genbank_file",
    "metadata_source",
]

for c in required_columns:
    if c not in final.columns:
        final[c] = ""


final = final[required_columns].copy()


# ------------------------------------------------------------------
# Sort numerically by KJU
# ------------------------------------------------------------------

final["_num"] = (
    final["kju_id"]
    .str.extract(r"KJU(\d{4})", expand=False)
    .astype(int)
)

final = final.sort_values("_num").drop(columns="_num")


# ------------------------------------------------------------------
# Strict validation
# ------------------------------------------------------------------

errors = []

if final["kju_id"].duplicated().any():
    errors.append("Duplicate KJU IDs")

if final["accession"].duplicated().any():
    duplicated_accessions = final.loc[
        final["accession"].duplicated(keep=False),
        "accession"
    ].tolist()

    errors.append(
        "Duplicate accessions: " +
        ", ".join(duplicated_accessions)
    )

if (final["kju_id"].str.strip() == "").any():
    errors.append("Missing KJU ID")

if (final["accession"].str.strip() == "").any():
    errors.append("Missing accession")

if (final["organism"].str.strip() == "").any():
    errors.append("Missing organism")

if final["kju_id"].duplicated().any():
    errors.append("KJU ID duplication after sorting")


# ------------------------------------------------------------------
# Determine unresolved canonical KJUs
# ------------------------------------------------------------------

canonical_ids = set()

if not canonical.empty:

    canonical_ids = set(
        canonical["_kju"]
        .dropna()
        .astype(str)
        .str.strip()
    )

final_ids = set(final["kju_id"])

unresolved_ids = sorted(
    canonical_ids - final_ids,
    key=lambda x: int(x[3:])
)

unresolved_rows = []

for kju in unresolved_ids:

    cr = canonical[
        canonical["_kju"] == kju
    ]

    accession = ""
    organism = ""

    if not cr.empty:

        accession = clean(
            cr.iloc[0].get("accession", "")
        )

        organism = clean(
            cr.iloc[0].get("organism", "")
        )

    unresolved_rows.append({
        "kju_id": kju,
        "canonical_accession": accession,
        "canonical_organism": organism,
        "resolution_status": "UNRESOLVED_METADATA",
    })


unresolved = pd.DataFrame(unresolved_rows)


# ------------------------------------------------------------------
# Audit
# ------------------------------------------------------------------

audit_rows = []

for _, r in final.iterrows():

    audit_rows.append({
        "kju_id": r["kju_id"],
        "accession": r["accession"],
        "organism": r["organism"],
        "metadata_source": r["metadata_source"],
        "valid_kju": bool(re.fullmatch(r"KJU\d{4}", r["kju_id"])),
        "valid_accession": bool(r["accession"].strip()),
        "valid_organism": bool(r["organism"].strip()),
        "status": "VALID"
    })

audit_df = pd.DataFrame(audit_rows)


# ------------------------------------------------------------------
# HARD STOP
# ------------------------------------------------------------------

if errors:

    print("\n" + "=" * 110)
    print("VALIDATION FAILED")
    print("=" * 110)

    for e in errors:
        print("ERROR:", e)

    print("\nNO FINAL DATASET WRITTEN.")
    print("NO MASTER DATA MODIFIED.")

    raise SystemExit(1)


# ------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------

final.to_csv(OUT, index=False)
audit_df.to_csv(AUDIT, index=False)

if unresolved.empty:

    pd.DataFrame(
        columns=[
            "kju_id",
            "canonical_accession",
            "canonical_organism",
            "resolution_status"
        ]
    ).to_csv(UNRESOLVED, index=False)

else:

    unresolved.to_csv(UNRESOLVED, index=False)


# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

summary = []

summary.append(
    "DODDAGGUBBI - VALIDATED FINAL KJU DATASET"
)

summary.append("=" * 110)
summary.append("")

summary.append(
    f"Final validated records: {len(final)}"
)

summary.append(
    f"Unique KJU IDs: {final['kju_id'].nunique()}"
)

summary.append(
    f"Unique accessions: {final['accession'].nunique()}"
)

summary.append("")

summary.append("VALIDATION")
summary.append("-" * 110)

summary.append(
    f"Missing KJU IDs: {(final['kju_id'].str.strip() == '').sum()}"
)

summary.append(
    f"Missing accessions: {(final['accession'].str.strip() == '').sum()}"
)

summary.append(
    f"Missing organisms: {(final['organism'].str.strip() == '').sum()}"
)

summary.append(
    f"Duplicate KJU IDs: {final['kju_id'].duplicated().sum()}"
)

summary.append(
    f"Duplicate accessions: {final['accession'].duplicated().sum()}"
)

summary.append("")

summary.append("SOURCE DISTRIBUTION")
summary.append("-" * 110)

for source, count in final["metadata_source"].value_counts().items():

    summary.append(
        f"{source}: {count}"
    )

summary.append("")

summary.append("CANONICAL COVERAGE")
summary.append("-" * 110)

summary.append(
    f"Canonical KJUs: {len(canonical_ids)}"
)

summary.append(
    f"Validated KJUs: {len(final_ids)}"
)

summary.append(
    f"Unresolved canonical KJUs: {len(unresolved_ids)}"
)

summary.append("")

if unresolved_ids:

    summary.append("UNRESOLVED KJUs")
    summary.append("-" * 110)

    summary.append(
        ", ".join(unresolved_ids)
    )

    summary.append("")

summary.append(
    "IMPORTANT: No metadata was fabricated."
)

summary.append(
    "Records enter the final dataset only when KJU ID, accession and organism are available."
)

summary.append(
    "The incomplete canonical reconciliation table was NOT used to invent missing metadata."
)

summary.append("")

summary.append(
    f"Final dataset: {OUT}"
)

summary.append(
    f"Validation audit: {AUDIT}"
)

summary.append(
    f"Unresolved records: {UNRESOLVED}"
)

SUMMARY.write_text(
    "\n".join(summary),
    encoding="utf-8"
)


# ------------------------------------------------------------------
# Final console report
# ------------------------------------------------------------------

print("\n" + "=" * 110)
print("FINAL VALIDATION PASSED")
print("=" * 110)

print("Final rows              :", len(final))
print("Unique KJU IDs          :", final["kju_id"].nunique())
print("Unique accessions       :", final["accession"].nunique())
print("Missing KJU IDs         :", (final["kju_id"].str.strip() == "").sum())
print("Missing accessions      :", (final["accession"].str.strip() == "").sum())
print("Missing organisms       :", (final["organism"].str.strip() == "").sum())
print("Duplicate KJU IDs       :", final["kju_id"].duplicated().sum())
print("Duplicate accessions    :", final["accession"].duplicated().sum())
print("Unresolved canonical    :", len(unresolved_ids))

print("\n" + "=" * 110)
print("OUTPUT FILES")
print("=" * 110)

print("FINAL     :", OUT)
print("AUDIT     :", AUDIT)
print("UNRESOLVED:", UNRESOLVED)
print("SUMMARY   :", SUMMARY)

print("\nMASTER DATA WAS NOT MODIFIED.")
print("=" * 110)
