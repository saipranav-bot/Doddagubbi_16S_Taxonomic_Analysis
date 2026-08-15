from pathlib import Path
import pandas as pd
import re
import sys

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

CANONICAL = META / "MASTER_KJU_RECONCILIATION_001_183.csv"
DIRECT = META / "direct_genbank_kju_metadata.csv"
EXTENDED = META / "MASTER_KJU_16S_DATASET.csv"
ORPHAN = META / "ORPHAN_NCBI_RECONCILIATION_V1.csv"

OUT = META / "FINAL_KJU_16S_DATASET_V1.csv"
AUDIT = META / "FINAL_KJU_16S_DATASET_AUDIT_V1.csv"
SUMMARY = META / "FINAL_KJU_16S_DATASET_SUMMARY_V1.txt"

print("=" * 100)
print(" DODDAGGUBBI - FINAL KJU DATASET BUILDER")
print("=" * 100)

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------

canonical = pd.read_csv(CANONICAL, dtype=str).fillna("")
direct = pd.read_csv(DIRECT, dtype=str).fillna("")
extended = pd.read_csv(EXTENDED, dtype=str).fillna("")
orphan = pd.read_csv(ORPHAN, dtype=str).fillna("")

# ------------------------------------------------------------
# Normalize
# ------------------------------------------------------------

def norm(x):
    return str(x).strip()

def norm_kju(x):
    x = norm(x).upper()
    m = re.search(r"KJU\d{4}", x)
    return m.group(0) if m else ""

def norm_acc(x):
    x = norm(x).upper()
    if not x:
        return ""
    return x.split()[0]

def norm_org(x):
    x = norm(x).lower()
    x = re.sub(r"\s+", " ", x)
    return x

for df in [canonical, direct, extended, orphan]:

    if "kju_id" in df.columns:
        df["kju_norm"] = df["kju_id"].map(norm_kju)

    if "accession" in df.columns:
        df["acc_norm"] = df["accession"].map(norm_acc)

    if "ncbi_accession" in df.columns:
        df["ncbi_acc_norm"] = df["ncbi_accession"].map(norm_acc)

    if "organism" in df.columns:
        df["org_norm"] = df["organism"].map(norm_org)

# ------------------------------------------------------------
# Original canonical records
# ------------------------------------------------------------

canonical_records = {}

for _, r in canonical.iterrows():

    kju = r["kju_norm"]

    if not kju:
        continue

    canonical_records[kju] = {
        "kju_id": kju,
        "accession": norm_acc(r.get("accession", "")),
        "organism": norm(r.get("organism", "")),
        "sequence_type": norm(r.get("sequence_type", "")),
        "length": norm(r.get("length", "")),
        "source": "ORIGINAL_CANONICAL"
    }

# ------------------------------------------------------------
# Direct GenBank records
# ------------------------------------------------------------

direct_records = {}

for _, r in direct.iterrows():

    kju = r["kju_norm"]

    if not kju:
        continue

    acc = norm_acc(r.get("accession", ""))

    # Keep first complete record, but audit duplicates
    if kju not in direct_records:
        direct_records[kju] = {
            "kju_id": kju,
            "accession": acc,
            "organism": norm(r.get("organism", "")),
            "strain": norm(r.get("strain", "")),
            "sequence_type": norm(r.get("sequence_type", "")),
            "length": norm(r.get("length", "")),
            "isolation_source": norm(r.get("isolation_source", "")),
            "geo_loc_name": norm(r.get("geo_loc_name", "")),
            "collection_date": norm(r.get("collection_date", "")),
            "collected_by": norm(r.get("collected_by", "")),
            "genbank_file": norm(r.get("genbank_file", "")),
        }

# ------------------------------------------------------------
# Orphan NCBI records
# ------------------------------------------------------------

orphan_records = {}

for _, r in orphan.iterrows():

    kju = norm_kju(r.get("kju_id", ""))

    if not kju:
        continue

    orphan_records[kju] = {
        "accession": norm_acc(
            r.get("ncbi_accession",
                  r.get("accession", ""))
        ),
        "organism": norm(
            r.get("organism_in_file",
                  r.get("organism", ""))
        )
    }

# ------------------------------------------------------------
# Candidate set
#
# ORIGINAL CANONICAL
# +
# 19 NCBI-backed orphan KJUs
# ------------------------------------------------------------

candidate_ids = set(canonical_records)

candidate_ids.update(orphan_records)

candidate_ids = sorted(
    candidate_ids,
    key=lambda x: int(x[3:])
)

print()
print("Original canonical records :", len(canonical_records))
print("NCBI orphan candidates    :", len(orphan_records))
print("Candidate final records   :", len(candidate_ids))

# ------------------------------------------------------------
# Build final records
# ------------------------------------------------------------

final_rows = []
audit_rows = []

fatal_errors = []

for kju in candidate_ids:

    c = canonical_records.get(kju, {})
    d = direct_records.get(kju, {})
    o = orphan_records.get(kju, {})

    is_original = bool(c)
    is_orphan = bool(o)

    # --------------------------------------------------------
    # Determine accession
    # --------------------------------------------------------

    acc_c = norm_acc(c.get("accession", ""))
    acc_d = norm_acc(d.get("accession", ""))
    acc_o = norm_acc(o.get("accession", ""))

    accession_values = {
        x for x in [acc_c, acc_d, acc_o] if x
    }

    accession_conflict = len(accession_values) > 1

    # Prefer direct GenBank, then orphan, then canonical
    accession = acc_d or acc_o or acc_c

    # --------------------------------------------------------
    # Determine organism
    # --------------------------------------------------------

    org_c = norm(c.get("organism", ""))
    org_d = norm(d.get("organism", ""))
    org_o = norm(o.get("organism", ""))

    organism_values = {
        norm_org(x)
        for x in [org_c, org_d, org_o]
        if x
    }

    organism_conflict = len(organism_values) > 1

    organism = org_d or org_o or org_c

    # --------------------------------------------------------
    # Source classification
    # --------------------------------------------------------

    if is_original:
        source_class = "ORIGINAL_CANONICAL"
    elif is_orphan and d:
        source_class = "NCBI_ORPHAN_CONFIRMED_BY_DIRECT_GENBANK"
    elif is_orphan:
        source_class = "NCBI_ORPHAN"
    else:
        source_class = "UNKNOWN"

    # --------------------------------------------------------
    # Evidence
    # --------------------------------------------------------

    direct_genbank_evidence = bool(d)
    orphan_ncbi_evidence = bool(o)
    canonical_evidence = bool(c)

    # --------------------------------------------------------
    # Final status
    # --------------------------------------------------------

    if accession_conflict:
        final_status = "REVIEW_ACCESSION_CONFLICT"
        fatal_errors.append(
            f"{kju}: accession conflict {sorted(accession_values)}"
        )

    elif organism_conflict:
        final_status = "REVIEW_ORGANISM_CONFLICT"
        fatal_errors.append(
            f"{kju}: organism conflict {sorted(organism_values)}"
        )

    elif is_original:
        final_status = "CANONICAL"
    elif is_orphan and direct_genbank_evidence:
        final_status = "CANONICAL_EXTENDED"
    else:
        final_status = "REVIEW"

    # --------------------------------------------------------
    # Final metadata
    # --------------------------------------------------------

    final_rows.append({
        "kju_id": kju,
        "accession": accession,
        "organism": organism,
        "strain": d.get("strain", ""),
        "sequence_type": d.get(
            "sequence_type",
            c.get("sequence_type", "")
        ),
        "length": d.get(
            "length",
            c.get("length", "")
        ),
        "isolation_source": d.get("isolation_source", ""),
        "geo_loc_name": d.get("geo_loc_name", ""),
        "collection_date": d.get("collection_date", ""),
        "collected_by": d.get("collected_by", ""),
        "genbank_file": d.get("genbank_file", ""),
        "source_class": source_class,
        "final_status": final_status
    })

    audit_rows.append({
        "kju_id": kju,
        "accession": accession,
        "canonical_original": canonical_evidence,
        "direct_genbank_evidence": direct_genbank_evidence,
        "orphan_ncbi_evidence": orphan_ncbi_evidence,
        "accession_conflict": accession_conflict,
        "organism_conflict": organism_conflict,
        "final_status": final_status
    })

# ------------------------------------------------------------
# Sort
# ------------------------------------------------------------

final_df = pd.DataFrame(final_rows)
audit_df = pd.DataFrame(audit_rows)

final_df = final_df.sort_values(
    "kju_id",
    key=lambda s: s.str.extract(r"(\d+)")[0].astype(int)
).reset_index(drop=True)

audit_df = audit_df.sort_values(
    "kju_id",
    key=lambda s: s.str.extract(r"(\d+)")[0].astype(int)
).reset_index(drop=True)

# ------------------------------------------------------------
# Hard validation
# ------------------------------------------------------------

print()
print("=" * 100)
print(" FINAL VALIDATION")
print("=" * 100)

print("Rows:", len(final_df))
print("Unique KJU IDs:", final_df["kju_id"].nunique())
print("Unique accessions:", final_df["accession"].nunique())

# KJU duplicates
dup_kju = final_df[
    final_df["kju_id"].duplicated(keep=False)
]

if not dup_kju.empty:
    fatal_errors.append(
        "Duplicate KJU IDs detected"
    )

# accession duplicates
dup_acc = final_df[
    final_df["accession"].ne("")
    & final_df["accession"].duplicated(keep=False)
]

if not dup_acc.empty:
    fatal_errors.append(
        "Duplicate accession values detected"
    )

# Missing critical fields
missing_kju = final_df[
    final_df["kju_id"].eq("")
]

missing_acc = final_df[
    final_df["accession"].eq("")
]

missing_org = final_df[
    final_df["organism"].eq("")
]

if not missing_kju.empty:
    fatal_errors.append("Missing KJU IDs detected")

if not missing_acc.empty:
    fatal_errors.append("Missing accession detected")

if not missing_org.empty:
    fatal_errors.append("Missing organism detected")

# Expected 202
if len(final_df) != 202:
    fatal_errors.append(
        f"Expected 202 records but generated {len(final_df)}"
    )

# Expected statuses
review_count = (
    ~final_df["final_status"]
    .isin(["CANONICAL", "CANONICAL_EXTENDED"])
).sum()

if review_count > 0:
    fatal_errors.append(
        f"{review_count} records remain unresolved"
    )

# ------------------------------------------------------------
# Write only if completely clean
# ------------------------------------------------------------

if fatal_errors:

    print()
    print("!!! VALIDATION FAILED !!!")

    for e in fatal_errors:
        print("ERROR:", e)

    print()
    print("NO FINAL DATASET WRITTEN.")
    print("NO MASTER DATA MODIFIED.")

    sys.exit(1)

# ------------------------------------------------------------
# Write final dataset
# ------------------------------------------------------------

final_df.to_csv(
    OUT,
    index=False
)

audit_df.to_csv(
    AUDIT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

original_count = (
    final_df["source_class"]
    .eq("ORIGINAL_CANONICAL")
    .sum()
)

extended_count = (
    final_df["source_class"]
    .eq("NCBI_ORPHAN_CONFIRMED_BY_DIRECT_GENBANK")
    .sum()
)

with open(SUMMARY, "w") as f:

    f.write(
        "DODDAGGUBBI - FINAL KJU DATASET VALIDATION\n"
    )
    f.write("=" * 100 + "\n\n")

    f.write(
        f"Final records: {len(final_df)}\n"
    )

    f.write(
        f"Original canonical: {original_count}\n"
    )

    f.write(
        f"Extended NCBI-backed records: {extended_count}\n"
    )

    f.write(
        f"Unique KJU IDs: {final_df['kju_id'].nunique()}\n"
    )

    f.write(
        f"Unique accessions: {final_df['accession'].nunique()}\n"
    )

    f.write("\nSTATUS\n")
    f.write("-" * 100 + "\n")

    for status, n in final_df[
        "final_status"
    ].value_counts().items():

        f.write(f"{status}: {n}\n")

    f.write("\nVALIDATION\n")
    f.write("-" * 100 + "\n")
    f.write("Duplicate KJU IDs: 0\n")
    f.write("Duplicate accessions: 0\n")
    f.write("Missing KJU IDs: 0\n")
    f.write("Missing accessions: 0\n")
    f.write("Missing organisms: 0\n")
    f.write("Accession conflicts: 0\n")
    f.write("Organism conflicts: 0\n")
    f.write("Unresolved records: 0\n")

    f.write("\n")
    f.write(
        "MASTER_KJU_RECONCILIATION_001_183.csv was NOT modified.\n"
    )

print()
print("=" * 100)
print(" FINAL DATASET CREATED SUCCESSFULLY")
print("=" * 100)

print()
print("Final records :", len(final_df))
print("Original      :", original_count)
print("Extended      :", extended_count)

print()
print("Output:")
print(OUT)

print()
print("Audit:")
print(AUDIT)

print()
print("Summary:")
print(SUMMARY)

print()
print("VALIDATION: PASSED")
print("NO DUPLICATE KJU IDs")
print("NO DUPLICATE ACCESSIONS")
print("NO ACCESSION CONFLICTS")
print("NO ORGANISM CONFLICTS")
print("NO MISSING CRITICAL FIELDS")
print("NO UNRESOLVED RECORDS")

print()
print("Original master was NOT modified.")
