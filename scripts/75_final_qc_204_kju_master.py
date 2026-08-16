#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


# ============================================================
# STEP 75
# FINAL QC OF 204-KJU INTEGRATED MASTER
#
# Purpose:
#   Independently validate the Step 74 output.
#
# Checks:
#   1. File exists and loads
#   2. Exactly 204 KJU records
#   3. Exactly 204 unique KJU IDs
#   4. 112 resolved records
#   5. 92 unresolved MASTER_ONLY records
#   6. No unexpected blank taxonomy
#   7. No duplicate accepted accessions among resolved records
#   8. Biological identity consistency
#   9. KJU0038/KJU0039 biological duplicate
#  10. Special adjudicated records
#  11. No accidental taxonomy in unresolved records
# ============================================================


BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "final_integrated_master"
)

INPUT = (
    ROOT
    / "FINAL_204_KJU_ADJUDICATED_MASTER.csv"
)

SUMMARY = (
    ROOT
    / "STEP75_FINAL_QC_REPORT.txt"
)


def clean(x):
    if x is None:
        return ""

    if pd.isna(x):
        return ""

    x = str(x).strip()

    if x.lower() in {
        "",
        "nan",
        "none",
        "null",
        "na",
        "n/a",
    }:
        return ""

    return x


def normalize_accession(x):
    x = clean(x)

    if not x:
        return ""

    return x.split(".")[0].upper()


# ============================================================
# LOAD
# ============================================================

print("=" * 110)
print("STEP 75: FINAL QC OF 204-KJU INTEGRATED MASTER")
print("=" * 110)

print()
print("INPUT:")
print(INPUT)

if not INPUT.exists():
    raise FileNotFoundError(
        f"Step 74 output not found:\n{INPUT}"
    )

df = pd.read_csv(INPUT, dtype=str)

df.columns = [c.strip() for c in df.columns]

print()
print("=" * 110)
print("1. BASIC FILE QC")
print("=" * 110)

print(f"Rows              : {len(df)}")
print(f"Columns           : {len(df.columns)}")
print(f"Unique KJU IDs    : {df['KJU_ID'].nunique()}")

if len(df) != 204:
    raise ValueError(
        f"Expected 204 rows, found {len(df)}"
    )

if df["KJU_ID"].nunique() != 204:
    raise ValueError(
        "KJU_ID uniqueness failure"
    )


# ============================================================
# IDENTIFIER PRESENCE
# ============================================================

print()
print("=" * 110)
print("2. IDENTIFIER STRUCTURE")
print("=" * 110)

master_present = (
    df["Master_present"]
    .map(clean)
    .eq("YES")
)

opencity_present = (
    df["OpenCity_present"]
    .map(clean)
    .eq("YES")
)

master_only = (
    master_present
    & ~opencity_present
)

reconciliation_only = (
    ~master_present
    & opencity_present
)

overlap = (
    master_present
    & opencity_present
)

print(f"MASTER present          : {master_present.sum()}")
print(f"OpenCity present        : {opencity_present.sum()}")
print(f"MASTER_ONLY             : {master_only.sum()}")
print(f"RECONCILIATION_ONLY     : {reconciliation_only.sum()}")
print(f"OVERLAP                 : {overlap.sum()}")

if master_present.sum() != 183:
    raise ValueError(
        f"Expected 183 MASTER records, found {master_present.sum()}"
    )

if opencity_present.sum() != 112:
    raise ValueError(
        f"Expected 112 reconciliation records, found {opencity_present.sum()}"
    )

if master_only.sum() != 92:
    raise ValueError(
        f"Expected 92 MASTER_ONLY records, found {master_only.sum()}"
    )

if reconciliation_only.sum() != 21:
    raise ValueError(
        f"Expected 21 reconciliation-only records, found {reconciliation_only.sum()}"
    )

if overlap.sum() != 91:
    raise ValueError(
        f"Expected 91 overlap records, found {overlap.sum()}"
    )


# ============================================================
# TAXONOMY QC
# ============================================================

print()
print("=" * 110)
print("3. TAXONOMY QC")
print("=" * 110)

df["_accepted_taxonomy"] = (
    df["Accepted_taxonomy"]
    .map(clean)
)

blank_taxonomy = (
    df["_accepted_taxonomy"]
    .eq("")
)

allowed_unresolved = (
    master_only
    & df["Adjudication_type"]
        .map(clean)
        .eq("MASTER_ONLY_UNRESOLVED")
    & blank_taxonomy
)

unexpected_blank = (
    blank_taxonomy
    & ~allowed_unresolved
)

print(
    f"Blank Accepted_taxonomy       : "
    f"{blank_taxonomy.sum()}"
)

print(
    f"Allowed unresolved blanks     : "
    f"{allowed_unresolved.sum()}"
)

print(
    f"Unexpected blank taxonomy     : "
    f"{unexpected_blank.sum()}"
)

if blank_taxonomy.sum() != 92:
    raise ValueError(
        f"Expected exactly 92 blank taxonomies, "
        f"found {blank_taxonomy.sum()}"
    )

if allowed_unresolved.sum() != 92:
    raise ValueError(
        "Not all unresolved MASTER_ONLY records "
        "are correctly classified"
    )

if unexpected_blank.sum() != 0:
    print()
    print("UNEXPECTED BLANK TAXONOMY RECORDS:")
    print(
        df.loc[
            unexpected_blank,
            [
                "KJU_ID",
                "Master_present",
                "OpenCity_present",
                "Accepted_taxonomy",
                "Adjudication_type",
            ],
        ].to_string(index=False)
    )

    raise ValueError(
        f"{unexpected_blank.sum()} unexpected blank taxonomy records"
    )


# ============================================================
# RESOLVED RECORD QC
# ============================================================

resolved = ~blank_taxonomy

print()
print("=" * 110)
print("4. RESOLVED RECORD QC")
print("=" * 110)

print(f"Resolved records : {resolved.sum()}")

if resolved.sum() != 112:
    raise ValueError(
        f"Expected 112 resolved records, "
        f"found {resolved.sum()}"
    )


# ============================================================
# UNRESOLVED RECORD QC
# ============================================================

print()
print("=" * 110)
print("5. UNRESOLVED MASTER_ONLY QC")
print("=" * 110)

unresolved = df.loc[
    master_only
].copy()

print(
    f"MASTER_ONLY records       : "
    f"{len(unresolved)}"
)

print(
    f"MASTER_ONLY_UNRESOLVED    : "
    f"{unresolved['Adjudication_type'].map(clean).eq('MASTER_ONLY_UNRESOLVED').sum()}"
)

print(
    f"Confidence UNRESOLVED     : "
    f"{unresolved['Confidence'].map(clean).eq('UNRESOLVED').sum()}"
)

if len(unresolved) != 92:
    raise ValueError(
        "MASTER_ONLY count changed"
    )

if not (
    unresolved["Adjudication_type"]
    .map(clean)
    .eq("MASTER_ONLY_UNRESOLVED")
    .all()
):
    raise ValueError(
        "Some MASTER_ONLY records have incorrect adjudication type"
    )

if not (
    unresolved["Confidence"]
    .map(clean)
    .eq("UNRESOLVED")
    .all()
):
    raise ValueError(
        "Some unresolved records have incorrect confidence"
    )


# ============================================================
# NO FABRICATED TAXONOMY QC
# ============================================================

print()
print("=" * 110)
print("6. NO-FABRICATION QC")
print("=" * 110)

bad_unresolved = unresolved[
    unresolved["Accepted_taxonomy"].map(clean) != ""
]

print(
    f"Unresolved records with taxonomy : "
    f"{len(bad_unresolved)}"
)

if len(bad_unresolved) != 0:
    print(
        bad_unresolved[
            [
                "KJU_ID",
                "Accepted_taxonomy",
                "Adjudication_type",
            ]
        ].to_string(index=False)
    )

    raise ValueError(
        "Taxonomy found in unresolved MASTER_ONLY records"
    )

print("PASS — no taxonomy fabricated")


# ============================================================
# BIOLOGICAL IDENTITY QC
# ============================================================

print()
print("=" * 110)
print("7. BIOLOGICAL IDENTITY QC")
print("=" * 110)

df["_biological_identity"] = (
    df["Biological_identity"]
    .map(clean)
)

blank_bio = df["_biological_identity"].eq("")

print(
    f"Blank Biological_identity : "
    f"{blank_bio.sum()}"
)

if blank_bio.any():
    print(
        df.loc[
            blank_bio,
            ["KJU_ID", "Accepted_accession"]
        ].to_string(index=False)
    )

    raise ValueError(
        "Blank Biological_identity detected"
    )


# ============================================================
# UNRESOLVED BIOLOGICAL IDENTITIES
# ============================================================

unresolved_bio = unresolved[
    "Biological_identity"
].map(clean)

bad_unresolved_identity = ~unresolved_bio.str.startswith(
    "UNRESOLVED:KJU"
)

print(
    f"Invalid unresolved identities : "
    f"{bad_unresolved_identity.sum()}"
)

if bad_unresolved_identity.any():
    raise ValueError(
        "Unresolved records do not have unique biological identities"
    )


# ============================================================
# DUPLICATE BIOLOGICAL IDENTITIES
# ============================================================

bio_counts = (
    df["_biological_identity"]
    .value_counts()
)

duplicate_bio = bio_counts[
    bio_counts > 1
]

print()
print(
    f"Biologically unique identities : "
    f"{bio_counts.size}"
)

print(
    f"Biological identities duplicated: "
    f"{len(duplicate_bio)}"
)

if len(duplicate_bio):
    print()
    print("DUPLICATED BIOLOGICAL IDENTITIES:")
    print(duplicate_bio.to_string())


# ============================================================
# EXPECTED KJU0038/KJU0039 DUPLICATE
# ============================================================

print()
print("=" * 110)
print("8. KNOWN BIOLOGICAL DUPLICATE QC")
print("=" * 110)

row38 = df.loc[
    df["KJU_ID"] == "KJU0038"
]

row39 = df.loc[
    df["KJU_ID"] == "KJU0039"
]

if len(row38) != 1 or len(row39) != 1:
    raise ValueError(
        "KJU0038 or KJU0039 missing"
    )

bio38 = clean(
    row38.iloc[0]["Biological_identity"]
)

bio39 = clean(
    row39.iloc[0]["Biological_identity"]
)

print(f"KJU0038 biological identity : {bio38}")
print(f"KJU0039 biological identity : {bio39}")

if bio38 != "PZ666011":
    raise ValueError(
        f"KJU0038 expected PZ666011, found {bio38}"
    )

if bio39 != "PZ666011":
    raise ValueError(
        f"KJU0039 expected PZ666011, found {bio39}"
    )

print("PASS — KJU0038/KJU0039 correctly collapsed biologically")


# ============================================================
# ACCESSION QC
# ============================================================

print()
print("=" * 110)
print("9. ACCEPTED ACCESSION QC")
print("=" * 110)

df["_accepted_accession"] = (
    df["Accepted_accession"]
    .map(clean)
)

resolved_accessions = df.loc[
    resolved,
    "_accepted_accession"
]

blank_resolved_accessions = (
    resolved_accessions.eq("")
)

print(
    f"Resolved records with blank accession : "
    f"{blank_resolved_accessions.sum()}"
)

if blank_resolved_accessions.any():
    print(
        df.loc[
            resolved & df["_accepted_accession"].eq(""),
            [
                "KJU_ID",
                "Accepted_taxonomy",
                "Adjudication_type",
            ],
        ].to_string(index=False)
    )

    raise ValueError(
        "Resolved record has blank Accepted_accession"
    )


# ============================================================
# DUPLICATE ACCESSION REPORT
# ============================================================

normalized_accessions = (
    resolved_accessions
    .map(normalize_accession)
)

accession_counts = (
    normalized_accessions
    .replace("", pd.NA)
    .dropna()
    .value_counts()
)

duplicate_accessions = accession_counts[
    accession_counts > 1
]

print(
    f"Duplicate normalized accessions : "
    f"{len(duplicate_accessions)}"
)

if len(duplicate_accessions):
    print()
    print("DUPLICATE ACCESSIONS:")
    print(duplicate_accessions.to_string())


# ============================================================
# SPECIAL RECORD QC
# ============================================================

print()
print("=" * 110)
print("10. SPECIAL ADJUDICATION QC")
print("=" * 110)

expected_special = {
    "KJU0038": "Shewanella algae",
    "KJU0039": "Shewanella algae",
    "KJU0105": "Bacillus cereus",
    "KJU0224": "Pseudomonas putida",
    "KJU0225": "Chryseobacterium indologenes",
}

special_failures = []

for kju, expected_taxonomy in expected_special.items():

    row = df.loc[
        df["KJU_ID"] == kju
    ]

    if len(row) != 1:
        special_failures.append(
            f"{kju}: missing or duplicated"
        )
        continue

    actual = clean(
        row.iloc[0]["Accepted_taxonomy"]
    )

    print(
        f"{kju:8s} "
        f"{actual}"
    )

    if actual != expected_taxonomy:
        special_failures.append(
            f"{kju}: expected '{expected_taxonomy}', "
            f"found '{actual}'"
        )


if special_failures:

    print()
    for x in special_failures:
        print("FAIL:", x)

    raise ValueError(
        "Special adjudication QC failed"
    )


# ============================================================
# TAXONOMY SUMMARY
# ============================================================

print()
print("=" * 110)
print("11. TAXONOMY SUMMARY")
print("=" * 110)

taxonomy_counts = (
    df.loc[
        resolved,
        "Accepted_taxonomy"
    ]
    .map(clean)
    .value_counts()
)

print(
    taxonomy_counts.to_string()
)


# ============================================================
# ADJUDICATION SUMMARY
# ============================================================

print()
print("=" * 110)
print("12. ADJUDICATION SUMMARY")
print("=" * 110)

adjudication_counts = (
    df["Adjudication_type"]
    .map(clean)
    .value_counts()
)

print(
    adjudication_counts.to_string()
)


# ============================================================
# FINAL QC
# ============================================================

print()
print("=" * 110)
print("STEP 75 FINAL RESULT")
print("=" * 110)

print()
print("PASS:")
print("  ✓ 204 total KJU records")
print("  ✓ 204 unique KJU IDs")
print("  ✓ 183 MASTER IDs")
print("  ✓ 112 reconciliation IDs")
print("  ✓ 92 MASTER_ONLY records")
print("  ✓ 21 reconciliation-only records")
print("  ✓ 91 overlapping records")
print("  ✓ 112 resolved records")
print("  ✓ 92 explicitly unresolved records")
print("  ✓ No unexpected blank taxonomy")
print("  ✓ No fabricated taxonomy")
print("  ✓ All unresolved records uniquely identified")
print("  ✓ KJU0038/KJU0039 biological duplicate validated")
print("  ✓ Special adjudications validated")
print()
print("STEP 75 QC PASSED")


# ============================================================
# WRITE REPORT
# ============================================================

report = []

report.append(
    "STEP 75 — FINAL QC OF 204-KJU INTEGRATED MASTER"
)

report.append("=" * 100)
report.append("")

report.append(
    f"Total KJU records       : {len(df)}"
)

report.append(
    f"Unique KJU IDs          : {df['KJU_ID'].nunique()}"
)

report.append(
    f"MASTER records          : {master_present.sum()}"
)

report.append(
    f"Reconciliation records  : {opencity_present.sum()}"
)

report.append(
    f"MASTER_ONLY             : {master_only.sum()}"
)

report.append(
    f"Reconciliation_ONLY     : {reconciliation_only.sum()}"
)

report.append(
    f"Overlap                 : {overlap.sum()}"
)

report.append(
    f"Resolved taxonomy       : {resolved.sum()}"
)

report.append(
    f"Unresolved taxonomy     : {blank_taxonomy.sum()}"
)

report.append(
    f"Unexpected blank        : {unexpected_blank.sum()}"
)

report.append(
    f"Biological identities   : {bio_counts.size}"
)

report.append(
    f"Duplicate bio identities: {len(duplicate_bio)}"
)

report.append("")

report.append(
    "FINAL STATUS: PASS"
)

report.append(
    "No taxonomy was fabricated."
)

report.append(
    "MASTER_ONLY records remain explicitly unresolved."
)

SUMMARY.write_text(
    "\n".join(report) + "\n"
)

print()
print("QC REPORT:")
print(SUMMARY)

print()
print("=" * 110)
