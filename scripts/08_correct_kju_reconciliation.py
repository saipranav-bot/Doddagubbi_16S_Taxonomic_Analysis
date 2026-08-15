from pathlib import Path
import csv
import re
import pandas as pd
from Bio import SeqIO

MASTER = Path.home() / "Doddagubbi_Master"
SOURCE = Path.home() / "Doddagubbi_16S"

GB_DIR = SOURCE / "genbank"
OUT = MASTER / "07_Master_Metadata"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print(" DODDAGGUBBI - CORRECT KJU / GENBANK RECONCILIATION")
print("=" * 80)

# ------------------------------------------------------------
# 1. READ ALL GENBANK FILES DIRECTLY
# ------------------------------------------------------------

print("\n[1/6] Reading GenBank files directly...")
print("-" * 80)

records = []

gb_files = sorted(GB_DIR.glob("*.gb"))

if not gb_files:
    raise SystemExit(f"ERROR: No GenBank files found in {GB_DIR}")

for gb in gb_files:

    try:
        rec = SeqIO.read(gb, "genbank")

        accession = rec.id

        organism = rec.annotations.get("organism", "")

        qualifiers = {}
        for feature in rec.features:
            for key, values in feature.qualifiers.items():
                if values:
                    qualifiers.setdefault(key, []).extend(values)

        # Extract strain / KJU ID
        strain_values = (
            qualifiers.get("strain", [])
            + qualifiers.get("isolate", [])
            + qualifiers.get("culture_collection", [])
        )

        strain_text = "; ".join(strain_values)

        kju_match = re.search(r"\bKJU\d{4}\b", strain_text, re.I)

        kju_id = kju_match.group(0).upper() if kju_match else ""

        # Search the whole record if strain qualifier did not work
        if not kju_id:
            text = str(rec)

            kju_match = re.search(r"\bKJU\d{4}\b", text, re.I)

            if kju_match:
                kju_id = kju_match.group(0).upper()

        # Detect sequence type
        text_upper = str(rec).upper()

        if "16S RIBOSOMAL RNA" in text_upper or "16S RRNA" in text_upper:
            seq_type = "16S"

        elif "INTERNAL TRANSCRIBED SPACER" in text_upper or "ITS" in text_upper:
            seq_type = "ITS"

        else:
            seq_type = "OTHER"

        source = ""
        geo = ""
        collection_date = ""
        collected_by = ""

        source_values = qualifiers.get("isolation_source", [])
        geo_values = qualifiers.get("geo_loc_name", [])
        date_values = qualifiers.get("collection_date", [])
        collector_values = qualifiers.get("collected_by", [])

        if source_values:
            source = source_values[0]

        if geo_values:
            geo = geo_values[0]

        if date_values:
            collection_date = date_values[0]

        if collector_values:
            collected_by = collector_values[0]

        records.append({
            "kju_id": kju_id,
            "accession": accession,
            "organism": organism,
            "strain": strain_text,
            "sequence_type": seq_type,
            "length": len(rec.seq),
            "isolation_source": source,
            "geo_loc_name": geo,
            "collection_date": collection_date,
            "collected_by": collected_by,
            "genbank_file": str(gb)
        })

    except Exception as e:
        print(f"ERROR reading {gb.name}: {e}")

df = pd.DataFrame(records)

print(f"GenBank files found     : {len(gb_files)}")
print(f"Records successfully read: {len(df)}")

# ------------------------------------------------------------
# 2. SAVE DIRECT GENBANK METADATA
# ------------------------------------------------------------

direct_file = OUT / "direct_genbank_kju_metadata.csv"

df.to_csv(direct_file, index=False)

print(f"\nSaved:")
print(direct_file)

# ------------------------------------------------------------
# 3. SHOW KJU EXTRACTION
# ------------------------------------------------------------

print("\n[2/6] KJU IDs extracted from GenBank")
print("-" * 80)

with_kju = df[df["kju_id"] != ""].copy()

print(f"Records containing KJU IDs: {len(with_kju)}")
print(f"Unique KJU IDs            : {with_kju['kju_id'].nunique()}")

if len(with_kju):
    print(
        with_kju[
            ["kju_id", "accession", "organism", "sequence_type", "length"]
        ]
        .sort_values("kju_id")
        .to_string(index=False)
    )

# ------------------------------------------------------------
# 4. READ TARGETED NCBI RESULTS
# ------------------------------------------------------------

print("\n[3/6] Reading targeted NCBI search results")
print("-" * 80)

hits_file = OUT / "KJU_NCBI_targeted_hits.csv"

if hits_file.exists():

    hits = pd.read_csv(hits_file)

    print(f"NCBI hit rows: {len(hits)}")

    # Display columns so we know exactly what the previous script produced
    print("\nNCBI result columns:")
    print(", ".join(hits.columns))

else:

    print("WARNING: targeted NCBI results not found.")

    hits = pd.DataFrame(columns=["kju_id"])

# ------------------------------------------------------------
# 5. BUILD EXPECTED KJU LIST
# ------------------------------------------------------------

print("\n[4/6] Building complete KJU reconciliation")
print("-" * 80)

expected = []

# Bacterial KJU0001-KJU0183
for i in range(1, 184):
    expected.append({
        "kju_id": f"KJU{i:04d}",
        "expected_type": "Bacterial"
    })

# Fungal KJU identifiers are kept as a separate expected range
# based on the project inventory. They are NOT assumed to have
# the same numbering unless confirmed by source data.
#
# Therefore, we only reconcile the explicitly known KJU0001-KJU0183
# series here.

expected_df = pd.DataFrame(expected)

# Existing KJU IDs
existing_map = (
    with_kju
    .drop_duplicates("kju_id")
    .set_index("kju_id")
)

# NCBI hit KJU IDs
if "kju_id" in hits.columns:

    hits["kju_id"] = (
        hits["kju_id"]
        .astype(str)
        .str.upper()
        .str.extract(r"(KJU\d{4})", expand=False)
    )

    hit_ids = set(hits["kju_id"].dropna())

else:

    hit_ids = set()

rows = []

for _, row in expected_df.iterrows():

    kju = row["kju_id"]

    in_existing = kju in existing_map.index
    in_ncbi = kju in hit_ids

    if in_existing:
        status = "DOWNLOADED"

    elif in_ncbi:
        status = "NCBI_FOUND_NOT_DOWNLOADED"

    else:
        status = "NOT_FOUND"

    accession = ""
    organism = ""
    sequence_type = ""
    length = ""

    if in_existing:

        x = existing_map.loc[kju]

        accession = x["accession"]
        organism = x["organism"]
        sequence_type = x["sequence_type"]
        length = x["length"]

    rows.append({
        "kju_id": kju,
        "expected_type": row["expected_type"],
        "status": status,
        "accession": accession,
        "organism": organism,
        "sequence_type": sequence_type,
        "length": length,
        "ncbi_hit": "YES" if in_ncbi else "NO"
    })

recon = pd.DataFrame(rows)

recon_file = OUT / "KJU_CORRECTED_RECONCILIATION.csv"

recon.to_csv(recon_file, index=False)

# ------------------------------------------------------------
# 6. FINAL SUMMARY
# ------------------------------------------------------------

print("\n[5/6] FINAL RECONCILIATION SUMMARY")
print("=" * 80)

print(f"GenBank files                         : {len(gb_files)}")
print(f"GenBank records                       : {len(df)}")

print(f"\nRecords with extracted KJU ID         : {len(with_kju)}")
print(f"Unique extracted KJU IDs              : {with_kju['kju_id'].nunique()}")

print(f"\nExpected KJU IDs                      : {len(expected_df)}")
print(f"Already downloaded                    : {(recon.status == 'DOWNLOADED').sum()}")
print(f"NCBI hit but not downloaded           : {(recon.status == 'NCBI_FOUND_NOT_DOWNLOADED').sum()}")
print(f"No NCBI hit                           : {(recon.status == 'NOT_FOUND').sum()}")

print("\nStatus breakdown:")
print(recon["status"].value_counts().to_string())

print("\nSequence type among downloaded KJU records:")
print(
    recon[
        recon["status"] == "DOWNLOADED"
    ]["sequence_type"]
    .value_counts()
    .to_string()
)

# ------------------------------------------------------------
# SAVE SUBSETS
# ------------------------------------------------------------

downloaded = recon[recon["status"] == "DOWNLOADED"]

ncbi_missing = recon[
    recon["status"] == "NCBI_FOUND_NOT_DOWNLOADED"
]

not_found = recon[
    recon["status"] == "NOT_FOUND"
]

downloaded.to_csv(
    OUT / "KJU_DOWNLOADED.csv",
    index=False
)

ncbi_missing.to_csv(
    OUT / "KJU_NCBI_FOUND_NOT_DOWNLOADED.csv",
    index=False
)

not_found.to_csv(
    OUT / "KJU_NOT_FOUND.csv",
    index=False
)

# ------------------------------------------------------------
# PRINT IMPORTANT TABLES
# ------------------------------------------------------------

print("\n[6/6] IMPORTANT RESULTS")
print("=" * 80)

print("\nKJU IDs already present locally:")
print(
    downloaded[
        ["kju_id", "accession", "organism", "sequence_type", "length"]
    ].to_string(index=False)
)

print("\nKJU IDs found by NCBI but not present locally:")
print(
    ncbi_missing[
        ["kju_id", "ncbi_hit"]
    ].to_string(index=False)
)

print("\nKJU IDs with no targeted NCBI hit:")
print(
    not_found[
        ["kju_id"]
    ].to_string(index=False)
)

print("\n" + "=" * 80)
print(" RECONCILIATION COMPLETE")
print("=" * 80)

print("\nOUTPUT FILES")
print("-" * 80)

print(direct_file)
print(recon_file)
print(OUT / "KJU_DOWNLOADED.csv")
print(OUT / "KJU_NCBI_FOUND_NOT_DOWNLOADED.csv")
print(OUT / "KJU_NOT_FOUND.csv")

