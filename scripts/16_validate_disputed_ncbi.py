from pathlib import Path
import pandas as pd
import requests
import time
import re

# ============================================================
# DODDAGGUBBI
# NCBI VALIDATION OF DISPUTED KJU RECORDS
# ============================================================

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

RECON_FILE = META / "OPENCITY_MASTER_RECONCILIATION_V2.csv"

OUT = META / "NCBI_VALIDATION_DISPUTED_RECORDS.csv"

NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

print("=" * 90)
print(" DODDAGGUBBI - NCBI VALIDATION OF DISPUTED RECORDS")
print("=" * 90)

# ------------------------------------------------------------
# 1. READ RECONCILIATION
# ------------------------------------------------------------

if not RECON_FILE.exists():
    raise SystemExit(
        f"ERROR: File not found:\n{RECON_FILE}"
    )

df = pd.read_csv(RECON_FILE, dtype=str)

# ------------------------------------------------------------
# 2. SELECT RECORDS REQUIRING VALIDATION
# ------------------------------------------------------------

# Explicitly important records
important_kjus = {
    "KJU0008",
    "KJU0071",
    "KJU0076",
    "KJU0084",
    "KJU0096",
    "KJU0105",
    "KJU0038",
    "KJU0039",
}

# All OpenCity-only records
opencity_only = set(
    df.loc[
        df["reconciliation_status"] == "OPENCITY_ONLY",
        "kju_id"
    ].dropna()
)

target_kjus = important_kjus.union(opencity_only)

target = df[
    df["kju_id"].isin(target_kjus)
].copy()

target = target.sort_values(
    "kju_id"
)

print(
    f"\nRecords selected for NCBI validation : "
    f"{len(target)}"
)

print(
    f"OpenCity-only records                 : "
    f"{len(opencity_only)}"
)

# ------------------------------------------------------------
# 3. GET ACCESSIONS
# ------------------------------------------------------------

def clean_accession(value):

    if pd.isna(value):
        return None

    value = str(value).strip()

    if not value:
        return None

    return value


def accession_with_version(value):

    value = clean_accession(value)

    if value is None:
        return None

    # NCBI EFetch works best with accession.version.
    # If source lacks version, first try accession as supplied.
    return value


# Collect unique accession values
accessions = {}

for _, row in target.iterrows():

    kju = row["kju_id"]

    candidates = [
        row.get("accession_original_master"),
        row.get("accession_original_opencity"),
    ]

    for accession in candidates:

        accession = accession_with_version(accession)

        if accession:

            if kju not in accessions:
                accessions[kju] = []

            if accession not in accessions[kju]:
                accessions[kju].append(accession)


# ------------------------------------------------------------
# 4. NCBI FETCH FUNCTION
# ------------------------------------------------------------

def fetch_ncbi(accession):

    accession = clean_accession(accession)

    if accession is None:
        return {
            "ncbi_status": "NO_ACCESSION",
            "ncbi_accession": None,
            "ncbi_definition": None,
            "ncbi_organism": None,
            "ncbi_taxid": None,
            "ncbi_length": None,
            "ncbi_moltype": None,
            "ncbi_raw": None,
        }

    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
    }

    try:

        response = requests.get(
            NCBI_URL,
            params=params,
            timeout=30,
        )

        if response.status_code != 200:

            return {
                "ncbi_status":
                    f"HTTP_{response.status_code}",
                "ncbi_accession": accession,
                "ncbi_definition": None,
                "ncbi_organism": None,
                "ncbi_taxid": None,
                "ncbi_length": None,
                "ncbi_moltype": None,
                "ncbi_raw": response.text[:500],
            }

        text = response.text

        if (
            not text.strip()
            or "Error" in text[:500]
            or "Nothing has been found" in text
        ):

            return {
                "ncbi_status": "NOT_FOUND",
                "ncbi_accession": accession,
                "ncbi_definition": None,
                "ncbi_organism": None,
                "ncbi_taxid": None,
                "ncbi_length": None,
                "ncbi_moltype": None,
                "ncbi_raw": text[:500],
            }

        # ----------------------------------------------------
        # Parse LOCUS
        # ----------------------------------------------------

        locus_match = re.search(
            r"^LOCUS\s+(\S+)\s+(\d+)\s+bp\s+(\S+)",
            text,
            re.MULTILINE,
        )

        ncbi_accession = None
        ncbi_length = None
        ncbi_moltype = None

        if locus_match:

            ncbi_accession = locus_match.group(1)
            ncbi_length = locus_match.group(2)
            ncbi_moltype = locus_match.group(3)

        # ----------------------------------------------------
        # Parse DEFINITION
        # ----------------------------------------------------

        definition_match = re.search(
            r"^DEFINITION\s+(.+?)(?=\nACCESSION|\nVERSION)",
            text,
            re.MULTILINE | re.DOTALL,
        )

        ncbi_definition = None

        if definition_match:

            ncbi_definition = re.sub(
                r"\s+",
                " ",
                definition_match.group(1)
            ).strip()

        # ----------------------------------------------------
        # Parse VERSION
        # ----------------------------------------------------

        version_match = re.search(
            r"^VERSION\s+(\S+)",
            text,
            re.MULTILINE,
        )

        if version_match:
            ncbi_accession = version_match.group(1)

        # ----------------------------------------------------
        # Parse ORGANISM
        # ----------------------------------------------------

        organism_match = re.search(
            r"^\s+ORGANISM\s+(.+)$",
            text,
            re.MULTILINE,
        )

        ncbi_organism = None

        if organism_match:
            ncbi_organism = (
                organism_match.group(1).strip()
            )

        # ----------------------------------------------------
        # Parse TaxID
        # ----------------------------------------------------

        taxid_match = re.search(
            r"/db_xref=\"taxon:(\d+)\"",
            text,
        )

        ncbi_taxid = None

        if taxid_match:
            ncbi_taxid = taxid_match.group(1)

        return {
            "ncbi_status": "FOUND",
            "ncbi_accession": ncbi_accession,
            "ncbi_definition": ncbi_definition,
            "ncbi_organism": ncbi_organism,
            "ncbi_taxid": ncbi_taxid,
            "ncbi_length": ncbi_length,
            "ncbi_moltype": ncbi_moltype,
            "ncbi_raw": None,
        }

    except Exception as e:

        return {
            "ncbi_status": "ERROR",
            "ncbi_accession": accession,
            "ncbi_definition": None,
            "ncbi_organism": None,
            "ncbi_taxid": None,
            "ncbi_length": None,
            "ncbi_moltype": None,
            "ncbi_raw": str(e),
        }


# ------------------------------------------------------------
# 5. VALIDATE
# ------------------------------------------------------------

results = []

print("\n")
print("=" * 90)
print(" NCBI ACCESSION VALIDATION")
print("=" * 90)

for kju in sorted(accessions):

    for accession in accessions[kju]:

        print(
            f"\nValidating {kju} -> {accession}"
        )

        result = fetch_ncbi(accession)

        result["kju_id"] = kju
        result["source_accession"] = accession

        results.append(result)

        print(
            f"Status   : {result['ncbi_status']}"
        )

        print(
            f"Accession: {result['ncbi_accession']}"
        )

        print(
            f"Organism : {result['ncbi_organism']}"
        )

        print(
            f"TaxID    : {result['ncbi_taxid']}"
        )

        time.sleep(0.35)


# ------------------------------------------------------------
# 6. WRITE NCBI RESULTS
# ------------------------------------------------------------

ncbi_df = pd.DataFrame(results)

# ------------------------------------------------------------
# 7. Attach source information
# ------------------------------------------------------------

source_cols = [
    "kju_id",
    "accession_original_master",
    "accession_original_opencity",
    "organism_original_master",
    "organism_original_opencity",
    "reconciliation_status",
]

source_info = target[
    [
        c for c in source_cols
        if c in target.columns
    ]
].copy()

ncbi_df = ncbi_df.merge(
    source_info,
    on="kju_id",
    how="left",
)

# ------------------------------------------------------------
# 8. NCBI COMPARISON
# ------------------------------------------------------------

def norm(value):

    if pd.isna(value):
        return None

    value = str(value).strip().lower()

    value = re.sub(
        r"\s+",
        " ",
        value
    )

    return value


ncbi_df["master_organism_normalized"] = (
    ncbi_df[
        "organism_original_master"
    ].apply(norm)
)

ncbi_df["opencity_organism_normalized"] = (
    ncbi_df[
        "organism_original_opencity"
    ].apply(norm)
)

ncbi_df["ncbi_organism_normalized"] = (
    ncbi_df[
        "ncbi_organism"
    ].apply(norm)
)

ncbi_df["master_matches_ncbi"] = (
    ncbi_df["master_organism_normalized"].notna()
    &
    ncbi_df["ncbi_organism_normalized"].notna()
    &
    (
        ncbi_df["master_organism_normalized"]
        ==
        ncbi_df["ncbi_organism_normalized"]
    )
)

ncbi_df["opencity_matches_ncbi"] = (
    ncbi_df["opencity_organism_normalized"].notna()
    &
    ncbi_df["ncbi_organism_normalized"].notna()
    &
    (
        ncbi_df["opencity_organism_normalized"]
        ==
        ncbi_df["ncbi_organism_normalized"]
    )
)

def determine_resolution(row):

    if row["ncbi_status"] != "FOUND":
        return "NCBI_NOT_VALIDATED"

    master_match = row["master_matches_ncbi"]
    opencity_match = row["opencity_matches_ncbi"]

    if master_match and opencity_match:
        return "BOTH_AGREE_WITH_NCBI"

    if master_match and not opencity_match:
        return "MASTER_SUPPORTED_BY_NCBI"

    if opencity_match and not master_match:
        return "OPENCITY_SUPPORTED_BY_NCBI"

    return "NEITHER_SOURCE_MATCHES_NCBI"


ncbi_df["validation_resolution"] = ncbi_df.apply(
    determine_resolution,
    axis=1
)

# ------------------------------------------------------------
# 9. Sort
# ------------------------------------------------------------

ncbi_df["kju_number"] = pd.to_numeric(
    ncbi_df["kju_id"].str.extract(
        r"KJU(\d+)",
        expand=False
    ),
    errors="coerce"
)

ncbi_df = ncbi_df.sort_values(
    [
        "kju_number",
        "source_accession",
    ]
)

# ------------------------------------------------------------
# 10. Save
# ------------------------------------------------------------

ncbi_df.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# 11. Print summary
# ------------------------------------------------------------

print("\n")
print("=" * 90)
print(" NCBI VALIDATION SUMMARY")
print("=" * 90)

print(
    f"\nTotal accession checks : {len(ncbi_df)}"
)

print(
    "\nNCBI status:"
)

print(
    ncbi_df[
        "ncbi_status"
    ].value_counts(dropna=False).to_string()
)

print(
    "\nValidation resolution:"
)

print(
    ncbi_df[
        "validation_resolution"
    ].value_counts(dropna=False).to_string()
)

# ------------------------------------------------------------
# 12. Print important conflicts
# ------------------------------------------------------------

print("\n")
print("=" * 90)
print(" RECORDS REQUIRING ATTENTION")
print("=" * 90)

attention = ncbi_df[
    ncbi_df["validation_resolution"].isin(
        [
            "MASTER_SUPPORTED_BY_NCBI",
            "OPENCITY_SUPPORTED_BY_NCBI",
            "NEITHER_SOURCE_MATCHES_NCBI",
            "NCBI_NOT_VALIDATED",
        ]
    )
]

if len(attention) == 0:

    print(
        "\nNo unresolved NCBI validation records."
    )

else:

    print(
        attention[
            [
                "kju_id",
                "source_accession",
                "organism_original_master",
                "organism_original_opencity",
                "ncbi_accession",
                "ncbi_organism",
                "ncbi_taxid",
                "validation_resolution",
            ]
        ].to_string(index=False)
    )

# ------------------------------------------------------------
# 13. Final
# ------------------------------------------------------------

print("\n")
print("=" * 90)
print(" NCBI VALIDATION COMPLETE")
print("=" * 90)

print(
    f"\nOutput:"
    f"\n{OUT}"
)

print("\nNo master data was modified.")

