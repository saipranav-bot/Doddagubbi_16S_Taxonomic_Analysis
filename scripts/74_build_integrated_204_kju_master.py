#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


# ============================================================
# STEP 74
# BUILD INTEGRATED 204-KJU ADJUDICATED MASTER
#
# IMPORTANT:
#
# There are 204 unique KJU IDs:
#   183 MASTER IDs
#   112 reconciliation IDs
#   204 combined unique IDs
#
# Of these, 92 are MASTER_ONLY records.
#
# Diagnostic confirmed that these 92 MASTER_ONLY records have:
#   accession     = blank
#   organism      = blank
#   sequence_type = blank
#   length        = blank
#
# They are therefore genuinely unresolved.
#
# NO TAXONOMY IS FABRICATED.
#
# For these records:
#   Accepted_taxonomy      = ""
#   Adjudication_type      = "MASTER_ONLY_UNRESOLVED"
#   Confidence             = "UNRESOLVED"
#   Taxonomic_evidence_source = "NO_TAXONOMY_AVAILABLE"
#   Biological_identity    = "UNRESOLVED:KJUxxxx"
#
# Blank Accepted_taxonomy is therefore permitted ONLY for
# MASTER_ONLY_UNRESOLVED records.
# ============================================================


from pathlib import Path
import pandas as pd


# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
)

MASTER_FILE = (
    BASE
    / "07_Master_Metadata"
    / "MASTER_KJU_RECONCILIATION_001_183.csv"
)

RECON_FILE = (
    ROOT
    / "OPENCITY_111_VS_FINAL_109_COMPLETE_RECONCILIATION.csv"
)

ADJ_FILE = (
    ROOT
    / "final_adjudication"
    / "FINAL_KJU_TAXONOMIC_ADJUDICATION.csv"
)

OUTDIR = (
    ROOT
    / "final_integrated_master"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTFILE = (
    OUTDIR
    / "FINAL_204_KJU_ADJUDICATED_MASTER.csv"
)

BIOFILE = (
    OUTDIR
    / "FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER.csv"
)

SUMMARY = (
    OUTDIR
    / "FINAL_204_KJU_ADJUDICATED_MASTER_SUMMARY.txt"
)


# ============================================================
# HELPERS
# ============================================================

def clean(value):
    """
    Convert NaN/None/empty-like values to "".
    """
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    value = str(value).strip()

    if value.lower() in {
        "",
        "nan",
        "none",
        "null",
        "na",
        "n/a",
    }:
        return ""

    return value


def normalize_accession(value):
    """
    Normalize GenBank accession by removing version suffix.

    Example:
        AB123456.1 -> AB123456
    """
    value = clean(value)

    if not value:
        return ""

    return value.split(".")[0].upper()


def get_value(row, column):
    """
    Safely extract a value from a pandas Series.
    """
    if row is None:
        return ""

    if column not in row.index:
        return ""

    return clean(row[column])


def first_nonempty(*values):
    """
    Return first non-empty value.
    """
    for value in values:
        value = clean(value)
        if value:
            return value

    return ""


# ============================================================
# LOAD INPUT FILES
# ============================================================

print("=" * 110)
print("STEP 74: BUILD INTEGRATED 204-KJU ADJUDICATED MASTER")
print("=" * 110)

print()
print("INPUT FILES")
print("-" * 110)
print(f"MASTER : {MASTER_FILE}")
print(f"RECON  : {RECON_FILE}")
print(f"ADJ    : {ADJ_FILE}")
print()


master = pd.read_csv(
    MASTER_FILE,
    dtype=str,
    keep_default_na=False
)

recon = pd.read_csv(
    RECON_FILE,
    dtype=str,
    keep_default_na=False
)

adj = pd.read_csv(
    ADJ_FILE,
    dtype=str,
    keep_default_na=False
)


# ============================================================
# CLEAN COLUMN NAMES
# ============================================================

master.columns = [
    str(c).strip()
    for c in master.columns
]

recon.columns = [
    str(c).strip()
    for c in recon.columns
]

adj.columns = [
    str(c).strip()
    for c in adj.columns
]


# ============================================================
# BASIC INPUT VALIDATION
# ============================================================

required_master = {
    "kju_id",
    "status",
    "accession",
    "organism",
    "sequence_type",
    "length",
}

required_recon = {
    "_KJU",
    "OpenCity_GenBank Acc. No.",
    "OpenCity_Species",
    "Final_accession",
    "Final_organism",
    "Final_sequence_type",
    "Final_length",
    "RECONCILIATION_STATUS",
}

required_adj = {
    "KJU_ID",
    "Accepted_taxonomy",
    "Final_sequence_type",
    "Final_length",
    "Reconciliation_status",
    "Adjudication_type",
    "Evidence",
    "Confidence",
    "Adjudication_note",
}


missing_master = required_master - set(master.columns)
missing_recon = required_recon - set(recon.columns)
missing_adj = required_adj - set(adj.columns)


if missing_master:
    raise ValueError(
        f"MASTER FILE missing columns: {sorted(missing_master)}"
    )

if missing_recon:
    raise ValueError(
        f"RECONCILIATION FILE missing columns: {sorted(missing_recon)}"
    )

if missing_adj:
    raise ValueError(
        f"ADJUDICATION FILE missing columns: {sorted(missing_adj)}"
    )


# ============================================================
# NORMALIZE IDENTIFIERS
# ============================================================

master["kju_id"] = master["kju_id"].map(clean)
recon["_KJU"] = recon["_KJU"].map(clean)
adj["KJU_ID"] = adj["KJU_ID"].map(clean)


# Remove completely blank IDs
master = master[
    master["kju_id"].ne("")
].copy()

recon = recon[
    recon["_KJU"].ne("")
].copy()

adj = adj[
    adj["KJU_ID"].ne("")
].copy()


# ============================================================
# IDENTIFIER INVENTORY
# ============================================================

master_ids = set(
    master["kju_id"]
)

recon_ids = set(
    recon["_KJU"]
)

adj_ids = set(
    adj["KJU_ID"]
)

all_ids = sorted(
    master_ids | recon_ids
)


print("INPUT SIZES")
print("-" * 110)
print(f"Master rows              : {len(master)}")
print(f"Reconciliation rows      : {len(recon)}")
print(f"Adjudication rows        : {len(adj)}")
print()

print("IDENTIFIER COUNTS")
print("-" * 110)
print(f"Master unique KJU IDs    : {len(master_ids)}")
print(f"Reconciliation unique IDs: {len(recon_ids)}")
print(f"Adjudication unique IDs  : {len(adj_ids)}")
print(f"Combined unique IDs      : {len(all_ids)}")
print()


# ============================================================
# EXPECTED 204 RECORDS
# ============================================================

if len(master_ids) != 183:
    raise ValueError(
        f"Expected 183 MASTER KJU IDs, found {len(master_ids)}"
    )

if len(recon_ids) != 112:
    raise ValueError(
        f"Expected 112 reconciliation KJU IDs, found {len(recon_ids)}"
    )

if len(adj_ids) != 112:
    raise ValueError(
        f"Expected 112 adjudication KJU IDs, found {len(adj_ids)}"
    )

if len(all_ids) != 204:
    raise ValueError(
        f"Expected 204 combined KJU IDs, found {len(all_ids)}"
    )


# ============================================================
# MASTER-ONLY / RECONCILIATION IDs
# ============================================================

master_only_ids = master_ids - recon_ids
recon_only_ids = recon_ids - master_ids
overlap_ids = master_ids & recon_ids


print("RECONCILIATION STRUCTURE")
print("-" * 110)
print(f"MASTER_ONLY IDs          : {len(master_only_ids)}")
print(f"RECONCILIATION_ONLY IDs  : {len(recon_only_ids)}")
print(f"OVERLAP IDs              : {len(overlap_ids)}")
print()


if len(master_only_ids) != 92:
    raise ValueError(
        f"Expected 92 MASTER_ONLY IDs, found {len(master_only_ids)}"
    )


# ============================================================
# INDEX TABLES
# ============================================================

master_idx = (
    master
    .drop_duplicates(
        "kju_id",
        keep="first"
    )
    .set_index("kju_id")
)

recon_idx = (
    recon
    .drop_duplicates(
        "_KJU",
        keep="first"
    )
    .set_index("_KJU")
)

adj_idx = (
    adj
    .drop_duplicates(
        "KJU_ID",
        keep="first"
    )
    .set_index("KJU_ID")
)


# ============================================================
# BUILD INTEGRATED RECORDS
# ============================================================

records = []


for kju in all_ids:

    m = (
        master_idx.loc[kju]
        if kju in master_idx.index
        else None
    )

    r = (
        recon_idx.loc[kju]
        if kju in recon_idx.index
        else None
    )

    a = (
        adj_idx.loc[kju]
        if kju in adj_idx.index
        else None
    )


    # ========================================================
    # MASTER INFORMATION
    # ========================================================

    master_accession = get_value(
        m,
        "accession"
    )

    master_organism = get_value(
        m,
        "organism"
    )

    master_sequence_type = get_value(
        m,
        "sequence_type"
    )

    master_length = get_value(
        m,
        "length"
    )

    master_status = get_value(
        m,
        "status"
    )


    # ========================================================
    # OPENCITY / RECONCILIATION INFORMATION
    # ========================================================

    opencity_accession = get_value(
        r,
        "OpenCity_GenBank Acc. No."
    )

    opencity_species = get_value(
        r,
        "OpenCity_Species"
    )

    opencity_status = get_value(
        r,
        "RECONCILIATION_STATUS"
    )


    # ========================================================
    # FINAL RECONCILIATION INFORMATION
    # ========================================================

    final_accession = get_value(
        r,
        "Final_accession"
    )

    final_organism = get_value(
        r,
        "Final_organism"
    )

    final_sequence_type = get_value(
        r,
        "Final_sequence_type"
    )

    final_length = get_value(
        r,
        "Final_length"
    )


    # ========================================================
    # ADJUDICATION INFORMATION
    # ========================================================

    accepted_taxonomy = get_value(
        a,
        "Accepted_taxonomy"
    )

    adjudication_type = get_value(
        a,
        "Adjudication_type"
    )

    confidence = get_value(
        a,
        "Confidence"
    )

    adjudication_note = get_value(
        a,
        "Adjudication_note"
    )

    evidence = get_value(
        a,
        "Evidence"
    )


    # ========================================================
    # DETERMINE WHETHER THIS IS MASTER-ONLY
    # ========================================================

    is_master_only = (
        kju in master_only_ids
    )


    # ========================================================
    # MASTER-ONLY RESOLUTION
    #
    # Diagnostic proved these 92 records contain no metadata
    # from which taxonomy can legitimately be derived.
    # ========================================================

    if is_master_only:

        # ----------------------------------------------------
        # These must remain unresolved.
        # ----------------------------------------------------

        if (
            master_accession
            or master_organism
            or master_sequence_type
            or master_length
        ):

            # This is a safety check. If this ever happens,
            # investigate the source rather than silently
            # changing the classification.
            raise ValueError(
                f"{kju} is MASTER_ONLY but contains unexpected "
                f"metadata: accession={master_accession!r}, "
                f"organism={master_organism!r}, "
                f"sequence_type={master_sequence_type!r}, "
                f"length={master_length!r}"
            )


        accepted_taxonomy = ""

        accepted_accession = ""

        accepted_sequence_type = ""

        accepted_length = ""

        adjudication_type = (
            "MASTER_ONLY_UNRESOLVED"
        )

        confidence = "UNRESOLVED"

        adjudication_note = (
            "MASTER_ONLY record; no accession, organism, "
            "sequence type, or sequence length available in "
            "the current master metadata or reconciliation "
            "sources. Taxonomy intentionally not assigned."
        )

        evidence_source = (
            "NO_TAXONOMY_AVAILABLE"
        )

        biological_identity = (
            f"UNRESOLVED:{kju}"
        )

        normalized_acc = ""


    # ========================================================
    # RECONCILIATION / ADJUDICATED RECORDS
    # ========================================================

    else:

        # ----------------------------------------------------
        # ACCEPTED ACCESSION
        # ----------------------------------------------------

        accepted_accession = first_nonempty(
            final_accession,
            master_accession,
            opencity_accession
        )


        # ----------------------------------------------------
        # ACCEPTED SEQUENCE TYPE
        # ----------------------------------------------------

        accepted_sequence_type = first_nonempty(
            final_sequence_type,
            get_value(a, "Final_sequence_type"),
            master_sequence_type
        )


        # ----------------------------------------------------
        # ACCEPTED LENGTH
        # ----------------------------------------------------

        accepted_length = first_nonempty(
            final_length,
            get_value(a, "Final_length"),
            master_length
        )


        # ----------------------------------------------------
        # TAXONOMY
        #
        # Priority:
        #   1. Step-72 accepted taxonomy
        #   2. Final organism
        #   3. Master organism
        #   4. OpenCity species
        #
        # This is ONLY used for records participating in
        # reconciliation.
        # ----------------------------------------------------

        accepted_taxonomy = first_nonempty(
            accepted_taxonomy,
            final_organism,
            master_organism,
            opencity_species
        )


        # ----------------------------------------------------
        # IF TAXONOMY STILL BLANK
        # ----------------------------------------------------

        if not accepted_taxonomy:

            raise ValueError(
                f"{kju} is a reconciliation record but has "
                f"no accepted taxonomy. "
                f"OpenCity={opencity_species!r}; "
                f"Final={final_organism!r}; "
                f"Master={master_organism!r}; "
                f"Adjudication={get_value(a, 'Accepted_taxonomy')!r}"
            )


        # ----------------------------------------------------
        # EVIDENCE SOURCE
        # ----------------------------------------------------

        if adjudication_type == (
            "OPENCITY_ONLY_NCBI_VALIDATED"
        ):

            evidence_source = (
                "NCBI_BLAST_PHYLOGENY_VALIDATION"
            )

        elif adjudication_type == (
            "OPEN_CITY_TAXONOMIC_ERROR"
        ):

            evidence_source = (
                "FINAL_GENBANK_TAXONOMY"
            )

        elif adjudication_type == (
            "OPEN_CITY_IDENTIFIER_DISCREPANCY"
        ):

            evidence_source = (
                "GENBANK_IDENTIFIER_RESOLUTION"
            )

        elif adjudication_type == (
            "FINAL_ONLY_GENBANK_RECORD"
        ):

            evidence_source = (
                "GENBANK_RECORD"
            )

        elif adjudication_type == (
            "OPEN_CITY_TAXONOMY_SPELLING_ERROR"
        ):

            evidence_source = (
                "FINAL_METADATA_TAXONOMY"
            )

        elif adjudication_type == (
            "MATCHED_FINAL_TAXONOMY"
        ):

            evidence_source = (
                "FINAL_METADATA"
            )

        elif evidence:

            evidence_source = (
                "STEP72_ADJUDICATION"
            )

        else:

            evidence_source = (
                "RECONCILIATION_METADATA"
            )


        # ----------------------------------------------------
        # BIOLOGICAL IDENTITY
        # ----------------------------------------------------

        normalized_acc = normalize_accession(
            accepted_accession
        )

        if normalized_acc:

            biological_identity = normalized_acc

        else:

            # This should be uncommon for reconciliation
            # records, but avoid collapsing records with
            # missing accession into one group.
            biological_identity = (
                f"NO_ACCESSION:{kju}"
            )


        # ----------------------------------------------------
        # KNOWN BIOLOGICAL DUPLICATE
        # ----------------------------------------------------

        # KJU0038 and KJU0039 both resolve to PZ666011.
        if kju in {
            "KJU0038",
            "KJU0039"
        }:

            biological_identity = (
                "PZ666011"
            )


    # ========================================================
    # FINAL RECORD
    # ========================================================

    records.append({

        "KJU_ID": kju,

        "Master_present":
            "YES"
            if kju in master_ids
            else "NO",

        "OpenCity_present":
            "YES"
            if kju in recon_ids
            else "NO",

        "Master_status":
            master_status,

        "Master_accession":
            master_accession,

        "Master_taxonomy":
            master_organism,

        "OpenCity_accession":
            opencity_accession,

        "OpenCity_taxonomy":
            opencity_species,

        "Final_accession":
            final_accession,

        "Final_taxonomy":
            final_organism,

        "Final_sequence_type":
            final_sequence_type,

        "Final_length":
            final_length,

        "Accepted_accession":
            accepted_accession,

        "Accepted_taxonomy":
            accepted_taxonomy,

        "Accepted_sequence_type":
            accepted_sequence_type,

        "Accepted_length":
            accepted_length,

        "Adjudication_type":
            adjudication_type,

        "Confidence":
            confidence,

        "Adjudication_reason":
            adjudication_note,

        "Adjudication_evidence":
            evidence,

        "Taxonomic_evidence_source":
            evidence_source,

        "Normalized_accession":
            normalized_acc,

        "Biological_identity":
            biological_identity,

        "Reconciliation_status":
            opencity_status,
    })


# ============================================================
# CREATE DATAFRAME
# ============================================================

df = pd.DataFrame(
    records
)


# ============================================================
# SORT KJU IDs NUMERICALLY
# ============================================================

df = (
    df
    .sort_values(
        "KJU_ID",
        key=lambda s:
            s.str.extract(
                r"(\d+)$"
            )[0].astype(int)
    )
    .reset_index(drop=True)
)


# ============================================================
# FINAL IDENTIFIER QC
# ============================================================

print()
print("=" * 110)
print("FINAL IDENTIFIER QC")
print("=" * 110)

print(
    f"Rows                     : {len(df)}"
)

print(
    f"Unique KJU IDs            : "
    f"{df['KJU_ID'].nunique()}"
)

print(
    f"MASTER_ONLY records       : "
    f"{(df['Adjudication_type'] == 'MASTER_ONLY_UNRESOLVED').sum()}"
)

print(
    f"Resolved records          : "
    f"{(df['Adjudication_type'] != 'MASTER_ONLY_UNRESOLVED').sum()}"
)


if len(df) != 204:
    raise ValueError(
        f"Final integrated dataset has {len(df)} rows; "
        f"expected 204"
    )


if df["KJU_ID"].nunique() != 204:
    raise ValueError(
        "Duplicate KJU IDs detected"
    )


# ============================================================
# TAXONOMY QC
#
# BLANK TAXONOMY IS ALLOWED ONLY FOR THE 92 UNRESOLVED
# MASTER_ONLY RECORDS.
# ============================================================

blank_taxonomy = (
    df["Accepted_taxonomy"]
    .fillna("")
    .astype(str)
    .str.strip()
    .eq("")
)


allowed_blank = (
    df["Adjudication_type"]
    .eq("MASTER_ONLY_UNRESOLVED")
)


unexpected_blank = (
    blank_taxonomy
    & ~allowed_blank
)


print()
print("TAXONOMY QC")
print("-" * 110)

print(
    f"Total blank Accepted_taxonomy : "
    f"{blank_taxonomy.sum()}"
)

print(
    f"Allowed unresolved blanks     : "
    f"{(blank_taxonomy & allowed_blank).sum()}"
)

print(
    f"Unexpected blank taxonomy     : "
    f"{unexpected_blank.sum()}"
)


if unexpected_blank.any():

    print()
    print(
        "UNEXPECTED BLANK TAXONOMY RECORDS"
    )
    print("-" * 110)

    print(
        df.loc[
            unexpected_blank,
            [
                "KJU_ID",
                "Master_present",
                "OpenCity_present",
                "Master_taxonomy",
                "OpenCity_taxonomy",
                "Final_taxonomy",
                "Accepted_taxonomy",
                "Adjudication_type",
                "Confidence",
            ]
        ].to_string(index=False)
    )

    raise ValueError(
        f"{unexpected_blank.sum()} resolved records "
        f"have blank Accepted_taxonomy"
    )


# ============================================================
# VERIFY ALL 92 MASTER-ONLY ARE UNRESOLVED
# ============================================================

master_only_rows = df[
    df["Master_present"].eq("YES")
    & df["OpenCity_present"].eq("NO")
].copy()


print()
print("MASTER-ONLY QC")
print("-" * 110)

print(
    f"MASTER_ONLY rows          : "
    f"{len(master_only_rows)}"
)

print(
    f"MASTER_ONLY unresolved    : "
    f"{(
        master_only_rows['Adjudication_type']
        == 'MASTER_ONLY_UNRESOLVED'
    ).sum()}"
)


if len(master_only_rows) != 92:
    raise ValueError(
        f"Expected 92 MASTER_ONLY rows, "
        f"found {len(master_only_rows)}"
    )


if not (
    master_only_rows["Adjudication_type"]
    .eq("MASTER_ONLY_UNRESOLVED")
    .all()
):

    raise ValueError(
        "Not all MASTER_ONLY records are marked "
        "MASTER_ONLY_UNRESOLVED"
    )


# ============================================================
# BIOLOGICALLY UNIQUE DATASET
# ============================================================

bio = (
    df
    .sort_values(
        [
            "Biological_identity",
            "KJU_ID"
        ]
    )
    .drop_duplicates(
        "Biological_identity",
        keep="first"
    )
    .copy()
)


# ============================================================
# KJU ALIASES
# ============================================================

alias_map = (
    df
    .groupby(
        "Biological_identity"
    )["KJU_ID"]
    .apply(
        lambda x:
            ";".join(
                sorted(
                    x.tolist()
                )
            )
    )
    .to_dict()
)


bio["KJU_ID_ALIASES"] = (
    bio["Biological_identity"]
    .map(alias_map)
)


# ============================================================
# BIOLOGICAL RECORD COUNT
# ============================================================

count_map = (
    df[
        "Biological_identity"
    ]
    .value_counts()
    .to_dict()
)


bio["BIOLOGICAL_RECORD_COUNT"] = (
    bio["Biological_identity"]
    .map(count_map)
)


# ============================================================
# BIOLOGICAL QC
# ============================================================

if (
    bio["Biological_identity"]
    .isna()
    .any()
):

    raise ValueError(
        "Blank Biological_identity detected"
    )


if (
    bio["Biological_identity"]
    .astype(str)
    .str.strip()
    .eq("")
    .any()
):

    raise ValueError(
        "Empty Biological_identity detected"
    )


# ============================================================
# SAVE FINAL 204 MASTER
# ============================================================

df.to_csv(
    OUTFILE,
    index=False
)


# ============================================================
# SAVE BIOLOGICALLY UNIQUE MASTER
# ============================================================

bio.to_csv(
    BIOFILE,
    index=False
)


# ============================================================
# TAXONOMIC SUMMARY
# ============================================================

resolved_taxonomy = df[
    df["Accepted_taxonomy"]
    .fillna("")
    .astype(str)
    .str.strip()
    .ne("")
]


taxonomy_counts = (
    resolved_taxonomy[
        "Accepted_taxonomy"
    ]
    .value_counts()
)


# ============================================================
# ADJUDICATION SUMMARY
# ============================================================

adjudication_counts = (
    df[
        "Adjudication_type"
    ]
    .fillna("")
    .astype(str)
    .value_counts()
)


confidence_counts = (
    df[
        "Confidence"
    ]
    .fillna("")
    .astype(str)
    .value_counts()
)


# ============================================================
# WRITE SUMMARY
# ============================================================

summary_lines = []


summary_lines.append(
    "FINAL 204-KJU ADJUDICATED MASTER DATASET"
)

summary_lines.append(
    "=" * 110
)

summary_lines.append("")

summary_lines.append(
    f"Original master KJU IDs       : "
    f"{len(master_ids)}"
)

summary_lines.append(
    f"Reconciliation KJU IDs        : "
    f"{len(recon_ids)}"
)

summary_lines.append(
    f"Combined unique KJU IDs       : "
    f"{len(df)}"
)

summary_lines.append(
    f"MASTER_ONLY unresolved        : "
    f"{len(master_only_rows)}"
)

summary_lines.append(
    f"Resolved records              : "
    f"{len(resolved_taxonomy)}"
)

summary_lines.append(
    f"Biologically unique records   : "
    f"{len(bio)}"
)

summary_lines.append("")


summary_lines.append(
    "TAXONOMY RESOLUTION"
)

summary_lines.append(
    "-" * 110
)

summary_lines.append(
    f"Blank Accepted_taxonomy       : "
    f"{blank_taxonomy.sum()}"
)

summary_lines.append(
    f"Allowed unresolved blanks     : "
    f"{(blank_taxonomy & allowed_blank).sum()}"
)

summary_lines.append(
    f"Unexpected blank taxonomy     : "
    f"{unexpected_blank.sum()}"
)

summary_lines.append("")


summary_lines.append(
    "ADJUDICATION TYPES"
)

summary_lines.append(
    "-" * 110
)

for key, value in adjudication_counts.items():

    summary_lines.append(
        f"{key}: {value}"
    )

summary_lines.append("")


summary_lines.append(
    "CONFIDENCE"
)

summary_lines.append(
    "-" * 110
)

for key, value in confidence_counts.items():

    summary_lines.append(
        f"{key}: {value}"
    )

summary_lines.append("")


summary_lines.append(
    "BIOLOGICALLY UNIQUE DATASET"
)

summary_lines.append(
    "-" * 110
)

summary_lines.append(
    f"Total biological identities : "
    f"{len(bio)}"
)

summary_lines.append(
    f"Total KJU records           : "
    f"{len(df)}"
)

summary_lines.append("")


summary_lines.append(
    "KNOWN BIOLOGICAL DUPLICATE"
)

summary_lines.append(
    "-" * 110
)

summary_lines.append(
    "PZ666011 -> KJU0038;KJU0039"
)

summary_lines.append("")


summary_lines.append(
    "OUTPUT FILES"
)

summary_lines.append(
    "-" * 110
)

summary_lines.append(
    str(OUTFILE)
)

summary_lines.append(
    str(BIOFILE)
)

summary_lines.append(
    str(SUMMARY)
)


SUMMARY.write_text(
    "\n".join(summary_lines)
    + "\n"
)


# ============================================================
# FINAL CONSOLE REPORT
# ============================================================

print()
print("=" * 110)
print("STEP 74 COMPLETE")
print("=" * 110)

print()
print("FINAL DATASET")
print("-" * 110)

print(
    f"Total KJU records          : {len(df)}"
)

print(
    f"Unique KJU IDs             : "
    f"{df['KJU_ID'].nunique()}"
)

print(
    f"Resolved taxonomy records  : "
    f"{len(resolved_taxonomy)}"
)

print(
    f"Unresolved MASTER_ONLY     : "
    f"{len(master_only_rows)}"
)

print(
    f"Biologically unique        : "
    f"{len(bio)}"
)


print()
print("TAXONOMY QC")
print("-" * 110)

print(
    f"Blank taxonomy             : "
    f"{blank_taxonomy.sum()}"
)

print(
    f"Allowed unresolved blanks  : "
    f"{(blank_taxonomy & allowed_blank).sum()}"
)

print(
    f"Unexpected blanks          : "
    f"{unexpected_blank.sum()}"
)


print()
print("MASTER-ONLY RECORDS")
print("-" * 110)

for kju in sorted(
    master_only_ids,
    key=lambda x: int(x.replace("KJU", ""))
):

    row = df[
        df["KJU_ID"] == kju
    ]

    if len(row) == 1:

        rr = row.iloc[0]

        print(
            f"{kju:<10} "
            f"{rr['Adjudication_type']:<30} "
            f"{rr['Confidence']}"
        )


print()
print("SPECIAL RECORDS")
print("-" * 110)

special_ids = [
    "KJU0038",
    "KJU0039",
    "KJU0105",
    "KJU0224",
    "KJU0225",
]

for kju in special_ids:

    row = df[
        df["KJU_ID"] == kju
    ]

    if len(row):

        rr = row.iloc[0]

        print(
            f"{kju} -> "
            f"{rr['Accepted_taxonomy']} | "
            f"{rr['Adjudication_type']} | "
            f"{rr['Confidence']} | "
            f"accession={rr['Accepted_accession']}"
        )


print()
print("OUTPUT")
print("-" * 110)

print(
    f"FINAL 204 MASTER:"
)

print(
    OUTFILE
)

print()

print(
    f"BIOLOGICALLY UNIQUE MASTER:"
)

print(
    BIOFILE
)

print()

print(
    f"SUMMARY:"
)

print(
    SUMMARY
)

print()
print("=" * 110)
print("NO TAXONOMY WAS FABRICATED")
print("MASTER_ONLY records are retained as explicitly unresolved.")
print("=" * 110)

