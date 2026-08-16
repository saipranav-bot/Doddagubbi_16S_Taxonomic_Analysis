#!/usr/bin/env python3

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

RECON_FILE = (
    ROOT
    / "OPENCITY_111_VS_FINAL_109_COMPLETE_RECONCILIATION.csv"
)

BIO_FILE = (
    ROOT
    / "biological_reconciliation"
    / "SHARED_KJU_BIOLOGICAL_RECONCILIATION.csv"
)

TAX_FILE = (
    ROOT
    / "taxonomic_validation"
    / "FINAL_KJU0224_KJU0225_TAXONOMIC_VALIDATION.csv"
)

OUTDIR = ROOT / "final_adjudication"
OUTDIR.mkdir(parents=True, exist_ok=True)

FINAL_FILE = (
    OUTDIR
    / "FINAL_KJU_TAXONOMIC_ADJUDICATION.csv"
)

BIOLOGICAL_FILE = (
    OUTDIR
    / "FINAL_UNIQUE_BIOLOGICAL_RECORDS.csv"
)

SUMMARY_FILE = (
    OUTDIR
    / "FINAL_KJU_TAXONOMIC_ADJUDICATION_SUMMARY.txt"
)


# ============================================================
# MANUAL SCIENTIFIC ADJUDICATIONS
# ============================================================

# These are the 13 records where the OpenCity name is clearly
# misspelled while the Final/GenBank taxonomy gives the
# corrected scientific name.

NAME_CORRECTIONS = {

    "KJU0008": {
        "accepted_taxonomy": "Citrobacter braakii",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Citrobacter brakii; corrected to Citrobacter braakii."
    },

    "KJU0071": {
        "accepted_taxonomy": "Lwoffella lincolnii",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity genus spelling differs from GenBank-backed record."
    },

    "KJU0076": {
        "accepted_taxonomy": "Paenibacillus polymyxa",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Paenibaillus polymyxa; corrected to Paenibacillus polymyxa."
    },

    "KJU0084": {
        "accepted_taxonomy": "Shewanella putrefaciens",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity contains a spelling error in the species name."
    },

    "KJU0096": {
        "accepted_taxonomy": "Staphylococcus epidermidis",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Staphlyococcus epidermidis; corrected to Staphylococcus epidermidis."
    },

    "KJU0193": {
        "accepted_taxonomy": "Glutamicibacter nicotianae",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity contains a spelling error."
    },

    "KJU0195": {
        "accepted_taxonomy": "Microbacterium oxydans",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Microbacterium oxidans; corrected to Microbacterium oxydans."
    },

    "KJU0205": {
        "accepted_taxonomy": "Delftia acidovorans",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity contains a spelling error."
    },

    "KJU0210": {
        "accepted_taxonomy": "Shewanella putrefaciens",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Shwanella putrefecians; corrected to Shewanella putrefaciens."
    },

    "KJU0213": {
        "accepted_taxonomy": "Paenibacillus amylolyticus",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity contains a spelling error."
    },

    "KJU0222": {
        "accepted_taxonomy": "Pseudomonas aeruginosa",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Pseudomonas aerugonisa; corrected to Pseudomonas aeruginosa."
    },

    "KJU0232": {
        "accepted_taxonomy": "Shewanella putrefaciens",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Shwanella putrefaciens; corrected genus spelling."
    },

    "KJU0234": {
        "accepted_taxonomy": "Shewanella putrefaciens",
        "issue": "OPEN_CITY_TAXONOMY_SPELLING_ERROR",
        "evidence": "Final GenBank-backed metadata",
        "confidence": "HIGH",
        "note": "OpenCity: Shwanella putrefaciens; corrected genus spelling."
    },
}


# KJU0105 is a genuine species discrepancy and was resolved
# directly against its GenBank record PZ710590.1.

SPECIAL_ADJUDICATIONS = {

    "KJU0105": {
        "accepted_taxonomy": "Bacillus cereus",
        "issue": "OPEN_CITY_TAXONOMIC_ERROR",
        "evidence": "GenBank PZ710590.1",
        "confidence": "VERY_HIGH",
        "note": (
            "GenBank explicitly identifies PZ710590.1 as "
            "Bacillus cereus strain KJU0105."
        )
    },

    "KJU0224": {
        "accepted_taxonomy": "Pseudomonas putida",
        "issue": "OPENCITY_ONLY_NCBI_VALIDATED",
        "evidence": (
            "NCBI BLAST + NCBI-backed phylogenetic placement; "
            "Step 70 validation"
        ),
        "confidence": "VERY_HIGH",
        "note": (
            "99.21% identity, 100% coverage; "
            "10/10 genus and 8/10 species phylogenetic agreement."
        )
    },

    "KJU0225": {
        "accepted_taxonomy": "Chryseobacterium indologenes",
        "issue": "OPENCITY_ONLY_NCBI_VALIDATED",
        "evidence": (
            "NCBI BLAST + NCBI-backed phylogenetic placement; "
            "Step 70 validation"
        ),
        "confidence": "VERY_HIGH",
        "note": (
            "99.55% identity, 100% coverage; "
            "10/10 genus and 6/10 species phylogenetic agreement."
        )
    },

    "KJU0038": {
        "accepted_taxonomy": "Shewanella algae",
        "issue": "FINAL_ONLY_GENBANK_RECORD",
        "evidence": "GenBank PZ666011.1",
        "confidence": "VERY_HIGH",
        "note": (
            "GenBank PZ666011.1 explicitly identifies strain KJU0038 "
            "as Shewanella algae."
        )
    },

    "KJU0039": {
        "accepted_taxonomy": "Shewanella algae",
        "issue": "OPEN_CITY_IDENTIFIER_DISCREPANCY",
        "evidence": "GenBank PZ666011.1",
        "confidence": "VERY_HIGH",
        "note": (
            "OpenCity assigns KJU0039 to PZ666011, while the "
            "authoritative GenBank record PZ666011.1 identifies "
            "strain KJU0038. Treat as identifier discrepancy, "
            "not as a second biological sequence."
        )
    },
}


# ============================================================
# LOAD DATA
# ============================================================

print("=" * 100)
print("STEP 72: FINAL TAXONOMIC ADJUDICATION")
print("=" * 100)

recon = pd.read_csv(RECON_FILE)

bio = pd.read_csv(BIO_FILE)

print(f"Reconciliation rows: {len(recon)}")
print(f"Biological reconciliation rows: {len(bio)}")


# ============================================================
# PREPARE BIOLOGICAL RECONCILIATION LOOKUP
# ============================================================

bio_lookup = {}

for _, row in bio.iterrows():

    kju = str(row["_KJU"]).strip()

    bio_lookup[kju] = {
        "accession_status": row["ACCESSION_RECONCILIATION"],
        "species_status": row["SPECIES_RECONCILIATION"],
        "overall_status": row["OVERALL_BIOLOGICAL_STATUS"],
    }


# ============================================================
# PREPARE STEP 70 TAXONOMIC VALIDATION
# ============================================================

tax_lookup = {}

if TAX_FILE.exists():

    tax = pd.read_csv(TAX_FILE)

    for _, row in tax.iterrows():

        kju = str(row.iloc[0]).strip()

        tax_lookup[kju] = row.to_dict()


# ============================================================
# BUILD FINAL ADJUDICATION
# ============================================================

records = []


for _, row in recon.iterrows():

    kju = str(row["_KJU"]).strip()

    open_city_accession = (
        "" if pd.isna(row["OpenCity_GenBank Acc. No."])
        else str(row["OpenCity_GenBank Acc. No."]).strip()
    )

    final_accession = (
        "" if pd.isna(row["Final_accession"])
        else str(row["Final_accession"]).strip()
    )

    open_city_species = (
        "" if pd.isna(row["OpenCity_Species"])
        else str(row["OpenCity_Species"]).strip()
    )

    final_organism = (
        "" if pd.isna(row["Final_organism"])
        else str(row["Final_organism"]).strip()
    )

    merge_status = str(row["_merge"]).strip()

    # --------------------------------------------------------
    # DEFAULT ADJUDICATION
    # --------------------------------------------------------

    accepted_taxonomy = ""
    issue_type = ""
    evidence_source = ""
    confidence = ""
    note = ""

    # --------------------------------------------------------
    # SPECIAL RECORDS FIRST
    # --------------------------------------------------------

    if kju in SPECIAL_ADJUDICATIONS:

        info = SPECIAL_ADJUDICATIONS[kju]

        accepted_taxonomy = info["accepted_taxonomy"]
        issue_type = info["issue"]
        evidence_source = info["evidence"]
        confidence = info["confidence"]
        note = info["note"]

    elif kju in NAME_CORRECTIONS:

        info = NAME_CORRECTIONS[kju]

        accepted_taxonomy = info["accepted_taxonomy"]
        issue_type = info["issue"]
        evidence_source = info["evidence"]
        confidence = info["confidence"]
        note = info["note"]

    # --------------------------------------------------------
    # NORMAL MATCHED RECORDS
    # --------------------------------------------------------

    elif merge_status == "both":

        if final_organism:

            accepted_taxonomy = final_organism

            issue_type = "MATCHED_FINAL_TAXONOMY"

            evidence_source = (
                "Final GenBank-backed metadata"
            )

            confidence = "HIGH"

            note = (
                "OpenCity and Final metadata share the same KJU ID "
                "and the accession/species comparison supports the "
                "final taxonomy."
            )

        elif open_city_species:

            accepted_taxonomy = open_city_species

            issue_type = "MATCHED_OPENCITY_TAXONOMY"

            evidence_source = "OpenCity metadata"

            confidence = "MODERATE"

            note = (
                "Final organism field unavailable; retained "
                "OpenCity taxonomy."
            )

        else:

            issue_type = "MISSING_TAXONOMY"
            evidence_source = "No taxonomic field available"
            confidence = "UNRESOLVED"
            note = "No usable taxonomy found."

    # --------------------------------------------------------
    # FINAL-ONLY
    # --------------------------------------------------------

    elif merge_status == "right_only":

        if final_organism:

            accepted_taxonomy = final_organism

            issue_type = "FINAL_ONLY"

            evidence_source = (
                "Final GenBank-backed metadata"
            )

            confidence = "HIGH"

            note = (
                "KJU exists in Final metadata but not in original "
                "OpenCity inventory."
            )

        else:

            issue_type = "FINAL_ONLY_MISSING_TAXONOMY"
            evidence_source = "Final metadata"
            confidence = "UNRESOLVED"
            note = "Final-only record without taxonomy."

    # --------------------------------------------------------
    # OPENCITY-ONLY
    # --------------------------------------------------------

    elif merge_status == "left_only":

        if open_city_species:

            accepted_taxonomy = open_city_species

            issue_type = "OPENCITY_ONLY"

            evidence_source = "OpenCity metadata"

            confidence = "MODERATE"

            note = (
                "OpenCity-only record; no corresponding Final "
                "metadata record."
            )

        else:

            issue_type = "OPENCITY_ONLY_MISSING_TAXONOMY"
            evidence_source = "OpenCity metadata"
            confidence = "UNRESOLVED"
            note = "OpenCity-only record without taxonomy."

    # --------------------------------------------------------
    # UNKNOWN
    # --------------------------------------------------------

    else:

        issue_type = "UNRESOLVED"
        evidence_source = "Unknown"
        confidence = "UNRESOLVED"
        note = "Unexpected reconciliation state."

    # --------------------------------------------------------
    # BIOLOGICAL RECONCILIATION INFORMATION
    # --------------------------------------------------------

    if kju in bio_lookup:

        accession_status = bio_lookup[kju]["accession_status"]
        species_status = bio_lookup[kju]["species_status"]
        overall_bio = bio_lookup[kju]["overall_status"]

    else:

        accession_status = "NOT_COMPARABLE"
        species_status = "NOT_COMPARABLE"
        overall_bio = "NOT_COMPARABLE"

    # --------------------------------------------------------
    # APPEND
    # --------------------------------------------------------

    records.append({

        "KJU_ID": kju,

        "OpenCity_accession": open_city_accession,

        "Final_accession": final_accession,

        "OpenCity_taxonomy": open_city_species,

        "Final_taxonomy": final_organism,

        "Accepted_taxonomy": accepted_taxonomy,

        "Final_sequence_type": row["Final_sequence_type"],

        "Final_length": row["Final_length"],

        "Reconciliation_status": (
            "MATCHED"
            if merge_status == "both"
            else (
                "OPENCITY_ONLY"
                if merge_status == "left_only"
                else (
                    "FINAL_ONLY"
                    if merge_status == "right_only"
                    else merge_status
                )
            )
        ),

        "Accession_reconciliation": accession_status,

        "Species_reconciliation": species_status,

        "Biological_reconciliation": overall_bio,

        "Adjudication_type": issue_type,

        "Evidence": evidence_source,

        "Confidence": confidence,

        "Adjudication_note": note,
    })


final_df = pd.DataFrame(records)


# ============================================================
# SORT
# ============================================================

final_df = final_df.sort_values(
    by="KJU_ID"
).reset_index(drop=True)


# ============================================================
# SAVE FINAL ADJUDICATION
# ============================================================

final_df.to_csv(
    FINAL_FILE,
    index=False
)


# ============================================================
# UNIQUE BIOLOGICAL RECORD VIEW
# ============================================================
#
# KJU0039 is an identifier discrepancy pointing to the
# KJU0038 GenBank record. Therefore it must not be counted
# as a second biological sequence.
#
# Keep KJU0039 in the full audit table but exclude it from
# the unique biological record count.

unique_df = final_df[
    final_df["KJU_ID"] != "KJU0039"
].copy()

unique_df["Biological_record_note"] = ""

unique_df.loc[
    unique_df["KJU_ID"] == "KJU0038",
    "Biological_record_note"
] = (
    "Authoritative GenBank record PZ666011.1; "
    "KJU0039 is an OpenCity identifier discrepancy."
)

unique_df.to_csv(
    BIOLOGICAL_FILE,
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

with open(SUMMARY_FILE, "w") as f:

    f.write("=" * 100 + "\n")
    f.write("DODDAGUBBI FINAL KJU TAXONOMIC ADJUDICATION\n")
    f.write("=" * 100 + "\n\n")

    f.write(f"Total KJU rows in reconciliation: {len(final_df)}\n")
    f.write(
        f"Unique biological records after KJU0039 identifier "
        f"deduplication: {len(unique_df)}\n\n"
    )

    # --------------------------------------------------------
    # RECONCILIATION
    # --------------------------------------------------------

    f.write("RECONCILIATION STATUS\n")
    f.write("-" * 100 + "\n")

    for status, count in (
        final_df["Reconciliation_status"]
        .value_counts()
        .items()
    ):
        f.write(f"{status}: {count}\n")

    f.write("\n")

    # --------------------------------------------------------
    # BIOLOGICAL RECONCILIATION
    # --------------------------------------------------------

    f.write("BIOLOGICAL RECONCILIATION\n")
    f.write("-" * 100 + "\n")

    for status, count in (
        final_df["Biological_reconciliation"]
        .value_counts()
        .items()
    ):
        f.write(f"{status}: {count}\n")

    f.write("\n")

    # --------------------------------------------------------
    # ADJUDICATION TYPES
    # --------------------------------------------------------

    f.write("ADJUDICATION TYPES\n")
    f.write("-" * 100 + "\n")

    for status, count in (
        final_df["Adjudication_type"]
        .value_counts()
        .items()
    ):
        f.write(f"{status}: {count}\n")

    f.write("\n")

    # --------------------------------------------------------
    # CONFIDENCE
    # --------------------------------------------------------

    f.write("CONFIDENCE\n")
    f.write("-" * 100 + "\n")

    for status, count in (
        final_df["Confidence"]
        .value_counts()
        .items()
    ):
        f.write(f"{status}: {count}\n")

    f.write("\n")

    # --------------------------------------------------------
    # EXPLICIT SPECIAL CASES
    # --------------------------------------------------------

    f.write("SPECIAL ADJUDICATIONS\n")
    f.write("-" * 100 + "\n")

    for kju in [
        "KJU0038",
        "KJU0039",
        "KJU0105",
        "KJU0224",
        "KJU0225",
    ]:

        row = final_df[
            final_df["KJU_ID"] == kju
        ]

        if len(row) == 1:

            r = row.iloc[0]

            f.write(
                f"{kju} -> {r['Accepted_taxonomy']} | "
                f"{r['Adjudication_type']} | "
                f"{r['Confidence']}\n"
            )

    f.write("\n")

    # --------------------------------------------------------
    # TAXONOMY CORRECTIONS
    # --------------------------------------------------------

    f.write("TAXONOMY NAME CORRECTIONS\n")
    f.write("-" * 100 + "\n")

    for kju in sorted(NAME_CORRECTIONS):

        row = final_df[
            final_df["KJU_ID"] == kju
        ]

        if len(row) == 1:

            r = row.iloc[0]

            f.write(
                f"{kju}: "
                f"{r['OpenCity_taxonomy']} "
                f"-> "
                f"{r['Accepted_taxonomy']}\n"
            )

    f.write("\n")

    f.write("=" * 100 + "\n")
    f.write("STEP 72 COMPLETE\n")
    f.write("=" * 100 + "\n")


# ============================================================
# CONSOLE OUTPUT
# ============================================================

print("\nFINAL ADJUDICATION CREATED")
print("-" * 100)

print(f"Final adjudication CSV:")
print(FINAL_FILE)

print(f"\nUnique biological records:")
print(BIOLOGICAL_FILE)

print(f"\nSummary:")
print(SUMMARY_FILE)

print("\nADJUDICATION TYPES:")
print(
    final_df["Adjudication_type"]
    .value_counts()
    .to_string()
)

print("\nCONFIDENCE:")
print(
    final_df["Confidence"]
    .value_counts()
    .to_string()
)

print("\nSPECIAL RECORDS:")

for kju in [
    "KJU0038",
    "KJU0039",
    "KJU0105",
    "KJU0224",
    "KJU0225",
]:

    row = final_df[
        final_df["KJU_ID"] == kju
    ]

    if len(row) == 1:

        r = row.iloc[0]

        print(
            f"{kju} -> "
            f"{r['Accepted_taxonomy']} | "
            f"{r['Adjudication_type']} | "
            f"{r['Confidence']}"
        )

print("\n" + "=" * 100)
print("STEP 72 COMPLETE")
print("=" * 100)

