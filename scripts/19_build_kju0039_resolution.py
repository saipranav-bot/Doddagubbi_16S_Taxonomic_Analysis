from pathlib import Path
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

OUT = META / "KJU0039_FINAL_RESOLUTION_V1.csv"
SUMMARY = META / "KJU0039_FINAL_RESOLUTION_SUMMARY_V1.txt"

print("=" * 90)
print(" DODDAGGUBBI - KJU0039 FINAL RESOLUTION")
print("=" * 90)

# ------------------------------------------------------------
# Evidence
# ------------------------------------------------------------

row = {
    "kju_id": "KJU0039",

    "opencity_accession": "PZ666011",
    "opencity_organism": "Shewanella algae",

    "ncbi_accession_claimed": "PZ666011.1",
    "ncbi_accession_owner": "KJU0038",
    "ncbi_organism": "Shewanella algae",

    "ncbi_kju0039_search_hits": 0,

    "accession_conflict": True,
    "independent_ncbi_record": False,

    "resolution_status": "UNRESOLVED_ACCESSION_COLLISION",

    "resolved_accession": "",
    "resolved_organism": "",

    "action": "DO_NOT_MERGE",

    "recommended_master_action":
        "RETAIN KJU0039 AS UNRESOLVED; REMOVE PZ666011 "
        "FROM KJU0039 AS A VALIDATED ACCESSION ONLY AFTER "
        "SOURCE-LEVEL REVIEW",

    "evidence_basis":
        "NCBI contains PZ666011.1 as Shewanella algae "
        "strain KJU0038. NCBI searches for KJU0039 returned "
        "zero nucleotide records.",

    "confidence": "HIGH"
}

df = pd.DataFrame([row])

df.to_csv(OUT, index=False)

# ------------------------------------------------------------
# Human-readable summary
# ------------------------------------------------------------

summary = f"""
DODDAGGUBBI - KJU0039 FINAL RESOLUTION
================================================================================

KJU ID
--------------------------------------------------------------------------------
KJU0039

SOURCE EVIDENCE
--------------------------------------------------------------------------------
OpenCity accession : PZ666011
OpenCity organism  : Shewanella algae

NCBI EVIDENCE
--------------------------------------------------------------------------------
Accession found    : PZ666011.1
NCBI organism      : Shewanella algae
NCBI strain        : KJU0038
NCBI KJU0039 hits  : 0

CONFLICT
--------------------------------------------------------------------------------
PZ666011.1 is already assigned by NCBI to KJU0038.

No independent NCBI nucleotide record for KJU0039 was found.

RESOLUTION
--------------------------------------------------------------------------------
Status             : UNRESOLVED_ACCESSION_COLLISION
Resolved accession : NONE
Resolved organism  : NONE
Action             : DO_NOT_MERGE
Confidence         : HIGH

RECOMMENDATION
--------------------------------------------------------------------------------
KJU0039 must remain a separate unresolved sample identifier.

Do NOT assign PZ666011.1 to KJU0039.

Do NOT copy the KJU0038 sequence into KJU0039.

Do NOT delete KJU0039 from the project.

KJU0039 should remain in the master metadata as an unresolved
record until independent source evidence identifies its accession
or sequence.

AUDIT INTERPRETATION
--------------------------------------------------------------------------------
This is an accession identity collision, not a confirmed organism
identity for KJU0039.

The organism name "Shewanella algae" may remain as a source-reported
value, but it must NOT be treated as NCBI-validated for KJU0039.

MASTER DATA
--------------------------------------------------------------------------------
No master data was modified by this script.
"""

SUMMARY.write_text(summary.strip() + "\n", encoding="utf-8")

print("\nRESOLUTION CREATED")
print(f"Table   : {OUT}")
print(f"Summary : {SUMMARY}")
print("\nMASTER DATA WAS NOT MODIFIED.")

