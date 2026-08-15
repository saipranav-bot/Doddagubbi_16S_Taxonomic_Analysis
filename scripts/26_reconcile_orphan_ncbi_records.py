from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

CANONICAL_FILE = META / "CANONICAL_KJU_METADATA_V1.csv"
AUDIT_FILE = META / "DEEP_SEQUENCE_INVENTORY_AUDIT_V1.csv"

OUT = META / "ORPHAN_NCBI_RECONCILIATION_V1.csv"
SUMMARY = META / "ORPHAN_NCBI_RECONCILIATION_SUMMARY_V1.txt"

print("=" * 100)
print(" DODDAGGUBBI - ORPHAN NCBI RECORD RECONCILIATION")
print("=" * 100)

canonical = pd.read_csv(
    CANONICAL_FILE,
    dtype=str
).fillna("")

audit = pd.read_csv(
    AUDIT_FILE,
    dtype=str
).fillna("")

canonical_ids = set(canonical["kju_id"].str.strip())

rows = []

kju_pattern = r"\bKJU\d{4}\b"

for _, r in audit.iterrows():

    if r["classification"] != "ACCESSION_NOT_IN_CANONICAL":
        continue

    file_name = r["file_name"]
    accession = r["filename_accession"]
    header = r["header_text"]

    kju_hits = sorted(set(
        x.upper()
        for x in re.findall(
            kju_pattern,
            header,
            flags=re.IGNORECASE
        )
    ))

    if kju_hits:

        for kju in kju_hits:

            rows.append({
                "file_type": r["file_type"],
                "file_name": file_name,
                "accession": accession,
                "kju_id": kju,
                "ncbi_header": header,
                "kju_in_canonical": kju in canonical_ids,
                "organism_from_header": header,
                "reconciliation_status":
                    "CANDIDATE_ADDITIONAL_KJU"
                    if kju not in canonical_ids
                    else "ALREADY_CANONICAL"
            })

    else:

        rows.append({
            "file_type": r["file_type"],
            "file_name": file_name,
            "accession": accession,
            "kju_id": "",
            "ncbi_header": header,
            "kju_in_canonical": False,
            "organism_from_header": header,
            "reconciliation_status":
                "NO_KJU_ID_IN_NCBI_HEADER"
        })

out = pd.DataFrame(rows)

out = out.sort_values(
    ["reconciliation_status", "kju_id", "accession"]
)

out.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

candidate = out[
    out["reconciliation_status"] ==
    "CANDIDATE_ADDITIONAL_KJU"
]

no_kju = out[
    out["reconciliation_status"] ==
    "NO_KJU_ID_IN_NCBI_HEADER"
]

with open(SUMMARY, "w", encoding="utf-8") as fh:

    fh.write(
        "DODDAGGUBBI - ORPHAN NCBI RECORD RECONCILIATION\n"
    )
    fh.write("=" * 100 + "\n\n")

    fh.write(f"Orphan NCBI records examined: {len(out)}\n")
    fh.write(
        f"KJU-linked candidate records: {len(candidate)}\n"
    )
    fh.write(
        f"Records without KJU ID: {len(no_kju)}\n\n"
    )

    fh.write("CANDIDATE ADDITIONAL KJUs\n")
    fh.write("-" * 100 + "\n")

    for _, r in candidate.iterrows():

        fh.write(
            f"{r['kju_id']} -> "
            f"{r['accession']} -> "
            f"{r['ncbi_header']}\n"
        )

    fh.write("\n")
    fh.write("RECORDS WITHOUT KJU ID\n")
    fh.write("-" * 100 + "\n")

    for _, r in no_kju.iterrows():

        fh.write(
            f"{r['accession']} -> "
            f"{r['ncbi_header']}\n"
        )

    fh.write("\n")
    fh.write("IMPORTANT\n")
    fh.write("-" * 100 + "\n")
    fh.write(
        "These records are candidates only.\n"
        "The canonical master was NOT modified.\n"
        "No KJU IDs were added.\n"
        "No organism assignments were changed.\n"
    )

print("\n" + "=" * 100)
print(" RECONCILIATION COMPLETE")
print("=" * 100)

print(f"Output : {OUT}")
print(f"Summary: {SUMMARY}")

print("\nCandidate additional KJUs:",
      len(candidate))

print("Records without KJU ID:",
      len(no_kju))

print("\nMASTER DATA WAS NOT MODIFIED.")
