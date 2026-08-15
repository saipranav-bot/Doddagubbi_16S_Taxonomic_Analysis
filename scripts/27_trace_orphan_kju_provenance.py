from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

ORPHAN = META / "ORPHAN_NCBI_RECONCILIATION_V1.csv"

OUT = META / "ORPHAN_KJU_PROVENANCE_TRACE_V1.csv"
SUMMARY = META / "ORPHAN_KJU_PROVENANCE_TRACE_SUMMARY_V1.txt"

print("=" * 100)
print(" DODDAGGUBBI - ORPHAN KJU PROVENANCE TRACE")
print("=" * 100)

orphan = pd.read_csv(
    ORPHAN,
    dtype=str
).fillna("")

candidates = orphan[
    orphan["reconciliation_status"] ==
    "CANDIDATE_ADDITIONAL_KJU"
].copy()

kju_ids = sorted(
    candidates["kju_id"].dropna().unique()
)

print(f"\nCandidate KJUs: {len(kju_ids)}")

# ------------------------------------------------------------
# Search all metadata CSV/TXT files
# ------------------------------------------------------------

metadata_files = []

for f in META.rglob("*"):

    if not f.is_file():
        continue

    if f.name == OUT.name:
        continue

    if f.suffix.lower() in [".csv", ".txt"]:
        metadata_files.append(f)

rows = []

for kju in kju_ids:

    print(f"Tracing {kju} ...")

    hits = []

    for f in metadata_files:

        try:
            text = f.read_text(
                encoding="utf-8",
                errors="ignore"
            )
        except Exception:
            continue

        if kju in text:

            hits.append(f.name)

    candidate_rows = candidates[
        candidates["kju_id"] == kju
    ]

    if len(candidate_rows) == 0:
        continue

    r = candidate_rows.iloc[0]

    rows.append({
        "kju_id": kju,
        "ncbi_accession": r["accession"],
        "ncbi_header": r["ncbi_header"],
        "metadata_files_containing_kju":
            ";".join(sorted(set(hits))),
        "metadata_file_count": len(set(hits)),
    })

out = pd.DataFrame(rows)

out = out.sort_values("kju_id")

out.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

with open(SUMMARY, "w", encoding="utf-8") as fh:

    fh.write(
        "DODDAGGUBBI - ORPHAN KJU PROVENANCE TRACE\n"
    )
    fh.write("=" * 100 + "\n\n")

    fh.write(
        f"Candidate KJUs examined: {len(kju_ids)}\n\n"
    )

    fh.write("PROVENANCE RESULTS\n")
    fh.write("-" * 100 + "\n")

    for _, r in out.iterrows():

        fh.write(
            f"{r['kju_id']} | "
            f"{r['ncbi_accession']} | "
            f"metadata files: "
            f"{r['metadata_file_count']}\n"
        )

        if r["metadata_files_containing_kju"]:

            fh.write(
                "  "
                + r["metadata_files_containing_kju"]
                + "\n"
            )

    fh.write("\n")
    fh.write("IMPORTANT\n")
    fh.write("-" * 100 + "\n")
    fh.write(
        "This is a provenance investigation only.\n"
        "The canonical master was NOT modified.\n"
        "No KJU records were added.\n"
    )

print("\n" + "=" * 100)
print(" PROVENANCE TRACE COMPLETE")
print("=" * 100)

print(f"Output : {OUT}")
print(f"Summary: {SUMMARY}")

print("\nMASTER DATA WAS NOT MODIFIED.")
