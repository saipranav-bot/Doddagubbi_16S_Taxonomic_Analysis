from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

PROVENANCE = META / "ORPHAN_KJU_PROVENANCE_TRACE_V1.csv"

OUT = META / "ORPHAN_KJU_SCOPE_CLASSIFICATION_V1.csv"
SUMMARY = META / "ORPHAN_KJU_SCOPE_CLASSIFICATION_SUMMARY_V1.txt"

print("=" * 100)
print(" DODDAGGUBBI - ORPHAN KJU SCOPE CLASSIFICATION")
print("=" * 100)

df = pd.read_csv(PROVENANCE, dtype=str).fillna("")

# ------------------------------------------------------------
# Candidate KJUs
# ------------------------------------------------------------

candidates = df["kju_id"].dropna().unique().tolist()

# ------------------------------------------------------------
# Search project metadata for independent evidence
# ------------------------------------------------------------

search_files = []

for root in [META, MASTER / "00_Input", MASTER / "01_NCBI",
             MASTER / "02_Sequence_Data", MASTER / "03_Analysis"]:

    if not root.exists():
        continue

    for f in root.rglob("*"):
        if f.is_file():
            search_files.append(f)

print("Files searched:", len(search_files))

rows = []

for kju in candidates:

    accession = df.loc[
        df["kju_id"] == kju, "ncbi_accession"
    ].iloc[0]

    matching_files = []
    matching_lines = []

    patterns = [
        kju,
        accession,
        accession.replace(".1", "")
    ]

    for f in search_files:

        try:
            text = f.read_text(
                encoding="utf-8",
                errors="ignore"
            )
        except Exception:
            continue

        found = False

        for pattern in patterns:
            if pattern and pattern in text:
                found = True
                break

        if found:
            matching_files.append(str(f.relative_to(MASTER)))

            # Capture a few informative lines
            lines = text.splitlines()

            for line in lines:
                if any(
                    p in line
                    for p in patterns
                    if p
                ):
                    clean = line.strip()

                    if clean:
                        matching_lines.append(clean[:500])

                    if len(matching_lines) >= 10:
                        break

    rows.append({
        "kju_id": kju,
        "ncbi_accession": accession,
        "files_with_evidence": len(set(matching_files)),
        "evidence_files": ";".join(sorted(set(matching_files))),
        "example_evidence": " || ".join(
            matching_lines[:10]
        ),
        "provisional_classification": "NEEDS_REVIEW"
    })

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

result = pd.DataFrame(rows)

result.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

with open(SUMMARY, "w", encoding="utf-8") as fh:

    fh.write(
        "DODDAGGUBBI - ORPHAN KJU SCOPE CLASSIFICATION\n"
    )
    fh.write("=" * 100 + "\n\n")

    fh.write(
        f"Candidate KJUs examined: {len(result)}\n\n"
    )

    fh.write(
        "IMPORTANT\n"
        "This analysis does NOT modify the canonical master.\n"
        "All records remain provisional until dataset scope is established.\n\n"
    )

    fh.write(
        "NEXT DECISION\n"
        "Each orphan KJU must be classified as one of:\n"
        "EXPECTED_MEMBER\n"
        "EXCLUDED_BY_DESIGN\n"
        "DUPLICATE\n"
        "LATER_NCBI_RECORD\n"
        "UNRELATED_RECORD\n"
        "NEEDS_REVIEW\n"
    )

print()
print("=" * 100)
print("SCOPE CLASSIFICATION COMPLETE")
print("=" * 100)

print("Output :", OUT)
print("Summary:", SUMMARY)
print()
print("Candidate KJUs:", len(result))
print()
print(result[
    [
        "kju_id",
        "ncbi_accession",
        "files_with_evidence",
        "provisional_classification"
    ]
].to_string(index=False))

print()
print("MASTER DATA WAS NOT MODIFIED.")

