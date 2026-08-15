from pathlib import Path
import pandas as pd
import re
import hashlib

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

UNRESOLVED = META / "KJU_UNRESOLVED_METADATA_V1.csv"

OUT = META / "KJU_DEEP_RESOLUTION_V1.csv"
SUMMARY = META / "KJU_DEEP_RESOLUTION_SUMMARY_V1.txt"

print("=" * 110)
print(" DODDAGGUBBI - DEEP RESOLUTION OF UNRESOLVED KJUs")
print("=" * 110)

unresolved = pd.read_csv(
    UNRESOLVED,
    dtype=str
).fillna("")

candidate_kjus = unresolved["kju_id"].tolist()

print("Unresolved KJUs:", len(candidate_kjus))


# ------------------------------------------------------------
# Search roots
# ------------------------------------------------------------

roots = [
    MASTER / "00_Input",
    MASTER / "01_NCBI",
    MASTER / "02_Sequence_Data",
    MASTER / "03_Analysis",
    META,
]

files = []

for root in roots:

    if not root.exists():
        continue

    for f in root.rglob("*"):

        if f.is_file():

            # Ignore generated summary/audit files where possible.
            if f.name.startswith("KJU_DEEP_RESOLUTION"):
                continue

            files.append(f)

print("Files searched:", len(files))


# ------------------------------------------------------------
# Patterns
# ------------------------------------------------------------

accession_pattern = re.compile(
    r"\b(?:[A-Z]{1,4}\d{5,9}(?:\.\d+)?)\b"
)

kju_pattern = re.compile(
    r"\bKJU\d{4}\b",
    re.IGNORECASE
)


# ------------------------------------------------------------
# Search
# ------------------------------------------------------------

results = []

for idx, kju in enumerate(candidate_kjus, 1):

    print(
        f"[{idx:03d}/{len(candidate_kjus):03d}] Searching {kju} ..."
    )

    evidence = []

    for f in files:

        try:

            # Avoid enormous binary files.
            if f.stat().st_size > 50_000_000:
                continue

            data = f.read_bytes()

            # Skip obvious binary files.
            if b"\x00" in data[:10000]:
                continue

            text = data.decode(
                "utf-8",
                errors="ignore"
            )

        except Exception:
            continue

        if not re.search(
            rf"\b{re.escape(kju)}\b",
            text,
            flags=re.IGNORECASE
        ):
            continue

        # ----------------------------------------------------
        # Extract local context around KJU
        # ----------------------------------------------------

        contexts = []

        for m in re.finditer(
            rf"\b{re.escape(kju)}\b",
            text,
            flags=re.IGNORECASE
        ):

            start = max(0, m.start() - 500)
            end = min(len(text), m.end() + 1000)

            context = text[start:end]

            contexts.append(context)

        combined = "\n".join(contexts[:5])

        # ----------------------------------------------------
        # Find accession candidates in context
        # ----------------------------------------------------

        accessions = sorted(
            set(
                accession_pattern.findall(
                    combined
                )
            )
        )

        # Remove obvious KJU IDs from accession candidates.
        accessions = [
            x for x in accessions
            if not x.upper().startswith("KJU")
        ]

        evidence.append({
            "file": str(f.relative_to(MASTER)),
            "file_name": f.name,
            "accessions_found": ";".join(accessions),
            "context": combined[:3000],
        })

    # --------------------------------------------------------
    # Consolidate
    # --------------------------------------------------------

    unique_files = sorted(
        set(x["file"] for x in evidence)
    )

    accession_set = sorted(
        set(
            acc
            for e in evidence
            for acc in e["accessions_found"].split(";")
            if acc
        )
    )

    if not evidence:

        status = "NO_PROJECT_EVIDENCE"

    elif accession_set:

        status = "KJU_WITH_ACCESSION_EVIDENCE"

    else:

        status = "KJU_FOUND_NO_ACCESSION"

    results.append({
        "kju_id": kju,
        "status": status,
        "evidence_file_count": len(unique_files),
        "evidence_files": ";".join(unique_files),
        "accessions_found": ";".join(accession_set),
        "evidence_context": "\n---\n".join(
            e["context"][:1000]
            for e in evidence[:5]
        )
    })


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

result_df = pd.DataFrame(results)

result_df.to_csv(
    OUT,
    index=False
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary = []

summary.append(
    "DODDAGGUBBI - DEEP RESOLUTION OF UNRESOLVED KJUs"
)

summary.append("=" * 110)
summary.append("")

summary.append(
    f"Candidate unresolved KJUs: {len(candidate_kjus)}"
)

summary.append(
    f"With project evidence: {(result_df['evidence_file_count'] > 0).sum()}"
)

summary.append(
    f"With accession evidence: "
    f"{(result_df['accessions_found'].str.strip() != '').sum()}"
)

summary.append(
    f"No project evidence: "
    f"{(result_df['evidence_file_count'] == 0).sum()}"
)

summary.append("")

summary.append("RESULTS")
summary.append("-" * 110)

for _, r in result_df.iterrows():

    summary.append(
        f"{r['kju_id']} | "
        f"{r['status']} | "
        f"files={r['evidence_file_count']} | "
        f"accessions={r['accessions_found']}"
    )

summary.append("")
summary.append("No master dataset was modified.")

SUMMARY.write_text(
    "\n".join(summary),
    encoding="utf-8"
)


print("\n" + "=" * 110)
print("DEEP RESOLUTION COMPLETE")
print("=" * 110)

print("Output :", OUT)
print("Summary:", SUMMARY)

print("\nSTATUS COUNTS")
print(result_df["status"].value_counts().to_string())

print("\nACCESSION EVIDENCE")
print(
    result_df[
        result_df["accessions_found"].str.strip() != ""
    ][
        ["kju_id", "status", "accessions_found"]
    ].to_string(index=False)
)

print("\nMASTER DATA WAS NOT MODIFIED.")
