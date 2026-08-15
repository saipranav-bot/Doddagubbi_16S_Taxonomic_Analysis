from pathlib import Path
import pandas as pd
import requests
import re
import time

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

OUT = META / "KJU0039_INVESTIGATION_V1.csv"
SUMMARY = META / "KJU0039_INVESTIGATION_SUMMARY_V1.txt"

NCBI_ESearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

print("=" * 90)
print(" DODDAGGUBBI - KJU0039 IDENTITY INVESTIGATION")
print("=" * 90)

records = []

queries = [
    '"KJU0039"',
    '"KJU0039"[All Fields]',
    '"Doddagubbi" AND "KJU0039"',
    '"KJU0039" AND 16S',
]

for query in queries:
    print(f"\nNCBI SEARCH: {query}")

    params = {
        "db": "nuccore",
        "term": query,
        "retmode": "json",
        "retmax": 100,
    }

    try:
        r = requests.get(NCBI_ESearch, params=params, timeout=30)
        r.raise_for_status()
        data = r.json()

        ids = data.get("esearchresult", {}).get("idlist", [])

        print(f"Hits: {len(ids)}")

        for uid in ids:
            records.append({
                "query": query,
                "ncbi_uid": uid
            })

        time.sleep(0.35)

    except Exception as e:
        print(f"ERROR: {e}")

# Deduplicate
records_df = pd.DataFrame(records).drop_duplicates()

# ------------------------------------------------------------
# Fetch summaries
# ------------------------------------------------------------

summary_rows = []

if not records_df.empty:

    unique_ids = records_df["ncbi_uid"].unique().tolist()

    for uid in unique_ids:

        try:
            params = {
                "db": "nuccore",
                "id": uid,
                "retmode": "json",
            }

            r = requests.get(
                NCBI_ESummary,
                params=params,
                timeout=30
            )

            r.raise_for_status()

            data = r.json()
            result = data.get("result", {})
            item = result.get(uid, {})

            summary_rows.append({
                "ncbi_uid": uid,
                "accession": item.get("caption", ""),
                "title": item.get("title", ""),
                "organism": item.get("organism", ""),
                "taxid": item.get("taxid", ""),
                "length": item.get("slen", ""),
                "status": item.get("status", ""),
            })

            time.sleep(0.35)

        except Exception as e:
            print(f"ERROR fetching UID {uid}: {e}")

summary_df = pd.DataFrame(summary_rows)

# ------------------------------------------------------------
# Save evidence
# ------------------------------------------------------------

if summary_df.empty:
    summary_df = pd.DataFrame(columns=[
        "ncbi_uid",
        "accession",
        "title",
        "organism",
        "taxid",
        "length",
        "status"
    ])

summary_df.to_csv(OUT, index=False)

# ------------------------------------------------------------
# Additional local evidence
# ------------------------------------------------------------

local_hits = []

for path in MASTER.rglob("*"):

    if not path.is_file():
        continue

    if "venv" in path.parts:
        continue

    try:
        text = path.read_text(
            encoding="utf-8",
            errors="ignore"
        )
    except:
        continue

    if "KJU0039" in text:
        local_hits.append(str(path))

# ------------------------------------------------------------
# Classification
# ------------------------------------------------------------

ncbi_kju0039_hits = summary_df[
    summary_df["title"].astype(str).str.contains(
        "KJU0039",
        case=False,
        na=False
    )
]

ncbi_pz_hits = summary_df[
    summary_df["accession"].astype(str).str.contains(
        "PZ666011",
        case=False,
        na=False
    )
]

if not ncbi_kju0039_hits.empty:
    conclusion = "KJU0039_HAS_NCBI_RECORD"

elif not ncbi_pz_hits.empty:
    conclusion = "PZ666011_EXISTS_BUT_IS_NOT_KJU0039"

else:
    conclusion = "NO_NCBI_KJU0039_RECORD_FOUND"

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

lines = []

lines.append("DODDAGGUBBI - KJU0039 IDENTITY INVESTIGATION")
lines.append("=" * 80)
lines.append("")
lines.append(f"NCBI nucleotide records found: {len(summary_df)}")
lines.append(f"Records explicitly mentioning KJU0039: {len(ncbi_kju0039_hits)}")
lines.append(f"Records involving PZ666011: {len(ncbi_pz_hits)}")
lines.append("")
lines.append("CONCLUSION")
lines.append("-" * 80)
lines.append(conclusion)
lines.append("")
lines.append("NCBI RECORDS")
lines.append("-" * 80)

if summary_df.empty:
    lines.append("No NCBI records returned.")

else:
    for _, row in summary_df.iterrows():
        lines.append(
            f"{row['accession']} | "
            f"{row['organism']} | "
            f"{row['title']}"
        )

lines.append("")
lines.append("LOCAL FILES CONTAINING KJU0039")
lines.append("-" * 80)

for p in local_hits:
    lines.append(p)

lines.append("")
lines.append("INTERPRETATION")
lines.append("-" * 80)

if conclusion == "KJU0039_HAS_NCBI_RECORD":
    lines.append(
        "KJU0039 has independent NCBI evidence. "
        "The accession should be reviewed against that record."
    )

elif conclusion == "PZ666011_EXISTS_BUT_IS_NOT_KJU0039":
    lines.append(
        "PZ666011 is associated with another KJU identifier "
        "(currently KJU0038). Do NOT assign PZ666011 to KJU0039."
    )

else:
    lines.append(
        "No independent NCBI nucleotide record for KJU0039 was found "
        "using the tested queries. KJU0039 remains unresolved."
    )

lines.append("")
lines.append("MASTER DATA WAS NOT MODIFIED.")

SUMMARY.write_text(
    "\n".join(lines),
    encoding="utf-8"
)

print("\n" + "=" * 90)
print(" INVESTIGATION COMPLETE")
print("=" * 90)

print(f"\nEvidence table:")
print(OUT)

print(f"\nSummary:")
print(SUMMARY)

print("\nMASTER DATA WAS NOT MODIFIED.")

