from pathlib import Path
from Bio import Entrez, SeqIO
import pandas as pd
import time
import re

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

recon_file = META / "MASTER_KJU_RECONCILIATION_001_183.csv"

if not recon_file.exists():
    raise SystemExit("ERROR: reconciliation file not found")

df = pd.read_csv(recon_file)

missing = df[df["status"] != "DOWNLOADED"]["kju_id"].tolist()

print("=" * 80)
print(" DODDAGGUBBI - BROAD NCBI RECOVERY")
print("=" * 80)

print(f"\nMissing KJU IDs to investigate: {len(missing)}")

# ------------------------------------------------------------
# NCBI configuration
# ------------------------------------------------------------

Entrez.email = "pranav.bioinformatics@example.com"

results = []

# ------------------------------------------------------------
# Search strategies
# ------------------------------------------------------------

strategies = [
    '"{kju}"',
    '"{kju}" AND 16S',
    '"{kju}" AND bacteria',
    '"{kju}" AND "Doddagubbi"',
    '"{kju}" AND Bengaluru',
    '"{kju}" AND India',
]

for n, kju in enumerate(missing, 1):

    print(f"\n[{n}/{len(missing)}] {kju}")

    seen = set()

    for strategy in strategies:

        query = strategy.format(kju=kju)

        try:

            handle = Entrez.esearch(
                db="nuccore",
                term=query,
                retmax=20
            )

            result = Entrez.read(handle)
            handle.close()

            ids = [str(x) for x in result["IdList"]]

            for ncbi_id in ids:
                seen.add(ncbi_id)

            if ids:
                print(
                    f"  Query: {query} -> "
                    f"{len(ids)} hit(s)"
                )

            time.sleep(0.35)

        except Exception as e:

            print(f"  Search error: {e}")

    # --------------------------------------------------------
    # Fetch summaries for all unique hits
    # --------------------------------------------------------

    for ncbi_id in sorted(seen):

        try:

            handle = Entrez.efetch(
                db="nuccore",
                id=ncbi_id,
                rettype="gb",
                retmode="text"
            )

            text = handle.read()
            handle.close()

            if not text.strip():
                continue

            # Save temporary in-memory record
            import io

            record = SeqIO.read(
                io.StringIO(text),
                "genbank"
            )

            full_text = text.lower()

            results.append({
                "kju_id": kju,
                "ncbi_id": ncbi_id,
                "accession": record.id,
                "organism": record.annotations.get(
                    "organism", ""
                ),
                "length": len(record.seq),
                "description": record.description,
                "kju_present": (
                    kju.lower() in full_text
                ),
                "16s_present": (
                    "16s" in full_text or
                    "rrs" in full_text
                ),
                "india_present": (
                    "india" in full_text
                ),
                "bengaluru_present": (
                    "bengaluru" in full_text
                ),
                "doddagubbi_present": (
                    "doddagubbi" in full_text
                )
            })

            time.sleep(0.35)

        except Exception as e:

            print(
                f"  Fetch error {ncbi_id}: {e}"
            )


# ------------------------------------------------------------
# Save all candidates
# ------------------------------------------------------------

out = META / "KJU_BROAD_NCBI_CANDIDATES.csv"

if results:

    result_df = pd.DataFrame(results)

    result_df = result_df.drop_duplicates(
        subset=["kju_id", "accession"]
    )

    result_df.to_csv(
        out,
        index=False
    )

else:

    result_df = pd.DataFrame()

    result_df.to_csv(
        out,
        index=False
    )


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("\n")
print("=" * 80)
print(" BROAD NCBI SEARCH COMPLETE")
print("=" * 80)

print(f"Missing KJU IDs investigated : {len(missing)}")
print(f"Candidate records found      : {len(result_df)}")

if len(result_df):

    print(
        f"KJU IDs with candidates     : "
        f"{result_df.kju_id.nunique()}"
    )

    print(
        f"Candidates with KJU text    : "
        f"{result_df.kju_present.sum()}"
    )

    print(
        f"Candidates with 16S/rrs     : "
        f"{result_df['16s_present'].sum()}"
    )

    print(
        f"Candidates mentioning India : "
        f"{result_df['india_present'].sum()}"
    )

    print(
        f"Candidates mentioning "
        f"Doddagubbi                : "
        f"{result_df['doddagubbi_present'].sum()}"
    )

print(f"\nSaved:")
print(out)

print("\n")
print("=" * 80)
print(" STEP 12 COMPLETE")
print("=" * 80)
