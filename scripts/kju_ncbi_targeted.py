
from Bio import Entrez
from pathlib import Path
import pandas as pd
import time

PROJECT = Path.home() / "Doddagubbi_Master"

# NCBI requires a real identifying email.
# This is only sent as an API contact field.
Entrez.email = "pranav.research@example.com"

outfile = (
    PROJECT /
    "07_Master_Metadata" /
    "KJU_NCBI_targeted_hits.csv"
)

rows = []

for i in range(1, 184):

    kju = f"KJU{i:04d}"

    query = f'"{kju}"'

    try:

        with Entrez.esearch(
            db="nucleotide",
            term=query,
            retmax=100
        ) as handle:

            result = Entrez.read(handle)

        ids = list(
            result["IdList"]
        )

        if ids:

            print(
                f"{kju}: {len(ids)} hit(s)"
            )

            for ncbi_id in ids:

                rows.append(
                    {
                        "kju_id": kju,
                        "ncbi_id": ncbi_id
                    }
                )

        else:

            print(
                f"{kju}: no hit"
            )

    except Exception as e:

        print(
            f"{kju}: ERROR: {e}"
        )

    time.sleep(0.34)

df = pd.DataFrame(rows)

if len(df):

    df = df.drop_duplicates()

else:

    df = pd.DataFrame(
        columns=[
            "kju_id",
            "ncbi_id"
        ]
    )

df.to_csv(
    outfile,
    index=False
)

print()
print("=" * 60)
print("TARGETED NCBI SEARCH COMPLETE")
print("=" * 60)

print(
    "NCBI hits:",
    len(df)
)

print(
    "KJU IDs with hits:",
    df["kju_id"].nunique()
    if len(df)
    else 0
)

print(
    "Saved:",
    outfile
)

