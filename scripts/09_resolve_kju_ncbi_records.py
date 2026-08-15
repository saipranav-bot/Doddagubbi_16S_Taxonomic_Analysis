from pathlib import Path
import pandas as pd
from Bio import Entrez, SeqIO
import time
import re

MASTER = Path.home() / "Doddagubbi_Master"
OUT = MASTER / "07_Master_Metadata"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# NCBI configuration
# ------------------------------------------------------------

Entrez.email = "pranav@example.com"

INPUT = OUT / "KJU_NCBI_FOUND_NOT_DOWNLOADED.csv"
OUTPUT = OUT / "KJU_NCBI_RESOLVED.csv"

df = pd.read_csv(INPUT)

print("=" * 80)
print(" DODDAGGUBBI - RESOLVE NCBI KJU RECORDS")
print("=" * 80)

print(f"\nKJU IDs to resolve: {len(df)}")

results = []

# ------------------------------------------------------------
# Search NCBI for each KJU ID
# ------------------------------------------------------------

for i, row in df.iterrows():

    kju = str(row["kju_id"]).strip()

    print(f"\n[{i+1}/{len(df)}] Searching {kju}...")

    try:

        # Search exact KJU ID
        handle = Entrez.esearch(
            db="nuccore",
            term=f'"{kju}"',
            retmax=20
        )

        search = Entrez.read(handle)
        handle.close()

        ids = search.get("IdList", [])

        if not ids:
            print("  No NCBI nucleotide record found.")
            results.append({
                "kju_id": kju,
                "ncbi_uid": "",
                "accession": "",
                "organism": "",
                "sequence_type": "",
                "length": "",
                "title": "",
                "isolation_source": "",
                "geo_loc_name": "",
                "collection_date": "",
                "strain": "",
                "matched": "NO"
            })
            continue

        # Retrieve records
        handle = Entrez.efetch(
            db="nuccore",
            id=",".join(ids),
            rettype="gb",
            retmode="text"
        )

        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()

        # Try to find the best KJU-matching record
        best = None

        for rec in records:

            text = (
                rec.description
                + " "
                + str(rec.annotations)
            )

            for feature in rec.features:
                text += " " + str(feature.qualifiers)

            if kju.lower() in text.lower():
                best = rec
                break

        if best is None:
            best = records[0]

        accession = best.id
        organism = best.annotations.get("organism", "")
        title = best.description

        # ----------------------------------------------------
        # Extract qualifiers
        # ----------------------------------------------------

        isolation_source = ""
        geo_loc_name = ""
        collection_date = ""
        strain = ""

        for feature in best.features:

            q = feature.qualifiers

            if "isolation_source" in q:
                isolation_source = "; ".join(q["isolation_source"])

            if "geo_loc_name" in q:
                geo_loc_name = "; ".join(q["geo_loc_name"])

            if "collection_date" in q:
                collection_date = "; ".join(q["collection_date"])

            if "strain" in q:
                strain = "; ".join(q["strain"])

        # ----------------------------------------------------
        # Determine sequence type
        # ----------------------------------------------------

        seq_type = ""

        full_text = (
            best.description.lower()
            + " "
            + str(best.annotations).lower()
        )

        for feature in best.features:

            if feature.type.lower() in [
                "rRNA",
                "rrna"
            ]:

                product = " ".join(
                    feature.qualifiers.get("product", [])
                ).lower()

                gene = " ".join(
                    feature.qualifiers.get("gene", [])
                ).lower()

                combined = product + " " + gene

                if "16s" in combined or "16s ribosomal" in combined:
                    seq_type = "16S"

                if "internal transcribed spacer" in combined or "its" in combined:
                    seq_type = "ITS"

        if not seq_type:

            if "16s" in full_text:
                seq_type = "16S"

            elif "internal transcribed spacer" in full_text:
                seq_type = "ITS"

        print(f"  Accession : {accession}")
        print(f"  Organism  : {organism}")
        print(f"  Type      : {seq_type}")
        print(f"  Length    : {len(best.seq)}")

        results.append({
            "kju_id": kju,
            "ncbi_uid": ids[0],
            "accession": accession,
            "organism": organism,
            "sequence_type": seq_type,
            "length": len(best.seq),
            "title": title,
            "isolation_source": isolation_source,
            "geo_loc_name": geo_loc_name,
            "collection_date": collection_date,
            "strain": strain,
            "matched": "YES"
        })

        time.sleep(0.4)

    except Exception as e:

        print(f"  ERROR: {e}")

        results.append({
            "kju_id": kju,
            "ncbi_uid": "",
            "accession": "",
            "organism": "",
            "sequence_type": "",
            "length": "",
            "title": "",
            "isolation_source": "",
            "geo_loc_name": "",
            "collection_date": "",
            "strain": "",
            "matched": "ERROR"
        })

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

result_df = pd.DataFrame(results)

result_df.to_csv(
    OUTPUT,
    index=False
)

print("\n")
print("=" * 80)
print(" NCBI RESOLUTION COMPLETE")
print("=" * 80)

print(f"\nRecords processed : {len(result_df)}")
print(f"Successfully matched : {(result_df['matched'] == 'YES').sum()}")
print(f"Unresolved : {(result_df['matched'] != 'YES').sum()}")

print("\nSequence types:")
print(
    result_df["sequence_type"]
    .replace("", "UNKNOWN")
    .value_counts()
)

print("\nResolved records:")
print(
    result_df[
        [
            "kju_id",
            "accession",
            "organism",
            "sequence_type",
            "length"
        ]
    ].to_string(index=False)
)

print("\nSaved:")
print(OUTPUT)

print("\n" + "=" * 80)
