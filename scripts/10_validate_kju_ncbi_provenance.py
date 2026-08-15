from pathlib import Path
import pandas as pd
from Bio import Entrez, SeqIO
import time

MASTER = Path.home() / "Doddagubbi_Master"
OUT = MASTER / "07_Master_Metadata"

INPUT = OUT / "KJU_NCBI_RESOLVED.csv"
OUTPUT = OUT / "KJU_NCBI_PROVENANCE_VALIDATION.csv"

Entrez.email = "pranav@example.com"

df = pd.read_csv(INPUT)

results = []

print("=" * 80)
print(" DODDAGGUBBI - NCBI KJU PROVENANCE VALIDATION")
print("=" * 80)

for i, row in df.iterrows():

    kju = row["kju_id"]
    accession = row["accession"]

    print(f"\n[{i+1}/{len(df)}] {kju} -> {accession}")

    try:

        handle = Entrez.efetch(
            db="nuccore",
            id=accession,
            rettype="gb",
            retmode="text"
        )

        rec = SeqIO.read(handle, "genbank")
        handle.close()

        # ----------------------------------------------------
        # Collect all record text
        # ----------------------------------------------------

        all_text = []

        all_text.append(rec.description)
        all_text.append(str(rec.annotations))

        for ref in rec.annotations.get("references", []):
            all_text.append(str(ref))

        for feature in rec.features:
            all_text.append(str(feature.qualifiers))

        complete_text = " ".join(all_text)

        # ----------------------------------------------------
        # Extract qualifiers
        # ----------------------------------------------------

        strain = ""
        isolation_source = ""
        geo_loc_name = ""
        collection_date = ""
        host = ""
        note = ""
        submitter = ""

        for feature in rec.features:

            q = feature.qualifiers

            if not strain and "strain" in q:
                strain = "; ".join(q["strain"])

            if not isolation_source and "isolation_source" in q:
                isolation_source = "; ".join(q["isolation_source"])

            if not geo_loc_name and "geo_loc_name" in q:
                geo_loc_name = "; ".join(q["geo_loc_name"])

            if not collection_date and "collection_date" in q:
                collection_date = "; ".join(q["collection_date"])

            if not host and "host" in q:
                host = "; ".join(q["host"])

            if not note and "note" in q:
                note = "; ".join(q["note"])

            if "db_xref" in q:
                pass

        submitter = rec.annotations.get("taxonomy", "")

        # ----------------------------------------------------
        # Check KJU occurrence
        # ----------------------------------------------------

        kju_present = kju.lower() in complete_text.lower()

        # ----------------------------------------------------
        # 16S annotation check
        # ----------------------------------------------------

        is_16s = False

        for feature in rec.features:

            q = feature.qualifiers

            product = " ".join(
                q.get("product", [])
            ).lower()

            gene = " ".join(
                q.get("gene", [])
            ).lower()

            note_text = " ".join(
                q.get("note", [])
            ).lower()

            combined = (
                product + " " +
                gene + " " +
                note_text
            )

            if (
                "16s" in combined
                or "16s ribosomal" in combined
                or "small subunit ribosomal" in combined
            ):
                is_16s = True

        # ----------------------------------------------------
        # Geographic relevance
        # ----------------------------------------------------

        india_present = "india" in complete_text.lower()
        bengaluru_present = (
            "bengaluru" in complete_text.lower()
            or "bangalore" in complete_text.lower()
        )

        doddagubbi_present = (
            "doddagubbi" in complete_text.lower()
        )

        # ----------------------------------------------------
        # Print
        # ----------------------------------------------------

        print("  Organism          :", rec.annotations.get("organism", ""))
        print("  KJU in record     :", kju_present)
        print("  16S annotation    :", is_16s)
        print("  India mentioned   :", india_present)
        print("  Bengaluru mention:", bengaluru_present)
        print("  Doddagubbi mention:", doddagubbi_present)
        print("  Strain            :", strain)
        print("  Isolation source  :", isolation_source)
        print("  Geo location      :", geo_loc_name)
        print("  Collection date   :", collection_date)

        results.append({
            "kju_id": kju,
            "accession": accession,
            "organism": rec.annotations.get("organism", ""),
            "length": len(rec.seq),
            "kju_present": kju_present,
            "is_16s": is_16s,
            "india_present": india_present,
            "bengaluru_present": bengaluru_present,
            "doddagubbi_present": doddagubbi_present,
            "strain": strain,
            "isolation_source": isolation_source,
            "geo_loc_name": geo_loc_name,
            "collection_date": collection_date,
            "host": host,
            "note": note,
            "description": rec.description,
            "taxonomy": "; ".join(
                rec.annotations.get("taxonomy", [])
            )
        })

        time.sleep(0.4)

    except Exception as e:

        print("  ERROR:", e)

        results.append({
            "kju_id": kju,
            "accession": accession,
            "organism": "",
            "length": "",
            "kju_present": False,
            "is_16s": False,
            "india_present": False,
            "bengaluru_present": False,
            "doddagubbi_present": False,
            "strain": "",
            "isolation_source": "",
            "geo_loc_name": "",
            "collection_date": "",
            "host": "",
            "note": "",
            "description": "",
            "taxonomy": ""
        })

out = pd.DataFrame(results)

out.to_csv(
    OUTPUT,
    index=False
)

print("\n")
print("=" * 80)
print(" PROVENANCE VALIDATION COMPLETE")
print("=" * 80)

print("\nTotal records:", len(out))

print("\nKJU explicitly present:")
print(out["kju_present"].value_counts())

print("\n16S annotation:")
print(out["is_16s"].value_counts())

print("\nIndia:")
print(out["india_present"].value_counts())

print("\nBengaluru:")
print(out["bengaluru_present"].value_counts())

print("\nDoddagubbi:")
print(out["doddagubbi_present"].value_counts())

print("\nComplete validation table:")
print(
    out[
        [
            "kju_id",
            "accession",
            "organism",
            "length",
            "kju_present",
            "is_16s",
            "india_present",
            "bengaluru_present",
            "doddagubbi_present",
            "strain",
            "isolation_source",
            "geo_loc_name"
        ]
    ].to_string(index=False)
)

print("\nSaved:")
print(OUTPUT)

print("=" * 80)
