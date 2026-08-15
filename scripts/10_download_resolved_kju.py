from pathlib import Path
from Bio import Entrez, SeqIO
import pandas as pd
import time
import shutil

MASTER = Path.home() / "Doddagubbi_Master"
SOURCE = Path.home() / "Doddagubbi_16S"

META = MASTER / "07_Master_Metadata"
GB_DIR = MASTER / "01_NCBI" / "genbank"
FASTA_DIR = MASTER / "01_NCBI" / "fasta"
LOG_DIR = MASTER / "logs"

GB_DIR.mkdir(parents=True, exist_ok=True)
FASTA_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

resolved = META / "KJU_NCBI_RESOLVED.csv"

if not resolved.exists():
    raise SystemExit(f"ERROR: {resolved} not found")

df = pd.read_csv(resolved)

print("=" * 80)
print(" DODDAGGUBBI - DOWNLOAD RESOLVED KJU RECORDS")
print("=" * 80)

print(f"\nResolved records found: {len(df)}")

# ------------------------------------------------------------
# NCBI configuration
# ------------------------------------------------------------

Entrez.email = "pranav.bioinformatics@example.com"

downloaded = []
failed = []

# ------------------------------------------------------------
# Download records
# ------------------------------------------------------------

for i, row in df.iterrows():

    kju = str(row["kju_id"])
    accession = str(row["accession"])

    print(f"\n[{i+1}/{len(df)}] {kju} -> {accession}")

    gb_out = GB_DIR / f"{accession}.gb"
    fa_out = FASTA_DIR / f"{accession}.fasta"

    try:

        # Fetch GenBank
        handle = Entrez.efetch(
            db="nuccore",
            id=accession,
            rettype="gb",
            retmode="text"
        )

        gb_text = handle.read()
        handle.close()

        if not gb_text.strip():
            raise RuntimeError("Empty GenBank response")

        gb_out.write_text(gb_text)

        # Parse downloaded GenBank
        record = SeqIO.read(gb_out, "genbank")

        # Write FASTA
        SeqIO.write(record, fa_out, "fasta")

        downloaded.append({
            "kju_id": kju,
            "accession": accession,
            "organism": row.get("organism", ""),
            "sequence_type": row.get("sequence_type", ""),
            "length": len(record.seq),
            "genbank": str(gb_out),
            "fasta": str(fa_out),
            "status": "DOWNLOADED"
        })

        print(f"  OK: {record.id}")
        print(f"  Organism: {record.annotations.get('organism', '')}")
        print(f"  Length: {len(record.seq)}")

        time.sleep(0.35)

    except Exception as e:

        print(f"  ERROR: {e}")

        failed.append({
            "kju_id": kju,
            "accession": accession,
            "error": str(e),
            "status": "FAILED"
        })


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

downloaded_df = pd.DataFrame(downloaded)
failed_df = pd.DataFrame(failed)

downloaded_out = META / "KJU_RESOLVED_DOWNLOADED.csv"
failed_out = META / "KJU_RESOLVED_DOWNLOAD_FAILED.csv"

downloaded_df.to_csv(downloaded_out, index=False)

if len(failed_df):
    failed_df.to_csv(failed_out, index=False)


# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

print("\n")
print("=" * 80)
print(" DOWNLOAD COMPLETE")
print("=" * 80)

print(f"Records requested : {len(df)}")
print(f"Successfully downloaded : {len(downloaded_df)}")
print(f"Failed : {len(failed_df)}")

print("\nSaved:")
print(downloaded_out)

if len(failed_df):
    print(failed_out)

print("=" * 80)
