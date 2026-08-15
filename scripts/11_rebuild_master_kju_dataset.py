from pathlib import Path
import pandas as pd
from Bio import SeqIO
import re

MASTER = Path.home() / "Doddagubbi_Master"

GB_DIR = MASTER / "01_NCBI" / "genbank"
FASTA_DIR = MASTER / "01_NCBI" / "fasta"
META = MASTER / "07_Master_Metadata"

META.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print(" DODDAGGUBBI - MASTER KJU DATASET REBUILD")
print("=" * 80)

# ------------------------------------------------------------
# 1. Read ALL GenBank files in master directory
# ------------------------------------------------------------

print("\n[1/7] Reading all master GenBank records")
print("-" * 80)

records = []

for gb in sorted(GB_DIR.glob("*.gb")):

    try:
        rec = SeqIO.read(gb, "genbank")

        text = gb.read_text(errors="ignore")

        # Extract KJU ID
        matches = re.findall(r"KJU\d{4}", text, flags=re.I)

        kju_ids = sorted(set(x.upper() for x in matches))

        kju_id = kju_ids[0] if kju_ids else ""

        organism = rec.annotations.get("organism", "")

        records.append({
            "kju_id": kju_id,
            "accession": rec.id,
            "organism": organism,
            "length": len(rec.seq),
            "molecule_type": rec.annotations.get("molecule_type", ""),
            "description": rec.description,
            "genbank_file": str(gb),
            "sequence_type": "16S"
        })

    except Exception as e:
        print(f"ERROR: {gb.name}: {e}")


df = pd.DataFrame(records)

print(f"GenBank files scanned : {len(list(GB_DIR.glob('*.gb')))}")
print(f"Records successfully parsed : {len(df)}")
print(f"Records with KJU IDs : {(df['kju_id'] != '').sum()}")
print(f"Unique KJU IDs : {df.loc[df.kju_id != '', 'kju_id'].nunique()}")


# ------------------------------------------------------------
# 2. Remove accidental duplicate accessions
# ------------------------------------------------------------

print("\n[2/7] Checking accession duplicates")
print("-" * 80)

dup_acc = df[df.duplicated("accession", keep=False)]

print(f"Duplicate accession rows: {len(dup_acc)}")

df = df.drop_duplicates("accession").copy()


# ------------------------------------------------------------
# 3. Check duplicate KJU IDs
# ------------------------------------------------------------

print("\n[3/7] Checking KJU duplicates")
print("-" * 80)

kju_df = df[df["kju_id"] != ""].copy()

dup_kju = kju_df[kju_df.duplicated("kju_id", keep=False)]

print(f"Duplicate KJU rows: {len(dup_kju)}")

if len(dup_kju):
    print(dup_kju[
        ["kju_id", "accession", "organism", "length"]
    ].to_string(index=False))


# ------------------------------------------------------------
# 4. Save complete downloaded dataset
# ------------------------------------------------------------

print("\n[4/7] Saving master downloaded KJU dataset")
print("-" * 80)

master_out = META / "MASTER_KJU_16S_DATASET.csv"

kju_df.sort_values("kju_id").to_csv(
    master_out,
    index=False
)

print(f"Saved: {master_out}")


# ------------------------------------------------------------
# 5. Generate FASTA containing all KJU sequences
# ------------------------------------------------------------

print("\n[5/7] Building combined KJU FASTA")
print("-" * 80)

combined_fasta = FASTA_DIR / "MASTER_KJU_16S.fasta"

written = 0

with open(combined_fasta, "w") as out:

    for gb in sorted(GB_DIR.glob("*.gb")):

        try:

            rec = SeqIO.read(gb, "genbank")

            text = gb.read_text(errors="ignore")

            matches = re.findall(
                r"KJU\d{4}",
                text,
                flags=re.I
            )

            kju_ids = sorted(
                set(x.upper() for x in matches)
            )

            if not kju_ids:
                continue

            kju = kju_ids[0]

            out.write(
                f">{kju}|{rec.id}|{rec.annotations.get('organism','')}\n"
            )

            seq = str(rec.seq)

            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

            written += 1

        except Exception:
            pass


print(f"Sequences written: {written}")
print(f"Saved: {combined_fasta}")


# ------------------------------------------------------------
# 6. Build KJU 1–183 reconciliation
# ------------------------------------------------------------

print("\n[6/7] Building final KJU 1-183 reconciliation")
print("-" * 80)

expected = [
    f"KJU{i:04d}"
    for i in range(1, 184)
]

lookup = {}

for _, row in kju_df.iterrows():

    kju = row["kju_id"]

    if kju:
        lookup[kju] = row


rows = []

for kju in expected:

    if kju in lookup:

        row = lookup[kju]

        rows.append({
            "kju_id": kju,
            "status": "DOWNLOADED",
            "accession": row["accession"],
            "organism": row["organism"],
            "sequence_type": row["sequence_type"],
            "length": row["length"]
        })

    else:

        rows.append({
            "kju_id": kju,
            "status": "NOT_CURRENTLY_DOWNLOADED",
            "accession": "",
            "organism": "",
            "sequence_type": "",
            "length": ""
        })


recon = pd.DataFrame(rows)

recon_out = META / "MASTER_KJU_RECONCILIATION_001_183.csv"

recon.to_csv(
    recon_out,
    index=False
)

print(f"Saved: {recon_out}")


# ------------------------------------------------------------
# 7. Final summary
# ------------------------------------------------------------

print("\n")
print("=" * 80)
print(" FINAL MASTER DATASET STATUS")
print("=" * 80)

print(f"Expected KJU IDs              : 183")
print(f"KJU IDs downloaded            : {sum(recon.status == 'DOWNLOADED')}")
print(f"KJU IDs not currently downloaded : {sum(recon.status != 'DOWNLOADED')}")
print(f"Unique organisms               : {kju_df.organism.nunique()}")

print("\nStatus:")
print(
    recon["status"]
    .value_counts()
    .to_string()
)

print("\nOutput files:")
print(master_out)
print(combined_fasta)
print(recon_out)

print("\n")
print("=" * 80)
print(" STEP 11 COMPLETE")
print("=" * 80)
