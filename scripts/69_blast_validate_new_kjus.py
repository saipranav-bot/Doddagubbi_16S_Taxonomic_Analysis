#!/usr/bin/env python3

from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import time
import sys

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

NEW_FASTA = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "ncbi_check"
    / "OPENCITY_MISSING_KJU0224_KJU0225.fasta"
)

OUTDIR = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
    / "taxonomic_validation"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "BLAST_NCBI_VALIDATION_KJU0224_KJU0225.csv"

XMLDIR = OUTDIR / "blast_xml"
XMLDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# SETTINGS
# ============================================================

DATABASE = "nt"

# 16S sequences
PROGRAM = "blastn"

# Request reasonably stringent hits
EXPECT = 0.001

# Maximum number of returned alignments
HITLIST_SIZE = 20

# ============================================================
# BLAST FUNCTION
# ============================================================

def blast_search(sequence, kju_id):

    print()
    print("=" * 80)
    print(f"BLAST VALIDATION: {kju_id}")
    print("=" * 80)

    print(f"Sequence length: {len(sequence)}")
    print("Submitting BLAST search to NCBI...")
    print()
    print("This may take several minutes.")

    try:

        handle = NCBIWWW.qblast(
            PROGRAM,
            DATABASE,
            sequence,
            expect=EXPECT,
            hitlist_size=HITLIST_SIZE,
            descriptions=HITLIST_SIZE,
            alignments=HITLIST_SIZE,
            format_type="XML"
        )

        xml_file = XMLDIR / f"{kju_id}_blast.xml"

        with open(xml_file, "w") as out:
            out.write(handle.read())

        handle.close()

        print(f"BLAST XML saved:")
        print(xml_file)

    except Exception as e:

        print()
        print("ERROR DURING NCBI BLAST SUBMISSION")
        print(str(e))

        return []

    # ========================================================
    # PARSE XML
    # ========================================================

    print()
    print("Parsing BLAST results...")

    try:

        with open(xml_file) as handle:
            blast_record = NCBIXML.read(handle)

    except Exception as e:

        print("ERROR PARSING BLAST XML")
        print(str(e))

        return []

    results = []

    # ========================================================
    # TOP HITS
    # ========================================================

    print()
    print("-" * 80)
    print("TOP BLAST HITS")
    print("-" * 80)

    for rank, alignment in enumerate(
        blast_record.alignments[:HITLIST_SIZE],
        start=1
    ):

        hsp = alignment.hsps[0]

        title = alignment.title
        accession = alignment.accession

        identity_percent = (
            hsp.identities / hsp.align_length
        ) * 100

        coverage_percent = (
            hsp.align_length / len(sequence)
        ) * 100

        evalue = hsp.expect

        print(
            f"{rank:2d}. "
            f"{accession:<15} "
            f"identity={identity_percent:6.2f}% "
            f"coverage={coverage_percent:6.2f}% "
            f"E={evalue:.2e}"
        )

        results.append({
            "query": kju_id,
            "rank": rank,
            "accession": accession,
            "title": title,
            "identity_percent": round(identity_percent, 4),
            "alignment_length": hsp.align_length,
            "query_length": len(sequence),
            "query_coverage_percent": round(
                coverage_percent,
                4
            ),
            "mismatches": hsp.align_length - hsp.identities,
            "gaps": hsp.gaps,
            "evalue": evalue,
            "bit_score": hsp.bits,
        })

    return results


# ============================================================
# LOAD NEW SEQUENCES
# ============================================================

records = list(
    SeqIO.parse(
        NEW_FASTA,
        "fasta"
    )
)

print("=" * 80)
print("STEP 69 — NCBI BLAST VALIDATION")
print("=" * 80)

print(f"Sequences found: {len(records)}")

if len(records) == 0:

    print()
    print("ERROR: No sequences found.")
    sys.exit(1)

# ============================================================
# RUN BLAST
# ============================================================

ALL_RESULTS = []

for record in records:

    results = blast_search(
        str(record.seq),
        record.id
    )

    ALL_RESULTS.extend(results)

    # Avoid hammering NCBI
    print()
    print("Waiting before next NCBI request...")
    time.sleep(10)

# ============================================================
# SAVE RESULTS
# ============================================================

if len(ALL_RESULTS) == 0:

    print()
    print("=" * 80)
    print("NO BLAST RESULTS WERE RETURNED")
    print("=" * 80)
    sys.exit(1)

df = pd.DataFrame(ALL_RESULTS)

df.to_csv(
    OUTFILE,
    index=False
)

# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 80)
print("STEP 69 COMPLETE")
print("=" * 80)

print()
print("BLAST results saved:")
print(OUTFILE)

print()
print("BEST NCBI HIT PER KJU")
print("=" * 80)

for kju in df["query"].unique():

    sub = df[df["query"] == kju]

    best = sub.sort_values(
        [
            "identity_percent",
            "query_coverage_percent"
        ],
        ascending=False
    ).iloc[0]

    print()
    print(kju)
    print("-" * 80)
    print(f"Accession : {best['accession']}")
    print(f"Identity  : {best['identity_percent']:.2f}%")
    print(f"Coverage  : {best['query_coverage_percent']:.2f}%")
    print(f"E-value   : {best['evalue']:.2e}")
    print(f"Bit score : {best['bit_score']}")
    print(f"Title     : {best['title']}")

print()
print("=" * 80)
print("DONE")
print("=" * 80)

