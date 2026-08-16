#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
from Bio import SeqIO

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Doddagubbi_Master"

ROOT = (
    BASE
    / "github_results"
    / "reconciliation_111_vs_109"
)

TAXDIR = ROOT / "taxonomic_validation"

NEW_FASTA = (
    ROOT
    / "ncbi_check"
    / "OPENCITY_MISSING_KJU0224_KJU0225.fasta"
)

BLAST_FILE = (
    TAXDIR
    / "BLAST_NCBI_VALIDATION_KJU0224_KJU0225.csv"
)

PHYLO_FILE = (
    TAXDIR
    / "NCBI_BACKED_PHYLOGENETIC_PLACEMENT.csv"
)

OUTFILE = (
    TAXDIR
    / "FINAL_KJU0224_KJU0225_TAXONOMIC_VALIDATION.csv"
)

SUMMARY = (
    TAXDIR
    / "FINAL_KJU0224_KJU0225_TAXONOMIC_VALIDATION_SUMMARY.txt"
)

# ============================================================
# CLAIMED TAXONOMY
# ============================================================

claimed = {
    "KJU0224": "Pseudomonas putida",
    "KJU0225": "Chryseobacterium indologenes",
}

# ============================================================
# LOAD SEQUENCE LENGTHS
# ============================================================

sequence_lengths = {}

for record in SeqIO.parse(NEW_FASTA, "fasta"):
    sequence_lengths[record.id] = len(record.seq)

# ============================================================
# LOAD BLAST
# ============================================================

blast = pd.read_csv(BLAST_FILE)

# ============================================================
# LOAD PHYLOGENETIC DATA
# ============================================================

phylo = pd.read_csv(PHYLO_FILE)

# ============================================================
# NORMALIZE BOOLEAN COLUMNS
# ============================================================

phylo["genus_match"] = (
    phylo["genus_match"]
    .astype(str)
    .str.lower()
    .map({"true": True, "false": False})
)

phylo["species_match"] = (
    phylo["species_match"]
    .astype(str)
    .str.lower()
    .map({"true": True, "false": False})
)

# ============================================================
# FINAL RESULTS
# ============================================================

rows = []

for kju in ["KJU0224", "KJU0225"]:

    claim = claimed[kju]

    print()
    print("=" * 90)
    print(f"FINAL VALIDATION: {kju}")
    print("=" * 90)

    # --------------------------------------------------------
    # Sequence length
    # --------------------------------------------------------

    length = sequence_lengths.get(kju)

    # --------------------------------------------------------
    # BLAST
    # --------------------------------------------------------

    b = blast[
        blast["query"] == kju
    ].copy()

    if len(b) == 0:
        raise RuntimeError(
            f"No BLAST result found for {kju}"
        )

    # Sort strongest BLAST evidence
    sort_cols = []

    if "identity_percent" in b.columns:
        sort_cols.append("identity_percent")

    if "query_coverage_percent" in b.columns:
        sort_cols.append("query_coverage_percent")

    if "bit_score" in b.columns:
        sort_cols.append("bit_score")

    if sort_cols:
        b = b.sort_values(
            sort_cols,
            ascending=[False] * len(sort_cols)
        )

    best_blast = b.iloc[0]

    blast_accession = best_blast.get(
        "accession",
        ""
    )

    blast_identity = float(
        best_blast.get(
            "identity_percent",
            0
        )
    )

    blast_coverage = float(
        best_blast.get(
            "query_coverage_percent",
            0
        )
    )

    blast_evalue = best_blast.get(
        "evalue",
        ""
    )

    blast_title = best_blast.get(
        "title",
        ""
    )

    title_lower = str(
        blast_title
    ).lower()

    if "pseudomonas putida" in title_lower:
        blast_organism = "Pseudomonas putida"

    elif "chryseobacterium indologenes" in title_lower:
        blast_organism = "Chryseobacterium indologenes"

    else:
        blast_organism = str(
            best_blast.get(
                "organism",
                ""
            )
        )

    # --------------------------------------------------------
    # PHYLOGENETIC TOP 10
    # --------------------------------------------------------

    p = phylo[
        phylo["query"] == kju
    ].copy()

    if len(p) == 0:
        raise RuntimeError(
            f"No phylogenetic results found for {kju}"
        )

    # Closest references first
    p = p.sort_values(
        "distance",
        ascending=True
    )

    top10 = p.head(10)

    best_phylo = top10.iloc[0]

    best_reference = best_phylo[
        "reference_accession"
    ]

    best_phylo_organism = best_phylo[
        "reference_organism"
    ]

    best_distance = float(
        best_phylo["distance"]
    )

    genus_matches = int(
        top10["genus_match"]
        .fillna(False)
        .sum()
    )

    species_matches = int(
        top10["species_match"]
        .fillna(False)
        .sum()
    )

    # --------------------------------------------------------
    # PHYLOGENETIC INTERPRETATION
    # --------------------------------------------------------

    if (
        genus_matches == 10
        and species_matches >= 6
    ):
        phylo_interpretation = (
            "STRONG SPECIES-LEVEL PHYLOGENETIC SUPPORT"
        )

    elif genus_matches >= 8:
        phylo_interpretation = (
            "GENUS-LEVEL PHYLOGENETIC SUPPORT"
        )

    else:
        phylo_interpretation = (
            "PHYLOGENETIC PLACEMENT REQUIRES REVIEW"
        )

    # --------------------------------------------------------
    # TAXONOMIC AGREEMENT
    # --------------------------------------------------------

    blast_species_match = (
        blast_organism.lower()
        == claim.lower()
    )

    # --------------------------------------------------------
    # FINAL VERDICT
    # --------------------------------------------------------

    if (
        blast_species_match
        and blast_identity >= 99.0
        and blast_coverage >= 99.0
        and genus_matches == 10
        and species_matches >= 6
    ):

        verdict = "SUPPORTED"
        confidence = "VERY HIGH"

    elif (
        blast_species_match
        and blast_identity >= 97.0
        and blast_coverage >= 97.0
        and genus_matches >= 8
    ):

        verdict = "SUPPORTED"
        confidence = "HIGH"

    else:

        verdict = "REQUIRES REVIEW"
        confidence = "LOW/MODERATE"

    # --------------------------------------------------------
    # PRINT
    # --------------------------------------------------------

    print()
    print("CLAIMED TAXONOMY")
    print("-" * 90)
    print(claim)

    print()
    print("SEQUENCE")
    print("-" * 90)
    print(f"Length                  : {length} bp")

    print()
    print("NCBI BLAST")
    print("-" * 90)
    print(f"Best accession           : {blast_accession}")
    print(f"Organism                 : {blast_organism}")
    print(f"Identity                 : {blast_identity:.2f}%")
    print(f"Coverage                 : {blast_coverage:.2f}%")
    print(f"E-value                  : {blast_evalue}")

    print()
    print("PHYLOGENETIC PLACEMENT")
    print("-" * 90)
    print(f"Best neighbor            : {best_reference}")
    print(f"Reference organism       : {best_phylo_organism}")
    print(f"Distance                 : {best_distance:.6f}")
    print(f"Genus agreement          : {genus_matches}/10")
    print(f"Species agreement        : {species_matches}/10")
    print(f"Interpretation           : {phylo_interpretation}")

    print()
    print("FINAL ASSESSMENT")
    print("-" * 90)
    print(f"FINAL VERDICT            : {verdict}")
    print(f"CONFIDENCE               : {confidence}")

    # --------------------------------------------------------
    # STORE
    # --------------------------------------------------------

    rows.append({

        "KJU_ID": kju,

        "claimed_taxonomy": claim,

        "sequence_length_bp": length,

        "best_NCBI_BLAST_accession":
            blast_accession,

        "best_NCBI_BLAST_organism":
            blast_organism,

        "best_NCBI_BLAST_identity_percent":
            blast_identity,

        "best_NCBI_BLAST_coverage_percent":
            blast_coverage,

        "best_NCBI_BLAST_Evalue":
            blast_evalue,

        "best_phylogenetic_reference":
            best_reference,

        "best_phylogenetic_organism":
            best_phylo_organism,

        "phylogenetic_distance":
            best_distance,

        "phylogenetic_genus_matches_top10":
            genus_matches,

        "phylogenetic_species_matches_top10":
            species_matches,

        "phylogenetic_interpretation":
            phylo_interpretation,

        "final_taxonomic_verdict":
            verdict,

        "confidence":
            confidence,
    })

# ============================================================
# SAVE CSV
# ============================================================

final_df = pd.DataFrame(rows)

final_df.to_csv(
    OUTFILE,
    index=False
)

# ============================================================
# SAVE SUMMARY
# ============================================================

with open(SUMMARY, "w") as f:

    f.write("=" * 90 + "\n")
    f.write("FINAL TAXONOMIC VALIDATION REPORT\n")
    f.write("Doddagubbi 16S — KJU0224 / KJU0225\n")
    f.write("=" * 90 + "\n\n")

    f.write(
        "Evidence integrated:\n"
        "1. NCBI BLAST validation\n"
        "2. NCBI-backed phylogenetic placement\n"
        "3. Top-10 taxonomic agreement\n"
        "4. Sequence-level taxonomic consistency\n\n"
    )

    for _, row in final_df.iterrows():

        f.write("-" * 90 + "\n")
        f.write(
            f"{row['KJU_ID']}\n"
        )
        f.write("-" * 90 + "\n")

        f.write(
            f"Claimed taxonomy      : "
            f"{row['claimed_taxonomy']}\n"
        )

        f.write(
            f"Sequence length       : "
            f"{row['sequence_length_bp']} bp\n"
        )

        f.write(
            f"Best BLAST accession  : "
            f"{row['best_NCBI_BLAST_accession']}\n"
        )

        f.write(
            f"BLAST organism        : "
            f"{row['best_NCBI_BLAST_organism']}\n"
        )

        f.write(
            f"BLAST identity        : "
            f"{row['best_NCBI_BLAST_identity_percent']:.2f}%\n"
        )

        f.write(
            f"BLAST coverage        : "
            f"{row['best_NCBI_BLAST_coverage_percent']:.2f}%\n"
        )

        f.write(
            f"Phylogenetic neighbor : "
            f"{row['best_phylogenetic_reference']}\n"
        )

        f.write(
            f"Phylogenetic organism : "
            f"{row['best_phylogenetic_organism']}\n"
        )

        f.write(
            f"Phylogenetic distance : "
            f"{row['phylogenetic_distance']:.6f}\n"
        )

        f.write(
            f"Genus agreement       : "
            f"{row['phylogenetic_genus_matches_top10']}/10\n"
        )

        f.write(
            f"Species agreement     : "
            f"{row['phylogenetic_species_matches_top10']}/10\n"
        )

        f.write(
            f"Interpretation        : "
            f"{row['phylogenetic_interpretation']}\n"
        )

        f.write(
            f"FINAL VERDICT         : "
            f"{row['final_taxonomic_verdict']}\n"
        )

        f.write(
            f"CONFIDENCE            : "
            f"{row['confidence']}\n\n"
        )

    f.write("=" * 90 + "\n")
    f.write("END OF REPORT\n")
    f.write("=" * 90 + "\n")

# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("=" * 90)
print("STEP 70 COMPLETE")
print("=" * 90)

print()
print("FINAL RESULTS")
print("=" * 90)

for _, row in final_df.iterrows():

    print()
    print(
        f"{row['KJU_ID']} → "
        f"{row['claimed_taxonomy']}"
    )

    print(
        f"BLAST: "
        f"{row['best_NCBI_BLAST_identity_percent']:.2f}% "
        f"identity / "
        f"{row['best_NCBI_BLAST_coverage_percent']:.2f}% "
        f"coverage"
    )

    print(
        f"Phylogeny: "
        f"{row['phylogenetic_genus_matches_top10']}/10 "
        f"genus, "
        f"{row['phylogenetic_species_matches_top10']}/10 "
        f"species"
    )

    print(
        f"VERDICT: "
        f"{row['final_taxonomic_verdict']} "
        f"({row['confidence']})"
    )

print()
print("CSV:")
print(OUTFILE)

print()
print("SUMMARY:")
print(SUMMARY)

