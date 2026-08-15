from pathlib import Path
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"

BASE = (
    MASTER
    / "07_Master_Metadata"
    / "FINAL_109_PHYLOGENY_V1"
)

REF = BASE / "REFERENCE_PHYLOGENY_V2"
OUTDIR = BASE / "FINAL_REPORT_V1"
OUTDIR.mkdir(parents=True, exist_ok=True)

ASSESSMENT = REF / "FINAL_INTEGRATED_9KJU_TAXONOMIC_ASSESSMENT_V1.csv"
TARGET_SENS = REF / "TARGET_SPECIES_MAIN_VS_SENSITIVITY_V1.csv"
CLADE_SUPPORT = REF / "DISCORDANT_CLADE_SUPPORT_REVIEW_V2.csv"

OUT = OUTDIR / "FINAL_DODDAGGUBBI_109_KJU_TAXONOMIC_REPORT.txt"

assessment = pd.read_csv(ASSESSMENT, dtype=str).fillna("")
sens = pd.read_csv(TARGET_SENS, dtype=str).fillna("")
support = pd.read_csv(CLADE_SUPPORT, dtype=str).fillna("")

lines = []

lines += [
    "=" * 100,
    "DODDAGGUBBI 16S rRNA - FINAL 109-KJU TAXONOMIC ASSESSMENT",
    "=" * 100,
    "",
    "DATASET",
    "-" * 100,
    "Final validated KJU records: 109",
    "GenBank validation: 109/109 PASS",
    "Sequence QC: 109/109 PASS",
    "Original MAFFT alignment: 109 x 1550 bp",
    "Trimmed alignment: 109 x 1345 bp",
    "Main phylogeny model: GTR+F+R4",
    "Main phylogeny support: 1000 UFBoot + 1000 SH-aLRT",
    "",
    "REFERENCE-AUGMENTED ANALYSIS",
    "-" * 100,
    "Independent NCBI reference sequences: 65",
    "Target species represented: 13",
    "Total reference-augmented taxa: 174",
    "Independent-reference trimmed alignment: 1369 bp",
    "Independent-reference phylogeny model: GTR+F+R3",
    "Independent-reference support: 1000 UFBoot + 1000 SH-aLRT",
    "",
    "SENSITIVITY",
    "-" * 100,
    "Target KJUs compared between main and sensitivity trees: 38",
    "Stable nearest-reference placements: 38",
    "Changed placements: 0",
    "",
    "FINAL 9-KJU REVIEW",
    "-" * 100,
]

for _, r in assessment.iterrows():
    lines.append(
        f"{r['kju_id']} | "
        f"{r['current_annotation']} | "
        f"alternative={r['best_alternative_species']} | "
        f"delta={r['alternative_minus_current_distance']} | "
        f"SH-aLRT={r['sh_alrt_percent']} | "
        f"UFBoot={r['ufboot_percent']} | "
        f"stable={r['placement_stable']} | "
        f"{r['final_category']} | "
        f"{r['final_recommendation']}"
    )

lines += [
    "",
    "INTERPRETATION",
    "-" * 100,
    "KJU0017 and KJU0217 are the highest-priority taxonomic review candidates.",
    "Both show strong alternative placement support and stable placement in the independent-reference sensitivity analysis.",
    "",
    "KJU0007, KJU0088, KJU0099, KJU0106 and KJU0221 show a repeated moderate alternative signal",
    "toward Citrobacter braakii and should be flagged for further taxonomic review without automatic reassignment.",
    "",
    "KJU0105 and KJU0208 show weak/unresolved alternative signal and should retain their current Bacillus cereus annotation.",
    "",
    "FINAL SCIENTIFIC CAUTION",
    "-" * 100,
    "16S rRNA evidence alone is not sufficient to establish definitive species reassignment.",
    "High-priority candidates should be checked against additional markers, whole-genome ANI/dDDH,",
    "and curated taxonomic/reference records before any change to the master annotation.",
    "",
    "MASTER DATA STATUS",
    "-" * 100,
    "Original master data were not modified.",
    "This report is an evidence/interpretation layer only.",
]

OUT.write_text("\n".join(lines) + "\n")

print("=" * 100)
print("FINAL TAXONOMIC REPORT GENERATED")
print("=" * 100)
print(f"Output: {OUT}")
print(f"9-KJU rows: {len(assessment)}")
print("Master dataset was NOT modified.")
