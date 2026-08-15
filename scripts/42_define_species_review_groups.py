from pathlib import Path
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

INPUT = (
    META
    / "FINAL_109_PHYLOGENY_V1"
    / "SPECIES_RESOLUTION_REVIEW_V1.csv"
)

OUT = (
    META
    / "FINAL_109_PHYLOGENY_V1"
    / "SPECIES_REVIEW_TARGETS_V1.csv"
)

TARGETS = {
    "Aeromonas caviae",
    "Aeromonas hydrophila",
    "Aeromonas veronii",
    "Bacillus cereus",
    "Citrobacter freundii",
    "Shewanella putrefaciens",
}

df = pd.read_csv(INPUT, dtype=str).fillna("")

targets = df[
    df["species"].isin(TARGETS)
].copy()

targets.to_csv(
    OUT,
    index=False
)

print("=" * 100)
print("SPECIES REVIEW TARGETS")
print("=" * 100)

print(
    targets[
        [
            "species",
            "kju_count",
            "kju_ids",
            "nearest_different_species",
            "nearest_different_species_kju",
            "nearest_different_species_distance",
            "classification"
        ]
    ].to_string(index=False)
)

print("")
print("Output:", OUT)
