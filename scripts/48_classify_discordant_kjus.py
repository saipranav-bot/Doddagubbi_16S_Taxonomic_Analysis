from pathlib import Path
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata" / "FINAL_109_PHYLOGENY_V1"

INPUT = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_KJU_REFERENCE_REVIEW_V1.csv"
)

OUT = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_KJU_CLASSIFICATION_V1.csv"
)

SUMMARY = (
    META /
    "REFERENCE_PHYLOGENY_V2/"
    "DISCORDANT_KJU_CLASSIFICATION_SUMMARY_V1.txt"
)

df = pd.read_csv(INPUT, dtype=str).fillna("")
df["distance"] = pd.to_numeric(df["distance"])

rows = []

for kju, sub in df.groupby("kju_id", sort=True):

    annotation = sub["current_annotation"].iloc[0]

    sub = sub.sort_values("distance").reset_index(drop=True)

    nearest = sub.iloc[0]

    current_species = annotation
    nearest_species = nearest["reference_species"]

    # Best reference matching the current annotation
    same = sub[
        sub["reference_species"].str.lower()
        == current_species.lower()
    ]

    # Best reference from a different species
    different = sub[
        sub["reference_species"].str.lower()
        != current_species.lower()
    ]

    if same.empty:
        same_distance = float("nan")
        same_accession = ""
    else:
        same_distance = same["distance"].min()
        same_accession = same.loc[
            same["distance"].idxmin(),
            "reference_accession"
        ]

    if different.empty:
        alt_distance = float("nan")
        alt_species = ""
        alt_accession = ""
    else:
        idx = different["distance"].idxmin()
        alt = different.loc[idx]

        alt_distance = float(alt["distance"])
        alt_species = alt["reference_species"]
        alt_accession = alt["reference_accession"]

    if pd.isna(same_distance) or pd.isna(alt_distance):
        delta = float("nan")
        ratio = float("nan")
        classification = "INSUFFICIENT_REFERENCE_EVIDENCE"

    else:
        delta = same_distance - alt_distance
        ratio = alt_distance / same_distance

        # Alternative species is actually closer
        if alt_distance < same_distance:

            if abs(delta) >= 0.002:
                classification = "STRONG_ALTERNATIVE_PLACEMENT"

            elif abs(delta) >= 0.0002:
                classification = "MODERATE_ALTERNATIVE_PLACEMENT"

            else:
                classification = "CLOSE_ALTERNATIVE_AMBIGUITY"

        else:

            if abs(delta) >= 0.002:
                classification = "CURRENT_ANNOTATION_SUPPORTED"

            elif abs(delta) >= 0.0002:
                classification = "CURRENT_ANNOTATION_PREFERRED"

            else:
                classification = "NEAR_EQUAL_SPECIES_DISTANCE"

    rows.append({
        "kju_id": kju,
        "current_annotation": current_species,
        "nearest_reference_species": nearest_species,
        "nearest_reference_accession": nearest["reference_accession"],
        "nearest_reference_distance": float(nearest["distance"]),
        "best_same_species_reference": current_species,
        "best_same_species_accession": same_accession,
        "best_same_species_distance": same_distance,
        "best_alternative_species": alt_species,
        "best_alternative_accession": alt_accession,
        "best_alternative_distance": alt_distance,
        "alternative_minus_current_distance": delta,
        "alternative_to_current_distance_ratio": ratio,
        "classification": classification
    })

out = pd.DataFrame(rows)

out.to_csv(OUT, index=False)

with open(SUMMARY, "w") as f:

    f.write("=" * 110 + "\n")
    f.write("DODDAGGUBBI - DISCORDANT KJU CLASSIFICATION\n")
    f.write("=" * 110 + "\n\n")

    f.write(f"Discordant KJUs reviewed: {len(out)}\n\n")

    f.write("CLASSIFICATION COUNTS\n")
    f.write("-" * 110 + "\n")

    for cls, count in (
        out["classification"]
        .value_counts()
        .items()
    ):
        f.write(f"{cls}: {count}\n")

    f.write("\n\n")

    f.write("DETAILED REVIEW\n")
    f.write("-" * 110 + "\n")

    for _, r in out.iterrows():

        f.write(
            f"{r['kju_id']} | "
            f"{r['current_annotation']} | "
            f"nearest={r['nearest_reference_species']} | "
            f"nearest_d={r['nearest_reference_distance']:.8f} | "
            f"current_species_d={r['best_same_species_distance']:.8f} | "
            f"alternative={r['best_alternative_species']} | "
            f"alternative_d={r['best_alternative_distance']:.8f} | "
            f"delta={r['alternative_minus_current_distance']:.8f} | "
            f"{r['classification']}\n"
        )

print("=" * 110)
print("DISCORDANT KJU CLASSIFICATION COMPLETE")
print("=" * 110)

print(
    out[
        [
            "kju_id",
            "current_annotation",
            "best_same_species_distance",
            "best_alternative_species",
            "best_alternative_distance",
            "alternative_minus_current_distance",
            "classification"
        ]
    ].to_string(index=False)
)

print("\nOutput:", OUT)
print("Summary:", SUMMARY)
