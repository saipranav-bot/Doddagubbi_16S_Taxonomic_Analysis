from pathlib import Path
import pandas as pd

MASTER = Path.home() / "Doddagubbi_Master"

BASE = (
    MASTER /
    "07_Master_Metadata/FINAL_109_PHYLOGENY_V1/"
    "REFERENCE_PHYLOGENY_V2"
)

CLASSIFICATION = BASE / "DISCORDANT_KJU_CLASSIFICATION_V1.csv"
SUPPORT = BASE / "DISCORDANT_CLADE_SUPPORT_REVIEW_V2.csv"
TARGETS = BASE / "TARGET_SPECIES_MAIN_VS_SENSITIVITY_V1.csv"

OUT = BASE / "FINAL_DISCORDANT_KJU_EVIDENCE_V1.csv"
SUMMARY = BASE / "FINAL_DISCORDANT_KJU_EVIDENCE_SUMMARY_V1.txt"

cls = pd.read_csv(CLASSIFICATION, dtype=str).fillna("")
sup = pd.read_csv(SUPPORT, dtype=str).fillna("")
targets = pd.read_csv(TARGETS, dtype=str).fillna("")

# Convert numeric columns
for df, cols in [
    (cls, [
        "best_same_species_distance",
        "best_alternative_distance",
        "alternative_minus_current_distance"
    ]),
    (sup, [
        "sh_alrt_percent",
        "ufboot_percent"
    ]),
    (targets, [
        "main_nearest_reference_distance",
        "sensitivity_nearest_reference_distance"
    ])
]:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

merged = cls.merge(
    sup[
        [
            "kju_id",
            "sh_alrt_percent",
            "ufboot_percent",
            "classification"
        ]
    ],
    on="kju_id",
    how="left",
    suffixes=("", "_support")
)

# Rename support classification
merged = merged.rename(
    columns={
        "classification_support":
            "support_classification"
    }
)

# Add stability information
if "placement_stable" in targets.columns:

    merged = merged.merge(
        targets[
            [
                "kju_id",
                "placement_stable",
                "annotation_supported_in_both"
            ]
        ],
        on="kju_id",
        how="left"
    )

def final_category(row):

    ub = row.get("ufboot_percent")

    if pd.isna(ub):
        return "SUPPORT_UNAVAILABLE"

    if ub >= 95:
        return "STRONG_ALTERNATIVE_SIGNAL"

    if ub >= 80:
        return "MODERATE_ALTERNATIVE_SIGNAL"

    return "WEAK_ALTERNATIVE_SIGNAL"


merged["final_category"] = merged.apply(
    final_category,
    axis=1
)

# More conservative interpretation
def interpretation(row):

    cat = row["final_category"]

    if cat == "STRONG_ALTERNATIVE_SIGNAL":
        return (
            "Strong 16S phylogenetic evidence for alternative "
            "species placement; retain annotation provisionally "
            "pending independent taxonomic evidence."
        )

    if cat == "MODERATE_ALTERNATIVE_SIGNAL":
        return (
            "Moderate 16S phylogenetic evidence for alternative "
            "placement; taxonomic annotation requires further "
            "validation."
        )

    if cat == "WEAK_ALTERNATIVE_SIGNAL":
        return (
            "Weak alternative signal; insufficient evidence to "
            "challenge current annotation."
        )

    return "Insufficient support information."

merged["final_interpretation"] = merged.apply(
    interpretation,
    axis=1
)

# Sort strongest first
order = {
    "STRONG_ALTERNATIVE_SIGNAL": 1,
    "MODERATE_ALTERNATIVE_SIGNAL": 2,
    "WEAK_ALTERNATIVE_SIGNAL": 3,
    "SUPPORT_UNAVAILABLE": 4
}

merged["_order"] = merged["final_category"].map(order)
merged = merged.sort_values(
    ["_order", "kju_id"]
).drop(columns="_order")

merged.to_csv(
    OUT,
    index=False
)

with open(SUMMARY, "w") as f:

    f.write("=" * 110 + "\n")
    f.write(
        "DODDAGGUBBI - FINAL DISCORDANT KJU EVIDENCE REVIEW\n"
    )
    f.write("=" * 110 + "\n\n")

    f.write(
        f"Discordant KJUs reviewed: {len(merged)}\n\n"
    )

    f.write("FINAL EVIDENCE COUNTS\n")
    f.write("-" * 110 + "\n")

    for category, count in (
        merged["final_category"].value_counts().items()
    ):
        f.write(f"{category}: {count}\n")

    f.write("\nDETAILED RESULTS\n")
    f.write("-" * 110 + "\n")

    for _, r in merged.iterrows():

        f.write(
            f"{r['kju_id']} | "
            f"{r['current_annotation']} | "
            f"alternative={r['best_alternative_species']} | "
            f"current_dist={r['best_same_species_distance']:.6f} | "
            f"alternative_dist={r['best_alternative_distance']:.6f} | "
            f"delta={r['alternative_minus_current_distance']:.6f} | "
            f"SH-aLRT={r['sh_alrt_percent']} | "
            f"UFBoot={r['ufboot_percent']} | "
            f"{r['final_category']}\n"
        )

    f.write("\nIMPORTANT\n")
    f.write("-" * 110 + "\n")
    f.write(
        "These classifications are 16S-based evidence categories.\n"
        "They do NOT automatically rename KJU records or modify the\n"
        "canonical 109-record dataset.\n"
    )

print("=" * 110)
print("FINAL DISCORDANT KJU EVIDENCE TABLE COMPLETE")
print("=" * 110)

print(
    merged[
        [
            "kju_id",
            "current_annotation",
            "best_alternative_species",
            "best_same_species_distance",
            "best_alternative_distance",
            "alternative_minus_current_distance",
            "sh_alrt_percent",
            "ufboot_percent",
            "final_category"
        ]
    ].to_string(index=False)
)

print("\nOutput:", OUT)
print("Summary:", SUMMARY)
