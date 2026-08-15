from pathlib import Path
import pandas as pd
import numpy as np

MASTER = Path.home() / "Doddagubbi_Master"

BASE = (
    MASTER
    / "07_Master_Metadata"
    / "FINAL_109_PHYLOGENY_V1"
    / "REFERENCE_PHYLOGENY_V2"
)

CLASSIFICATION = BASE / "DISCORDANT_KJU_CLASSIFICATION_V1.csv"
SUPPORT = BASE / "DISCORDANT_CLADE_SUPPORT_REVIEW_V2.csv"
SENSITIVITY = BASE / "TARGET_SPECIES_MAIN_VS_SENSITIVITY_V1.csv"
EVIDENCE = BASE / "FINAL_DISCORDANT_KJU_EVIDENCE_V1.csv"

OUT = BASE / "FINAL_INTEGRATED_9KJU_TAXONOMIC_ASSESSMENT_V1.csv"
SUMMARY = BASE / "FINAL_INTEGRATED_9KJU_TAXONOMIC_ASSESSMENT_SUMMARY_V1.txt"

TARGETS = {
    "KJU0007",
    "KJU0017",
    "KJU0088",
    "KJU0099",
    "KJU0105",
    "KJU0106",
    "KJU0208",
    "KJU0217",
    "KJU0221",
}


def load(path):
    if not path.exists():
        raise SystemExit(f"ERROR: Missing file:\n{path}")
    return pd.read_csv(path, dtype=str).fillna("")


def num(df, col):
    if col not in df.columns:
        return pd.Series(np.nan, index=df.index)
    return pd.to_numeric(
        df[col].astype(str).str.strip(),
        errors="coerce"
    )


print("=" * 110)
print(" DODDAGGUBBI - FINAL INTEGRATED 9-KJU TAXONOMIC ASSESSMENT")
print("=" * 110)

classification = load(CLASSIFICATION)
support = load(SUPPORT)
sensitivity = load(SENSITIVITY)
evidence = load(EVIDENCE)

# ------------------------------------------------------------
# Validate target set
# ------------------------------------------------------------

def validate_targets(df, name):
    found = set(df["kju_id"].astype(str))
    missing = TARGETS - found
    if missing:
        raise SystemExit(
            f"ERROR: {name} is missing KJUs: {sorted(missing)}"
        )

validate_targets(classification, "classification table")
validate_targets(support, "support table")
validate_targets(sensitivity, "sensitivity table")
validate_targets(evidence, "evidence table")

# ------------------------------------------------------------
# Select only target records
# ------------------------------------------------------------

classification = classification[
    classification["kju_id"].isin(TARGETS)
].copy()

support = support[
    support["kju_id"].isin(TARGETS)
].copy()

sensitivity = sensitivity[
    sensitivity["kju_id"].isin(TARGETS)
].copy()

evidence = evidence[
    evidence["kju_id"].isin(TARGETS)
].copy()

# ------------------------------------------------------------
# Merge evidence layers
# ------------------------------------------------------------

df = evidence.merge(
    support,
    on="kju_id",
    how="left",
    suffixes=("", "_support")
)

df = df.merge(
    sensitivity,
    on="kju_id",
    how="left",
    suffixes=("", "_sensitivity")
)

df = df.merge(
    classification[
        [
            "kju_id",
            "current_annotation",
            "best_same_species_distance",
            "best_alternative_species",
            "best_alternative_distance",
            "alternative_minus_current_distance",
        ]
    ],
    on="kju_id",
    how="left",
    suffixes=("", "_classification")
)

# ------------------------------------------------------------
# Resolve annotation fields
# ------------------------------------------------------------

for col in ["current_annotation", "best_alternative_species"]:
    c1 = f"{col}_classification"
    if c1 in df.columns:
        df[col] = df[col].replace("", np.nan).fillna(df.get(col, ""))

# ------------------------------------------------------------
# Numeric columns
# ------------------------------------------------------------

numeric_cols = [
    "best_same_species_distance",
    "best_alternative_distance",
    "alternative_minus_current_distance",
    "sh_alrt_percent",
    "ufboot_percent",
    "main_nearest_reference_distance",
    "sensitivity_nearest_reference_distance",
]

for col in numeric_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# ------------------------------------------------------------
# Sensitivity result
# ------------------------------------------------------------

if "placement_stable" not in df.columns:
    df["placement_stable"] = ""

if "annotation_supported_in_both" not in df.columns:
    df["annotation_supported_in_both"] = ""

df["placement_stable"] = (
    df["placement_stable"]
    .astype(str)
    .str.strip()
    .str.lower()
)

df["annotation_supported_in_both"] = (
    df["annotation_supported_in_both"]
    .astype(str)
    .str.strip()
    .str.lower()
)

# ------------------------------------------------------------
# Final evidence category
# ------------------------------------------------------------

def classify_final(row):

    alt = str(row.get("best_alternative_species", "")).strip()
    current = str(row.get("current_annotation", "")).strip()

    sh = row.get("sh_alrt_percent", np.nan)
    boot = row.get("ufboot_percent", np.nan)

    stable = str(row.get("placement_stable", "")).lower() == "true"
    annotation_supported = (
        str(row.get("annotation_supported_in_both", "")).lower() == "true"
    )

    classification = str(row.get("final_category", "")).strip()

    # Strongest cases already established by integrated evidence
    if (
        classification == "STRONG_ALTERNATIVE_SIGNAL"
        and stable
        and pd.notna(sh)
        and pd.notna(boot)
        and sh >= 95
        and boot >= 95
    ):
        return "HIGH_PRIORITY_TAXONOMIC_REVIEW"

    # Moderate signal, especially when reference placement is stable
    if (
        classification == "MODERATE_ALTERNATIVE_SIGNAL"
        and stable
        and pd.notna(sh)
        and pd.notna(boot)
        and sh >= 80
        and boot >= 80
    ):
        return "RETAIN_BUT_FLAG_FOR_REVIEW"

    # Weak signal should not trigger reclassification
    if classification == "WEAK_ALTERNATIVE_SIGNAL":
        return "RETAIN_CURRENT_ANNOTATION"

    # Conservative fallback
    if stable and annotation_supported:
        return "RETAIN_CURRENT_ANNOTATION"

    return "RETAIN_BUT_FLAG_FOR_REVIEW"


df["final_recommendation"] = df.apply(
    classify_final,
    axis=1
)

# ------------------------------------------------------------
# Add interpretation
# ------------------------------------------------------------

def interpretation(row):

    rec = row["final_recommendation"]

    current = row.get("current_annotation", "")
    alt = row.get("best_alternative_species", "")

    sh = row.get("sh_alrt_percent", np.nan)
    boot = row.get("ufboot_percent", np.nan)
    delta = row.get("alternative_minus_current_distance", np.nan)

    stable = str(row.get("placement_stable", "")).lower() == "true"

    if rec == "HIGH_PRIORITY_TAXONOMIC_REVIEW":
        return (
            f"Strong and sensitivity-stable placement toward {alt} "
            f"relative to current annotation {current}; "
            f"additional taxonomic evidence is warranted before any renaming."
        )

    if rec == "RETAIN_BUT_FLAG_FOR_REVIEW":
        return (
            f"Alternative placement toward {alt} is supported at a "
            f"moderate level but does not by itself establish misidentification; "
            f"retain current annotation and flag for review."
        )

    return (
        f"Alternative signal is weak or insufficient for reassignment; "
        f"retain current annotation."
    )


df["interpretation"] = df.apply(
    interpretation,
    axis=1
)

# ------------------------------------------------------------
# Reorder / clean final table
# ------------------------------------------------------------

columns = [
    "kju_id",
    "current_annotation",
    "best_alternative_species",
    "best_same_species_distance",
    "best_alternative_distance",
    "alternative_minus_current_distance",
    "sh_alrt_percent",
    "ufboot_percent",
    "main_nearest_reference",
    "sensitivity_nearest_reference",
    "main_nearest_reference_distance",
    "sensitivity_nearest_reference_distance",
    "placement_stable",
    "annotation_supported_in_both",
    "final_category",
    "final_recommendation",
    "interpretation",
]

# Keep only columns that actually exist
columns = [c for c in columns if c in df.columns]

df = df[columns].copy()

# Sort strongest cases first
priority_order = {
    "HIGH_PRIORITY_TAXONOMIC_REVIEW": 0,
    "RETAIN_BUT_FLAG_FOR_REVIEW": 1,
    "RETAIN_CURRENT_ANNOTATION": 2,
}

df["_priority"] = df["final_recommendation"].map(priority_order).fillna(9)

df = df.sort_values(
    by=["_priority", "kju_id"]
).drop(columns=["_priority"])

# ------------------------------------------------------------
# Strict validation
# ------------------------------------------------------------

expected = sorted(TARGETS)
observed = sorted(df["kju_id"].tolist())

if observed != expected:
    raise SystemExit(
        f"ERROR: Final table KJU mismatch.\n"
        f"Expected: {expected}\n"
        f"Observed: {observed}"
    )

if len(df) != 9:
    raise SystemExit(
        f"ERROR: Expected 9 rows, found {len(df)}"
    )

if df["kju_id"].duplicated().any():
    raise SystemExit("ERROR: Duplicate KJU IDs detected.")

# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

df.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

counts = (
    df["final_recommendation"]
    .value_counts()
    .to_dict()
)

summary = []

summary.append("=" * 110)
summary.append("DODDAGGUBBI - FINAL INTEGRATED 9-KJU TAXONOMIC ASSESSMENT")
summary.append("=" * 110)
summary.append("")
summary.append("Targets reviewed: 9")
summary.append("")
summary.append("EVIDENCE SOURCES")
summary.append("-" * 110)
summary.append(f"Classification: {CLASSIFICATION}")
summary.append(f"Clade support:  {SUPPORT}")
summary.append(f"Sensitivity:    {SENSITIVITY}")
summary.append(f"Evidence table: {EVIDENCE}")
summary.append("")
summary.append("FINAL RECOMMENDATION COUNTS")
summary.append("-" * 110)

for key in [
    "HIGH_PRIORITY_TAXONOMIC_REVIEW",
    "RETAIN_BUT_FLAG_FOR_REVIEW",
    "RETAIN_CURRENT_ANNOTATION",
]:
    summary.append(
        f"{key}: {counts.get(key, 0)}"
    )

summary.append("")
summary.append("FINAL TABLE")
summary.append("-" * 110)

for _, r in df.iterrows():

    sh = (
        f"{r['sh_alrt_percent']:.1f}"
        if pd.notna(r.get("sh_alrt_percent"))
        else ""
    )

    boot = (
        f"{r['ufboot_percent']:.1f}"
        if pd.notna(r.get("ufboot_percent"))
        else ""
    )

    summary.append(
        f"{r['kju_id']} | "
        f"{r['current_annotation']} | "
        f"alternative={r['best_alternative_species']} | "
        f"delta={r['alternative_minus_current_distance']:.6f} | "
        f"SH-aLRT={sh} | "
        f"UFBoot={boot} | "
        f"stable={r['placement_stable']} | "
        f"{r['final_recommendation']}"
    )

summary.append("")
summary.append("INTERPRETATION")
summary.append("-" * 110)
summary.append(
    "This table is an integrated taxonomic review layer."
)
summary.append(
    "It does NOT automatically rename or alter any KJU annotation."
)
summary.append(
    "16S rRNA evidence alone may be insufficient for definitive species reassignment."
)
summary.append(
    "High-priority cases should be confirmed using additional sequence markers, "
    "whole-genome evidence, or curated taxonomic references."
)
summary.append("")
summary.append("SOURCE DATA WAS NOT MODIFIED.")

SUMMARY.write_text(
    "\n".join(summary) + "\n",
    encoding="utf-8"
)

print("=" * 110)
print("FINAL INTEGRATED 9-KJU TAXONOMIC ASSESSMENT COMPLETE")
print("=" * 110)

print(df.to_string(index=False))

print("")
print("OUTPUT:")
print(OUT)

print("")
print("SUMMARY:")
print(SUMMARY)

print("")
print("RECOMMENDATION COUNTS:")
for k, v in counts.items():
    print(f"{k}: {v}")

print("")
print("SOURCE DATA WAS NOT MODIFIED.")
