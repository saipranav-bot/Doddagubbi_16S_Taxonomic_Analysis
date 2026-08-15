from pathlib import Path
import pandas as pd

# ============================================================
# DODDAGGUBBI
# GENERATE NCBI MASTER CHANGE PROPOSAL
# SAFE / NON-DESTRUCTIVE
# ============================================================

MASTER = Path.home() / "Doddagubbi_Master"
META = MASTER / "07_Master_Metadata"

MASTER_FILE = META / "MASTER_KJU_RECONCILIATION_001_183.csv"
RESOLUTION_FILE = META / "NCBI_RESOLUTION_TABLE_V1.csv"
KJU0039_FILE = META / "KJU0039_FINAL_RESOLUTION_V1.csv"

OUT = META / "NCBI_MASTER_CHANGE_PROPOSAL_V1.csv"
SUMMARY = META / "NCBI_MASTER_CHANGE_PROPOSAL_SUMMARY_V1.txt"

print("=" * 90)
print(" DODDAGGUBBI - NCBI MASTER CHANGE PROPOSAL")
print("=" * 90)

# ------------------------------------------------------------
# Load files
# ------------------------------------------------------------

for f in [MASTER_FILE, RESOLUTION_FILE, KJU0039_FILE]:
    if not f.exists():
        raise SystemExit(f"ERROR: Missing file:\n{f}")

master = pd.read_csv(MASTER_FILE, dtype=str).fillna("")
resolution = pd.read_csv(RESOLUTION_FILE, dtype=str).fillna("")
kju0039 = pd.read_csv(KJU0039_FILE, dtype=str).fillna("")

# ------------------------------------------------------------
# Validate master structure
# ------------------------------------------------------------

required_master = {
    "kju_id",
    "status",
    "accession",
    "organism",
    "sequence_type",
    "length",
}

missing = required_master - set(master.columns)

if missing:
    raise SystemExit(
        f"ERROR: Master file missing columns: {sorted(missing)}"
    )

if master["kju_id"].duplicated().any():
    dupes = master.loc[
        master["kju_id"].duplicated(keep=False),
        "kju_id"
    ].tolist()

    raise SystemExit(
        "ERROR: Duplicate KJU IDs in master:\n"
        + "\n".join(dupes)
    )

# ------------------------------------------------------------
# Validate resolution structure
# ------------------------------------------------------------

required_resolution = {
    "kju_id",
    "ncbi_accession",
    "ncbi_organism",
    "master_organism",
    "opencity_organism",
    "resolution",
    "classification",
    "recommended_organism",
    "recommended_action",
}

missing = required_resolution - set(resolution.columns)

if missing:
    raise SystemExit(
        f"ERROR: Resolution file missing columns: {sorted(missing)}"
    )

# ------------------------------------------------------------
# KJU0039 safety validation
# ------------------------------------------------------------

if len(kju0039) != 1 or kju0039.iloc[0]["kju_id"] != "KJU0039":
    raise SystemExit(
        "ERROR: KJU0039 resolution file is not exactly one KJU0039 record."
    )

k39 = kju0039.iloc[0]

if k39["action"] != "DO_NOT_MERGE":
    raise SystemExit(
        "ERROR: KJU0039 is not marked DO_NOT_MERGE. "
        "Refusing to generate proposal."
    )

# ------------------------------------------------------------
# Build master lookup
# ------------------------------------------------------------

master_lookup = master.set_index("kju_id", drop=False)

master_kjus = set(master["kju_id"])

# ------------------------------------------------------------
# Generate proposal records
# ------------------------------------------------------------

proposals = []

for _, r in resolution.iterrows():

    kju = r["kju_id"]

    # --------------------------------------------------------
    # KJU0039 special safety rule
    # --------------------------------------------------------

    if kju == "KJU0039":

        proposals.append({
            "kju_id": kju,
            "master_row_exists": kju in master_kjus,
            "current_master_status": (
                master_lookup.loc[kju, "status"]
                if kju in master_kjus else ""
            ),
            "current_master_accession": (
                master_lookup.loc[kju, "accession"]
                if kju in master_kjus else ""
            ),
            "current_master_organism": (
                master_lookup.loc[kju, "organism"]
                if kju in master_kjus else ""
            ),
            "ncbi_accession": r["ncbi_accession"],
            "ncbi_organism": r["ncbi_organism"],
            "classification": r["classification"],
            "resolution": r["resolution"],
            "proposed_accession": "",
            "proposed_organism": "",
            "action": "DO_NOT_MODIFY",
            "reason": (
                "KJU0039 has an accession/KJU collision: "
                "PZ666011 belongs to KJU0038 in NCBI and no independent "
                "KJU0039 NCBI record was found."
            ),
        })

        continue

    # --------------------------------------------------------
    # No corresponding master row
    # --------------------------------------------------------

    if kju not in master_kjus:

        proposals.append({
            "kju_id": kju,
            "master_row_exists": False,
            "current_master_status": "",
            "current_master_accession": "",
            "current_master_organism": "",
            "ncbi_accession": r["ncbi_accession"],
            "ncbi_organism": r["ncbi_organism"],
            "classification": r["classification"],
            "resolution": r["resolution"],
            "proposed_accession": "",
            "proposed_organism": r["recommended_organism"],
            "action": "NO_MASTER_CHANGE",
            "reason": (
                "KJU is absent from the 183-row master dataset; "
                "resolution applies to OpenCity/external reconciliation "
                "and must not create or modify a master row."
            ),
        })

        continue

    # --------------------------------------------------------
    # Existing master row
    # --------------------------------------------------------

    m = master_lookup.loc[kju]

    current_org = m["organism"]
    recommended_org = r["recommended_organism"]

    # --------------------------------------------------------
    # Only propose organism changes explicitly supported by NCBI
    # --------------------------------------------------------

    if r["classification"] == "MASTER_SUPPORTED":

        if (
            recommended_org
            and recommended_org.lower() != current_org.lower()
        ):
            action = "PROPOSE_ORGANISM_CHANGE"
            reason = (
                "NCBI validates the master organism. "
                "Recommended organism differs from current master value."
            )
        else:
            action = "NO_CHANGE"
            reason = (
                "Master organism already agrees with NCBI."
            )

        proposals.append({
            "kju_id": kju,
            "master_row_exists": True,
            "current_master_status": m["status"],
            "current_master_accession": m["accession"],
            "current_master_organism": current_org,
            "ncbi_accession": r["ncbi_accession"],
            "ncbi_organism": r["ncbi_organism"],
            "classification": r["classification"],
            "resolution": r["resolution"],
            "proposed_accession": m["accession"],
            "proposed_organism": (
                recommended_org if action == "PROPOSE_ORGANISM_CHANGE"
                else current_org
            ),
            "action": action,
            "reason": reason,
        })

    else:

        proposals.append({
            "kju_id": kju,
            "master_row_exists": True,
            "current_master_status": m["status"],
            "current_master_accession": m["accession"],
            "current_master_organism": current_org,
            "ncbi_accession": r["ncbi_accession"],
            "ncbi_organism": r["ncbi_organism"],
            "classification": r["classification"],
            "resolution": r["resolution"],
            "proposed_accession": m["accession"],
            "proposed_organism": current_org,
            "action": "NO_MASTER_CHANGE",
            "reason": (
                "Resolution does not independently justify changing "
                "the master record."
            ),
        })

# ------------------------------------------------------------
# Create proposal table
# ------------------------------------------------------------

proposal_df = pd.DataFrame(proposals)

proposal_df = proposal_df.sort_values(
    by="kju_id",
    key=lambda s: s.str.extract(r"(\d+)", expand=False).astype(int)
)

proposal_df.to_csv(
    OUT,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

lines = []

lines.append("DODDAGGUBBI - NCBI MASTER CHANGE PROPOSAL SUMMARY")
lines.append("=" * 80)
lines.append("")
lines.append(f"Master records: {len(master)}")
lines.append(f"NCBI resolution records reviewed: {len(resolution)}")
lines.append(f"Proposal records generated: {len(proposal_df)}")
lines.append("")

lines.append("ACTION COUNTS")
lines.append("-" * 80)

for action, count in proposal_df["action"].value_counts().items():
    lines.append(f"{action}: {count}")

lines.append("")
lines.append("MASTER ROW COVERAGE")
lines.append("-" * 80)

existing = proposal_df["master_row_exists"].sum()

lines.append(
    f"Resolution records with master rows: {existing}"
)

lines.append(
    f"Resolution records without master rows: "
    f"{len(proposal_df) - existing}"
)

lines.append("")
lines.append("PROPOSED MASTER ORGANISM CHANGES")
lines.append("-" * 80)

changes = proposal_df[
    proposal_df["action"] == "PROPOSE_ORGANISM_CHANGE"
]

if len(changes) == 0:
    lines.append("NONE")
else:
    for _, r in changes.iterrows():
        lines.append(
            f'{r["kju_id"]}: '
            f'{r["current_master_organism"]} -> '
            f'{r["proposed_organism"]}'
        )

lines.append("")
lines.append("KJU0039 SAFETY")
lines.append("-" * 80)
lines.append("KJU0039: DO_NOT_MODIFY")
lines.append(
    "Reason: PZ666011.1 is owned by KJU0038 in NCBI; "
    "no independent KJU0039 NCBI record was found."
)

lines.append("")
lines.append("IMPORTANT")
lines.append("-" * 80)
lines.append(
    "This file is a PROPOSAL ONLY."
)
lines.append(
    "MASTER_KJU_RECONCILIATION_001_183.csv was NOT modified."
)
lines.append(
    "No OpenCity-only KJU was added to the master."
)
lines.append(
    "KJU0039 was not merged with KJU0038."
)

SUMMARY.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8"
)

# ------------------------------------------------------------
# Final output
# ------------------------------------------------------------

print("")
print("=" * 90)
print(" PROPOSAL GENERATION COMPLETE")
print("=" * 90)

print(f"Output : {OUT}")
print(f"Summary: {SUMMARY}")

print("")
print("ACTION COUNTS:")
print(proposal_df["action"].value_counts().to_string())

print("")
print("MASTER DATA WAS NOT MODIFIED.")
