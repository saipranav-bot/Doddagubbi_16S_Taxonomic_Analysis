from pathlib import Path
import pandas as pd
import re

MASTER = Path.home() / "Doddagubbi_Master"
SOURCE = Path.home() / "Doddagubbi_16S"

META = MASTER / "07_Master_Metadata"

recon_file = META / "MASTER_KJU_RECONCILIATION_001_183.csv"

if not recon_file.exists():
    raise SystemExit("ERROR: reconciliation file not found")

recon = pd.read_csv(recon_file)

missing = recon[
    recon["status"] != "DOWNLOADED"
]["kju_id"].tolist()

print("=" * 80)
print(" DODDAGGUBBI - LOCAL KJU PROVENANCE TRACE")
print("=" * 80)

print(f"\nMissing KJU IDs: {len(missing)}")
print("Searching the complete original + master project...")
print("-" * 80)

# ------------------------------------------------------------
# Files to search
# ------------------------------------------------------------

roots = [
    SOURCE,
    MASTER
]

extensions = {
    ".csv",
    ".tsv",
    ".txt",
    ".xlsx",
    ".xls",
    ".json",
    ".xml",
    ".gb",
    ".gbk",
    ".fasta",
    ".fa",
    ".fas",
    ".md",
    ".html"
}

# Avoid massive/generated directories
skip_parts = {
    ".git",
    "venv",
    "__pycache__",
    "node_modules"
}

files = []

for root in roots:

    if not root.exists():
        continue

    for p in root.rglob("*"):

        if not p.is_file():
            continue

        if p.suffix.lower() not in extensions:
            continue

        if any(part in skip_parts for part in p.parts):
            continue

        files.append(p)

files = sorted(set(files))

print(f"Files identified for search: {len(files)}")

# ------------------------------------------------------------
# Search
# ------------------------------------------------------------

rows = []

for i, kju in enumerate(missing, 1):

    print(
        f"\n[{i}/{len(missing)}] Searching {kju}"
    )

    hits = []

    pattern = re.compile(
        re.escape(kju),
        re.IGNORECASE
    )

    for file in files:

        try:

            # XLS/XLSX handled separately
            if file.suffix.lower() in {".xlsx", ".xls"}:

                try:

                    sheets = pd.read_excel(
                        file,
                        sheet_name=None,
                        dtype=str
                    )

                    for sheet, data in sheets.items():

                        data = data.fillna("")

                        mask = data.apply(
                            lambda col:
                            col.str.contains(
                                pattern,
                                regex=True,
                                na=False
                            )
                        ).any(axis=1)

                        if mask.any():

                            hits.append({
                                "file": str(file),
                                "type": "spreadsheet",
                                "sheet": str(sheet),
                                "matches": int(mask.sum())
                            })

                except Exception:
                    pass

                continue

            # Text-like files
            text = file.read_text(
                errors="ignore"
            )

            count = len(
                pattern.findall(text)
            )

            if count:

                hits.append({
                    "file": str(file),
                    "type": "text",
                    "sheet": "",
                    "matches": count
                })

        except Exception:
            continue

    print(f"  Files containing {kju}: {len(hits)}")

    for hit in hits:

        rows.append({
            "kju_id": kju,
            **hit
        })


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

out = META / "KJU_LOCAL_PROVENANCE_TRACE.csv"

result = pd.DataFrame(rows)

if len(result):

    result = result.sort_values(
        ["kju_id", "file"]
    )

result.to_csv(
    out,
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print("\n")
print("=" * 80)
print(" LOCAL PROVENANCE TRACE COMPLETE")
print("=" * 80)

print(
    f"Missing KJU IDs investigated : "
    f"{len(missing)}"
)

if len(result):

    found_ids = result["kju_id"].nunique()

    print(
        f"KJU IDs found locally       : "
        f"{found_ids}"
    )

    print(
        f"KJU IDs still absent locally: "
        f"{len(missing) - found_ids}"
    )

    print(
        f"Total matching file records : "
        f"{len(result)}"
    )

    print("\nPer-KJU file counts:")

    counts = (
        result.groupby("kju_id")
        .size()
        .sort_index()
    )

    print(counts.to_string())

else:

    print(
        "NO LOCAL FILE CONTAINS ANY OF "
        "THE 93 MISSING KJU IDs."
    )

print("\nSaved:")
print(out)

print("\n")
print("=" * 80)
print(" STEP 13 COMPLETE")
print("=" * 80)
