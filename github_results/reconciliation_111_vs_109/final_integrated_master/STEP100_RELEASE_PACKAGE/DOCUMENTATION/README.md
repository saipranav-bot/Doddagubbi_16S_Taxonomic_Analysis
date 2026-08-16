# Doddagubbi 16S/ITS KJU Reconciliation — Final Release

## Overview

This dataset contains 204 KJU (Known JU sequence) records reconciled from
NCBI Sequence Read Archive data. Of these:

- **111** have canonical 16S rRNA gene sequences
- **1** has no sequence (KJU0039 — shares accession PZ666011 with KJU0038)
- **92** are unresolved (no sequence available)

## Files

### DATA/

| File | Description |
|------|-------------|
| `FINAL_204_KJU_ADJUDICATED_MASTER_V2.csv` | V2 frozen master — 204 rows, 29 columns |
| `FINAL_BIOLOGICALLY_UNIQUE_KJU_MASTER_V2.csv` | 111 unique biological identities |
| `FINAL_CANONICAL_KJU_16S_V2.fasta` | 111 canonical 16S sequences |

### PROVENANCE/

| File | Description |
|------|-------------|
| `STEP92_BLAST_VALIDATION/` | BLAST results, taxonomy comparison, STEP92b fix |
| `STEP93_PHYLOGENETIC_VALIDATION/` | MAFFT alignment, IQ-TREE, clustering audit |
| `STEP99_INDEPENDENT_V2_PROVENANCE_AUDIT/` | 18-check independent audit (PASS) |
| `STEP83_TRUE_SEQUENCE_RECORDS.csv` | Raw sequence inventory (3511 records) |

### QC_REPORTS/

Quality control reports from each processing step.

### MANIFESTS/

| File | Description |
|------|-------------|
| `FINAL_KJU_MASTER_FREEZE_MANIFEST_V2.txt` | V2 freeze manifest with SHA256 checksums |

### DOCUMENTATION/

| File | Description |
|------|-------------|
| `DATA_DICTIONARY.csv` | Column definitions for V2 master |
| `METHODS.md` | Computational methods |
| `KNOWN_LIMITATIONS.md` | Known limitations and caveats |

## Key Statistics

| Metric | Value |
|--------|-------|
| Total KJUs | 204 |
| Canonical sequences | 111 |
| NO_SEQUENCE | 1 |
| UNRESOLVED | 92 |
| Unique biological identities | 111 |
| BLAST species agreement | 47/111 (8 genus-agreement, 56 insufficient reference, 0 discordant) |
| Phylogenetic inclusion | 111/111 |
| Independent QC checks | 18/18 PASS |

## Frozen Versions

| Version | SHA256 | Status |
|---------|--------|--------|
| V1 | `7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8` | Frozen, immutable |
| V2 | `65185f7545f73fcf36ff5062bb85305affcae76458f8049130f052a811b79bee` | Frozen, current |

## How to Cite

If using this dataset, please cite the processing pipeline and the original
NCBI SRA data sources.

## Contact

Doddagubbi 16S/ITS KJU Reconciliation Project
