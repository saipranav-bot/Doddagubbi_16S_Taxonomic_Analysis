# Known Limitations

## 1. KJU0038/KJU0039 Shared Accession

KJU0038 (PZ666011.1) and KJU0039 (PZ666011) share the same NCBI accession
root. KJU0038 has a canonical sequence; KJU0039 has NO_SEQUENCE. Both have
biological identity "Shewanella algae". These records are documented as a
known duplicate and are NOT collapsed.

## 2. 92 Unresolved KJUs

92 KJUs have no available sequence data. Their provenance is preserved in
the V2 master, but no sequence analysis was possible.

## 3. BLAST Reference Database Limitations

The NCBI reference database contains 65 sequences covering only 13 species.
The canonical dataset spans 51 species. For 56 KJUs, the accepted species
has no reference entry, making species-level BLAST comparison impossible.

### BLAST Taxonomy Agreement Categories

| Category | Count | Description |
|----------|-------|-------------|
| SPECIES_AGREEMENT | 47 | Top BLAST hit matches accepted species |
| GENUS_AGREEMENT | 8 | Same genus, different species (>99% identity) |
| INSUFFICIENT_REFERENCE | 56 | Species not in reference DB |
| DISCORDANT | 0 | Genuine cross-genus mismatch |

The 8 GENUS_AGREEMENT cases are within-genus hits with >99% identity,
expected when the reference DB lacks the exact species.

## 4. BLAST_taxonomy Column Empty in V2

The `BLAST_taxonomy` column in V2 is empty due to a bug in STEP92 where
the species lookup used versioned accession keys but the code normalized
to root accessions. This was fixed in STEP92b. Corrected data is available
in `STEP92b_TAXONOMY_COMPARISON_CORRECTED.csv`.

## 5. Cross-KJU Identical Sequences

8 pairs of canonical KJUs share identical sequences. Biological distinctness
is assumed from metadata (different accessions, different source records).
These pairs are documented in `STEP97_CROSS_KJU_RELATIONSHIPS.csv`.

## 6. 16S Only

All canonical sequences are 16S rRNA gene sequences. No ITS sequences
are included in this dataset.

## 7. Local BLAST Only

BLAST validation was performed against a local database of 65 NCBI
reference sequences, not the full NCBI nt database. Results are
limited to the reference set.

## 8. Priority 2 KJUs

KJU0224 and KJU0225 were selected by Priority 2 (accession match with
no accepted length in master). Their canonical sequences are validated
but the selection criteria differ from Priority 1 KJUs.
