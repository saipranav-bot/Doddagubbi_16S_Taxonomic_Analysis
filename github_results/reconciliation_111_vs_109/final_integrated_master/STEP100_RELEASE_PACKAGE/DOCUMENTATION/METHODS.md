# Computational Methods

## Pipeline Steps (83B — 100)

### Step 83B: KJU Sequence Mapping Validation
- Validated all 3511 sequence records against 204 KJU IDs
- 2258 valid sequence-bearing rows, 1253 blank/invalid, 0 unknown KJU mappings
- Verified frozen V1 master SHA256

### Step 87: Deterministic Canonical Sequence Selection
- Priority 1: 109 KJUs with single best sequence match
- Priority 2: 2 KJUs (KJU0224, KJU0225) with accession match but no accepted length
- NO_SEQUENCE: 1 KJU (KJU0039 — shares accession with KJU0038)
- UNRESOLVED: 92 KJUs with no available sequence

### Step 88: Canonical Sequence QC
- Verified SHA256 integrity of all 111 canonical sequences
- Confirmed IUPAC nucleotide validity
- Verified source file provenance

### Step 89: Canonical FASTA Generation
- Generated 111-entry FASTA file
- Total bases: 149,522
- GC content: 54.6%
- Length range: 1232–1499 bp

### Step 90: Accession-Taxonomy Consistency Audit
- Verified 111/111 accession root consistency between canonical and accepted
- 0 discrepancies

### Step 91: Sequence Type Validation
- Verified all 111 sequences are 16S rRNA gene
- 0 anomalies detected

### Step 92: BLAST Taxonomic Validation
- Ran BLAST against 65 NCBI reference sequences
- 555 total BLAST hits
- **Known bug**: `Hit_species` column empty due to versioned vs root accession key mismatch
- **Fixed in STEP92b**: Species extracted from Subject_accession 3rd pipe field

### Step 92b: BLAST Taxonomy Species Lookup Correction
- Fixed ref_lookup to use root accession keys
- Extracted species from `Subject_accession` pipe-delimited format
- Re-classified taxonomy agreement into 4 categories:
  - SPECIES_AGREEMENT: 47 (accepted species matches top BLAST hit)
  - GENUS_AGREEMENT: 8 (same genus, >99% identity)
  - INSUFFICIENT_REFERENCE: 56 (species not in reference DB)
  - DISCORDANT: 0

### Step 93: Phylogenetic Validation
- MAFFT multiple sequence alignment
- IQ-TREE maximum likelihood phylogeny (GTR+G model)
- 111/111 canonical KJUs included in tree
- 34 genera represented

### Step 94: Final Biological Identity Adjudication
- 111 canonical adjudicated
- 1 NO_SEQUENCE (KJU0039)
- 92 UNRESOLVED
- 0 conflicts

### Step 95: V2 Integrated Master
- 204 rows, 29 columns
- All required evidence columns populated; `BLAST_taxonomy` empty (known STEP92 bug, fixed in 92b)

### Step 96: Independent QC
- 18/18 independent QC checks passed

### Step 97: Biological Uniqueness
- 111 unique biological identities
- 1 known duplicate: KJU0038/KJU0039 sharing PZ666011
- 8 cross-KJU identical-sequence pairs (canonical sequences)

### Step 98: V2 Freeze
- V2 frozen with SHA256 verification
- V1 remained unchanged

### Step 99: Independent V2 Evidence & Provenance Audit
- 18 independent checks of frozen V2
- All checks PASS
- Verified: structure, integrity, FASTA consistency, sequence hashes,
  source provenance, accession consistency, length, sequence type,
  BLAST evidence, phylogenetic evidence, biological identity,
  cross-KJU duplicates, special cases, unresolved records, data loss, V1 immutability

### Step 100: Final Release Package
- Assembled all artifacts into release package
- Created documentation (README, DATA_DICTIONARY, METHODS, KNOWN_LIMITATIONS)

## Tools Used

- **BLAST+**: Nucleotide alignment against reference database
- **MAFFT**: Multiple sequence alignment (FFT-NS-2 method)
- **IQ-TREE**: Maximum likelihood phylogeny (GTR+G model)
- **Biopython 1.88**: Sequence I/O, parsing, SHA256 hashing
- **pandas 3.0.5**: Data manipulation and analysis
- **Python 3.14**: Pipeline execution

## Reference Database

- 65 NCBI reference sequences
- 13 species represented
- Source: `NCBI_REFERENCE_METADATA.csv`
