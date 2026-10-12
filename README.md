# Doddagubbi 16S Taxonomic Analysis

Reproducible 16S rRNA sequence QC, MAFFT alignment, maximum-likelihood phylogeny, independent NCBI reference validation, sensitivity analysis, and taxonomic assessment of 109 Doddagubbi bacterial KJU isolates.

## Final dataset

- 109 Doddagubbi KJU isolates
- Main trimmed alignment: 1,345 positions
- 65 independent NCBI reference sequences
- Reference-augmented dataset: 174 taxa
- Independent-reference trimmed alignment: 1,369 positions

## Main phylogeny

IQ-TREE ModelFinder selected `GTR+F+R4` by BIC.

Support:
- 1,000 SH-aLRT replicates
- 1,000 ultrafast bootstrap replicates

## Independent-reference phylogeny

The reference-augmented analysis contained 109 KJU sequences and 65 independent NCBI references representing 13 target species.

Best-fit model by BIC:

`GTR+F+R3`

## Final taxonomic assessment

Nine KJUs were selected for detailed review.

High-priority candidates:
- KJU0017 — *Aeromonas hydrophila*; strong alternative signal toward *A. caviae*
- KJU0217 — *Bacillus cereus*; strong alternative signal toward *B. wiedmannii*

Moderate review candidates:
- KJU0007
- KJU0088
- KJU0099
- KJU0106
- KJU0221

These five *Citrobacter freundii* isolates showed a repeated alternative signal toward *Citrobacter braakii*.

Weak/unresolved candidates:
- KJU0105
- KJU0208

No master annotation was changed solely on the basis of 16S rRNA evidence.

## Reproducibility

The repository contains the curated analysis inputs, scripts, alignments, phylogenetic trees, taxonomic evidence tables, final figures, reports, and SHA-256 reproducibility manifest.

## Software

Python, Biopython, pandas, NumPy, MAFFT, IQ-TREE, and Matplotlib.

## Author

Sai Pranav

Repository:
https://github.com/saipranav-bot/Doddagubbi_16S_Taxonomic_Analysis


## Plain-English Explanation
See [EXPLAIN_PLAIN_ENGLISH.md](EXPLAIN_PLAIN_ENGLISH.md) for a non-technical summary of what was done, steps, and verification.
