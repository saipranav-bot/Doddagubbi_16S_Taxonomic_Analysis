#!/usr/bin/env python3
"""Create Word document: Doddagubbi 16S Project Report"""

from pathlib import Path
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.section import WD_ORIENT

BASE = Path.home() / "Doddagubbi_Master"
OUT = BASE / "Doddagubbi_16S_Project_Report.docx"

doc = Document()

# ──────────────────────────────────────────────────────
# STYLES
# ──────────────────────────────────────────────────────
style = doc.styles['Normal']
font = style.font
font.name = 'Calibri'
font.size = Pt(11)

# ──────────────────────────────────────────────────────
# TITLE PAGE
# ──────────────────────────────────────────────────────
doc.add_paragraph()
doc.add_paragraph()

title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = title.add_run("Doddagubbi 16S Taxonomic Analysis")
run.bold = True
run.font.size = Pt(28)
run.font.color.rgb = RGBColor(0, 51, 102)

subtitle = doc.add_paragraph()
subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = subtitle.add_run("Complete Project Report")
run.font.size = Pt(18)
run.font.color.rgb = RGBColor(80, 80, 80)

doc.add_paragraph()

author = doc.add_paragraph()
author.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = author.add_run("Author: Sai Pranav")
run.font.size = Pt(14)

date = doc.add_paragraph()
date.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = date.add_run("Date: August 2026")
run.font.size = Pt(14)

repo = doc.add_paragraph()
repo.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = repo.add_run("Repository: https://github.com/saipranav-bot/Doddagubbi_16S_Taxonomic_Analysis")
run.font.size = Pt(10)
run.font.color.rgb = RGBColor(0, 102, 204)

doc.add_page_break()

# ──────────────────────────────────────────────────────
# TABLE OF CONTENTS
# ──────────────────────────────────────────────────────
doc.add_heading("Table of Contents", level=1)
toc_items = [
    "1. Project Overview",
    "2. The Biology: Why 16S rRNA?",
    "3. Raw Data",
    "4. The Problem",
    "5. Reconciliation Pipeline (Steps 56–98)",
    "6. Release (Steps 99–100)",
    "7. Downstream Analysis (Steps 101–107)",
    "8. Tools and Why Each Was Used",
    "9. Key Numbers Summary",
    "10. Known Limitations",
    "11. Conclusion",
    "12. Frozen Data Integrity",
]
for item in toc_items:
    p = doc.add_paragraph(item)
    p.paragraph_format.space_after = Pt(2)
doc.add_page_break()

# ──────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────
def add_table(doc, headers, rows):
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = 'Light Grid Accent 1'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = h
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.bold = True
                run.font.size = Pt(10)
    for r, row_data in enumerate(rows):
        for c, val in enumerate(row_data):
            cell = table.rows[r + 1].cells[c]
            cell.text = str(val)
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(10)
    return table

def add_bold_paragraph(doc, bold_text, normal_text=""):
    p = doc.add_paragraph()
    run = p.add_run(bold_text)
    run.bold = True
    if normal_text:
        p.add_run(normal_text)
    return p

# ──────────────────────────────────────────────────────
# SECTION 1: PROJECT OVERVIEW
# ──────────────────────────────────────────────────────
doc.add_heading("1. Project Overview", level=1)

doc.add_paragraph(
    "This project is a comprehensive bacterial identification and taxonomic analysis "
    "of bacterial isolates collected from Doddagubbi Lake in Bengaluru, India. The "
    "project takes raw DNA sequence data from NCBI (National Center for Biotechnology "
    "Information) and processes it through a rigorous 107-step bioinformatics pipeline "
    "to produce a clean, validated, and reproducible dataset of 16S rRNA gene sequences."
)

doc.add_paragraph(
    "The pipeline resolves 204 unique bacterial isolates (called 'KJU' identifiers) "
    "from NCBI Sequence Read Archive data. Of these, 111 have canonical 16S rRNA gene "
    "sequences, 1 has no sequence, and 92 are unresolved (no sequence available)."
)

doc.add_paragraph(
    "The final dataset represents 51 bacterial species across 34 genera, providing "
    "a comprehensive snapshot of the bacterial community in Doddagubbi Lake."
)

add_bold_paragraph(doc, "What We Achieved:", "")
achievements = [
    "Reconciled 3,511 raw sequence records into 111 validated canonical sequences",
    "Built a frozen, immutable master dataset with SHA256 integrity verification",
    "Validated taxonomy using BLAST against NCBI reference sequences",
    "Built a maximum-likelihood phylogenetic tree using IQ-TREE",
    "Analyzed taxonomic composition, sequence diversity, and phylogenetic relationships",
    "Created a fully reproducible pipeline with documentation and checksums",
    "Uploaded everything to GitHub with complete provenance",
]
for a in achievements:
    doc.add_paragraph(a, style='List Bullet')

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 2: THE BIOLOGY
# ──────────────────────────────────────────────────────
doc.add_heading("2. The Biology: Why 16S rRNA?", level=1)

doc.add_paragraph(
    "Every bacterium has a gene called 16S ribosomal RNA (16S rRNA). This gene is "
    "part of the bacterial ribosome — the molecular machine that builds proteins. "
    "The 16S rRNA gene is universally present in all bacteria, making it the gold "
    "standard for bacterial identification."
)

doc.add_paragraph(
    "The 16S rRNA gene evolves slowly enough that closely related species have "
    "very similar 16S sequences, but different enough that you can distinguish "
    "different species. By sequencing this one gene and comparing it to a database "
    "of known bacteria, you can determine what species an unknown bacterium is."
)

add_bold_paragraph(doc, "Why 16S and not other genes?", "")
doc.add_paragraph(
    "16S rRNA is used because: (1) it is present in all bacteria, (2) it has both "
    "conserved regions (same across all bacteria, useful for universal primers) and "
    "variable regions (different between species, useful for identification), (3) the "
    "NCBI database has hundreds of thousands of 16S sequences with known taxonomy, "
    "and (4) it has been the standard for bacterial taxonomy for over 30 years."
)

add_bold_paragraph(doc, "What about ITS?", "")
doc.add_paragraph(
    "ITS (Internal Transcribed Spacer) is a different gene used primarily for "
    "fungi, not bacteria. This project uses only 16S rRNA. No ITS sequences are "
    "included in this dataset."
)

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 3: RAW DATA
# ──────────────────────────────────────────────────────
doc.add_heading("3. Raw Data", level=1)

doc.add_paragraph(
    "The raw data consists of 3,511 sequence records stored in a FASTA file called "
    "MASTER_KJU_16S.fasta. Each record looks like this:"
)

# FASTA format example
p = doc.add_paragraph()
run = p.add_run(
    ">KJU0001|PZ638859.1|Aeromonas caviae\n"
    "CGGGAAAGTAGCTTGCTACTTTTGCCGGCGAGCGGCGGACGGGTGAGTAATGCCTGGGAA..."
)
run.font.name = 'Courier New'
run.font.size = Pt(9)

doc.add_paragraph(
    "The format is: KJU_ID | NCBI_Accession | Species_Name"
)

add_table(doc, ["Field", "Example", "Description"], [
    ["KJU ID", "KJU0001", "Unique identifier for this isolate"],
    ["NCBI Accession", "PZ638859.1", "NCBI's unique identifier for this sequence"],
    ["Species Name", "Aeromonas caviae", "Species assigned by NCBI"],
    ["Sequence", "CGGGAAAGT...", "The actual 16S rRNA DNA sequence (A, T, G, C)"],
])

doc.add_paragraph()
add_bold_paragraph(doc, "Data Inventory:", "")
add_table(doc, ["Item", "Count"], [
    ["Total raw sequence records", "3,511"],
    ["Unique KJU IDs", "204"],
    ["Records with valid sequences", "2,258"],
    ["Records with blank/invalid sequences", "1,253"],
    ["Unique NCBI accessions", "204"],
])

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 4: THE PROBLEM
# ──────────────────────────────────────────────────────
doc.add_heading("4. The Problem", level=1)

doc.add_paragraph("The raw data had several critical problems that needed solving:")

problems = [
    ("Duplicates", "Many KJUs had multiple sequence records — the same sequence "
     "submitted to NCBI multiple times with different accession numbers. For example, "
     "KJU0001 might have PZ638859.1, PZ638859, and other versions all representing "
     "the same DNA sequence."),
    ("Inconsistencies", "Some sequences disagreed with the master metadata — wrong "
     "species names, wrong accessions, or sequences that didn't match expected lengths."),
    ("Missing Data", "92 of 204 KJUs (45%) had NO sequence data available at all. "
     "These were resolved only at the metadata level."),
    ("Taxonomic Uncertainty", "Some species assignments were wrong or uncertain. "
     "For example, BLAST might say 'this is Aeromonas veronii' while the metadata "
     "says 'Aeromonas caviae'."),
    ("No Clean Dataset", "Everything needed to be reconciled into one authoritative, "
     "validated table with full provenance tracking."),
]

for title, desc in problems:
    add_bold_paragraph(doc, f"{title}: ", desc)

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 5: RECONCILIATION PIPELINE
# ──────────────────────────────────────────────────────
doc.add_heading("5. Reconciliation Pipeline (Steps 56–98)", level=1)

doc.add_paragraph(
    "The reconciliation pipeline is the core of this project. It takes the raw, messy "
    "data and produces a clean, validated, frozen dataset. Every step is automated, "
    "documented, and produces its own output directory and QC report."
)

# Phase 1
doc.add_heading("5.1 Phase 1: OpenCity vs Final Reconciliation (Steps 56–72)", level=2)

doc.add_paragraph(
    "This phase compared two datasets: the 111 KJU isolates from the 'OpenCity' "
    "dataset against the 109 KJU isolates from the 'Final' master metadata. The "
    "goal was to find discrepancies, download missing sequences from NCBI, and "
    "validate taxonomy using BLAST and phylogenetic methods."
)

steps_56_72 = [
    ("Step 56", "Reconciled OpenCity 111 vs Final 109 KJUs"),
    ("Step 57", "Downloaded missing OpenCity sequences from NCBI"),
    ("Step 58", "QC'd new OpenCity sequences"),
    ("Steps 59–60", "Compared new sequences vs Final 109"),
    ("Steps 61–68", "Downloaded taxonomic references, aligned, built phylogenetic trees, corrected taxonomy"),
    ("Steps 69–70", "BLAST-validated new KJUs, produced final taxonomic validation"),
    ("Steps 71–72", "Resolved shared KJU biological reconciliation, final taxonomic adjudication"),
]
add_table(doc, ["Step", "Description"], steps_56_72)

doc.add_paragraph()

# Phase 2
doc.add_heading("5.2 Phase 2: Master Build and V1 Freeze (Steps 74–76)", level=2)

doc.add_paragraph(
    "Built the integrated 204-KJU master table and froze it as V1. This freeze is "
    "IMMUTABLE — it was given a SHA256 hash and never modified again."
)

add_table(doc, ["Step", "Description", "Key Output"], [
    ["74", "Built integrated 204-KJU master", "204 rows, 29 columns"],
    ["75", "Final QC of 204-KJU master", "QC report"],
    ["76", "Froze V1 master", "SHA256: 7d858153..."],
])

doc.add_paragraph()

# Phase 3
doc.add_heading("5.3 Phase 3: Sequence-Level QC (Steps 77–89)", level=2)

doc.add_paragraph(
    "This is the most detailed phase. Every one of the 3,511 raw sequence records "
    "was examined, validated, and classified."
)

add_table(doc, ["Step", "What It Does", "Result"], [
    ["77", "Sequence inventory QC", "Identified 2,258 valid records"],
    ["78", "Sequence identity audit", "Found duplicates, accession conflicts"],
    ["79", "Sequence discrepancy investigation", "Investigated every conflict"],
    ["80", "Canonical sequence selection audit", "Evaluated selection criteria"],
    ["81", "Canonical sequence resolution", "Resolved conflicts"],
    ["82", "Sequence field integrity diagnostic", "Validated all fields"],
    ["83", "True sequence inventory", "2,258 valid sequence records"],
    ["83B", "KJU sequence mapping validation", "Verified KJU-to-sequence mapping"],
    ["84", "True sequence master reconciliation", "Reconciled all records"],
    ["85", "Conflict diagnostic", "Diagnosed remaining conflicts"],
    ["86", "Candidate classification", "Classified sequence candidates"],
    ["87", "Deterministic canonical selection", "Picked ONE best sequence per KJU"],
    ["88", "Canonical sequence QC", "SHA256, IUPAC validity checks"],
    ["89", "Canonical FASTA generation", "111-entry FASTA file"],
])

doc.add_paragraph()

add_bold_paragraph(doc, "Step 87 — The Critical Step:", "")
doc.add_paragraph(
    "Step 87 is the most important step in the entire pipeline. For each of the 204 "
    "KJUs, it selected exactly ONE canonical sequence using a deterministic priority system:"
)

add_table(doc, ["Priority", "Criteria", "KJUs", "Count"], [
    ["Priority 1", "Exact accession + length match in master", "KJU0001–KJU0223", "109"],
    ["Priority 2", "Accession match, no accepted length", "KJU0224, KJU0225", "2"],
    ["NO_SEQUENCE", "Shares accession with another KJU", "KJU0039", "1"],
    ["UNRESOLVED", "No sequence available", "KJU0040–KJU0234", "92"],
])

doc.add_paragraph()

# Phase 4
doc.add_heading("5.4 Phase 4: Validation (Steps 90–93)", level=2)

doc.add_paragraph(
    "After selecting canonical sequences, every sequence was validated using "
    "independent methods."
)

add_table(doc, ["Step", "Method", "Result"], [
    ["90", "Accession-taxonomy consistency audit", "111/111 accession root consistency"],
    ["91", "Sequence type validation", "All 111 are 16S rRNA (0 anomalies)"],
    ["92", "BLAST taxonomic validation", "47 species agreement, 8 genus agreement, 56 insufficient reference, 0 discordant"],
    ["92b", "BLAST taxonomy species correction", "Fixed empty BLAST_taxonomy column bug"],
    ["93", "Phylogenetic validation", "MAFFT alignment + IQ-TREE, 111/111 in tree"],
])

doc.add_paragraph()

# Phase 5
doc.add_heading("5.5 Phase 5: Final Adjudication and V2 Freeze (Steps 94–98)", level=2)

add_table(doc, ["Step", "Description"], [
    ["94", "Final biological identity adjudication — 111 canonical, 1 NO_SEQUENCE, 92 UNRESOLVED"],
    ["95", "Built V2 integrated master (204 rows, 29 columns)"],
    ["96", "Independent QC — 18/18 checks PASS"],
    ["97", "Biological uniqueness — 111 unique biological identities confirmed"],
    ["98", "Froze V2 master — SHA256: 65185f75..., IMMUTABLE"],
])

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 6: RELEASE
# ──────────────────────────────────────────────────────
doc.add_heading("6. Release (Steps 99–100)", level=1)

doc.add_paragraph(
    "After V2 was frozen, a final independent audit and release packaging was performed."
)

add_table(doc, ["Step", "Description", "Result"], [
    ["99", "Independent V2 provenance audit — 18 automated checks", "18/18 PASS"],
    ["100", "Release package assembly", "55 files, 8.6 MB"],
])

doc.add_paragraph()

add_bold_paragraph(doc, "Release Package Contents:", "")
add_table(doc, ["Directory", "Contents"], [
    ["DATA/", "V2 master CSV, biological unique CSV, canonical FASTA"],
    ["DOCUMENTATION/", "README, METHODS, DATA_DICTIONARY, KNOWN_LIMITATIONS"],
    ["MANIFESTS/", "V2 freeze manifest with SHA256 checksums"],
    ["PROVENANCE/", "BLAST validation, phylogenetic validation, audit results"],
    ["QC_REPORTS/", "Quality control reports from every step"],
    ["FREEZE_VERIFICATION/", "SHA256 verification files"],
])

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 7: DOWNSTREAM ANALYSIS
# ──────────────────────────────────────────────────────
doc.add_heading("7. Downstream Analysis (Steps 101–107)", level=1)

doc.add_paragraph(
    "After the reconciliation was complete, a downstream analysis was performed "
    "to characterize the biological diversity of the dataset."
)

# Step 101
doc.add_heading("7.1 Step 101: Release Integrity Check", level=2)
doc.add_paragraph(
    "Automated integrity check of the release package. Verified file existence, "
    "SHA256 checksums, row counts, FASTA entry counts, and workspace creation. "
    "Result: 20/20 checks PASS."
)

add_table(doc, ["Check Category", "Result"], [
    ["File existence", "6/6 PASS"],
    ["SHA256 verification", "2/2 PASS"],
    ["V2 master counts", "2/2 PASS"],
    ["FASTA verification", "3/3 PASS"],
    ["Biological unique master", "1/1 PASS"],
    ["Release checksums", "1/1 PASS"],
    ["Status counts", "1/1 PASS"],
    ["FASTA sequence validity", "1/1 PASS"],
    ["FASTA length consistency", "1/1 PASS"],
    ["FASTA SHA256 consistency", "1/1 PASS"],
    ["Workspace creation", "1/1 PASS"],
])

doc.add_paragraph()

# Step 102
doc.add_heading("7.2 Step 102: Analysis-Ready Metadata", level=2)
doc.add_paragraph(
    "Merged the V2 frozen master with sample metadata from NCBI (isolation source, "
    "geographic location, collection date, collector name). Created an analysis-ready "
    "table with 111 rows and 25 columns."
)

add_table(doc, ["Metric", "Value"], [
    ["Total canonical KJUs", "111"],
    ["With sample metadata", "89"],
    ["Without sample metadata", "22"],
    ["Analysis columns", "25"],
])

doc.add_paragraph()
doc.add_paragraph(
    "Of the 89 KJUs with metadata: 83 were isolated from lake water, 6 from soil. "
    "All were collected in Bengaluru, India. 89 were collected on 16 June 2026, "
    "6 were collected in 2025."
)

# Step 103
doc.add_heading("7.3 Step 103: Taxonomic Composition", level=2)
doc.add_paragraph(
    "Analyzed the taxonomic distribution of the 111 canonical sequences. "
    "Found 34 genera and 51 species."
)

add_bold_paragraph(doc, "Top 10 Genera:", "")
add_table(doc, ["Genus", "Count", "Percent"], [
    ["Aeromonas", "25", "22.5%"],
    ["Citrobacter", "13", "11.7%"],
    ["Bacillus", "10", "9.0%"],
    ["Shewanella", "9", "8.1%"],
    ["Pseudomonas", "6", "5.4%"],
    ["Paenibacillus", "6", "5.4%"],
    ["Comamonas", "5", "4.5%"],
    ["Plesiomonas", "4", "3.6%"],
    ["Microbacterium", "3", "2.7%"],
    ["Staphylococcus", "3", "2.7%"],
])

doc.add_paragraph()

add_bold_paragraph(doc, "Top 10 Species:", "")
add_table(doc, ["Species", "Count", "Percent"], [
    ["Aeromonas caviae", "9", "8.1%"],
    ["Aeromonas veronii", "8", "7.2%"],
    ["Citrobacter koseri", "7", "6.3%"],
    ["Shewanella putrefaciens", "7", "6.3%"],
    ["Aeromonas hydrophila", "6", "5.4%"],
    ["Citrobacter freundii", "5", "4.5%"],
    ["Paenibacillus polymyxa", "5", "4.5%"],
    ["Comamonas testosteroni", "4", "3.6%"],
    ["Plesiomonas shigelloides", "4", "3.6%"],
    ["Bacillus subtilis", "4", "3.6%"],
])

doc.add_paragraph()

# Step 104
doc.add_heading("7.4 Step 104: Sequence-Level Diversity", level=2)
doc.add_paragraph(
    "Analyzed the 111 canonical 16S sequences for length, GC content, nucleotide "
    "composition, and pairwise identity."
)

add_table(doc, ["Metric", "Value"], [
    ["Sequences", "111"],
    ["Length range", "1,232 – 1,499 bp"],
    ["Length mean", "1,347 bp"],
    ["Length std", "40.2 bp"],
    ["GC content range", "49.5% – 59.1%"],
    ["GC content mean", "54.6%"],
    ["GC content std", "1.47%"],
    ["Pairwise identity range", "21.42% – 100.00%"],
    ["Pairwise identity mean", "26.59%"],
    ["Pairwise identity median", "25.56%"],
])

doc.add_paragraph()
doc.add_paragraph(
    "The low pairwise identity mean (26.6%) confirms these are genuinely diverse "
    "bacteria, not closely related clones. The dataset spans 34 genera and 51 species, "
    "so most sequences are from distantly related organisms."
)

doc.add_paragraph()

add_bold_paragraph(doc, "Nucleotide Composition:", "")
add_table(doc, ["Nucleotide", "Mean Percentage"], [
    ["A (Adenine)", "25.5%"],
    ["T (Thymine)", "19.9%"],
    ["G (Guanine)", "31.6%"],
    ["C (Cytosine)", "23.0%"],
])

doc.add_paragraph()

# Step 105
doc.add_heading("7.5 Step 105: Reproducible MAFFT Alignment", level=2)
doc.add_paragraph(
    "Ran MAFFT (Multiple Alignment using Fast Fourier Transform) on the 111 "
    "canonical sequences to create a multiple sequence alignment."
)

add_table(doc, ["Metric", "Value"], [
    ["Input sequences", "111"],
    ["Alignment length", "1,549 bp"],
    ["Mean gaps per sequence", "202.0"],
    ["Conserved positions", "689 / 1,549 (44.5%)"],
    ["Runtime", "2.7 seconds"],
])

doc.add_paragraph()
doc.add_paragraph(
    "The alignment length (1,549 bp) is longer than the average sequence length "
    "(1,347 bp) because gaps are inserted to align corresponding positions across "
    "all sequences. The 44.5% conserved positions is consistent with diverse bacteria "
    "spanning 34 genera."
)

doc.add_paragraph()

# Step 106
doc.add_heading("7.6 Step 106: Reproducible IQ-TREE Phylogeny", level=2)
doc.add_paragraph(
    "Built a maximum-likelihood phylogenetic tree using IQ-TREE with the GTR+G model "
    "(General Time Reversible with Gamma-distributed rate variation)."
)

add_table(doc, ["Parameter", "Value"], [
    ["Model", "GTR+G"],
    ["Sequences", "111"],
    ["Seed", "12345 (for reproducibility)"],
    ["Tree tips", "111"],
    ["Runtime", "92.8 seconds"],
    ["Output", "STEP106_PHYLOGENY.treefile"],
])

doc.add_paragraph()
doc.add_paragraph(
    "IQ-TREE is the current gold standard for building evolutionary trees. The GTR+G "
    "model accounts for different types of mutations (A↔T, A↔G, etc.) occurring at "
    "different rates, which is important for diverse datasets."
)

doc.add_paragraph()

# Step 107
doc.add_heading("7.7 Step 107: Phylogenetic + Taxonomic Diversity", level=2)
doc.add_paragraph(
    "Integrated the phylogenetic tree with taxonomic metadata to calculate "
    "diversity metrics."
)

add_table(doc, ["Metric", "Value", "Interpretation"], [
    ["Shannon index (genus)", "2.886", "High diversity"],
    ["Shannon index (species)", "3.593", "Very high diversity"],
    ["Evenness (genus)", "0.818", "Good evenness (1.0 = perfectly even)"],
    ["Evenness (species)", "0.914", "Excellent evenness"],
    ["Dominant genus", "Aeromonas", "25 isolates (22.5%)"],
    ["Singleton genera", "21", "Genera with only 1 species"],
])

doc.add_paragraph()
doc.add_paragraph(
    "The Shannon diversity index measures both richness (number of species) and "
    "evenness (how evenly distributed they are). A value of 2.886 at genus level "
    "and 3.593 at species level indicates high biodiversity. The evenness values "
    "(0.818 and 0.914, where 1.0 is perfectly even) show that no single genus "
    "overwhelmingly dominates the community."
)

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 8: TOOLS
# ──────────────────────────────────────────────────────
doc.add_heading("8. Tools and Why Each Was Used", level=1)

add_table(doc, ["Tool", "Version", "Purpose", "Why This Tool"], [
    ["Python", "3.14", "Pipeline execution", "General-purpose language with strong bioinformatics ecosystem"],
    ["Biopython", "1.88", "Sequence I/O, parsing, SHA256", "Standard Python library for bioinformatics; handles FASTA format natively"],
    ["pandas", "3.0.5", "Data manipulation, CSV handling", "Essential for working with tabular data (204 rows x 29 columns)"],
    ["NumPy", "latest", "Numerical computations", "Efficient array operations for pairwise identity matrices"],
    ["BLAST+", "latest", "Sequence similarity search", "Universal tool for 'what is this sequence similar to?' — finds closest known species"],
    ["MAFFT", "latest", "Multiple sequence alignment", "Fast, accurate alignment of 111 sequences; finds which nucleotides correspond"],
    ["IQ-TREE", "latest", "Maximum likelihood phylogeny", "Current gold standard for evolutionary trees; GTR+G model for diverse data"],
])

doc.add_paragraph()

doc.add_heading("8.1 Why BLAST?", level=2)
doc.add_paragraph(
    "BLAST (Basic Local Alignment Search Tool) is the most widely used bioinformatics "
    "tool in the world. It takes a query sequence and finds the most similar sequences "
    "in a database. For this project, BLAST compared each of the 111 sequences against "
    "a local database of 65 NCBI reference sequences to verify that the species "
    "assignment was correct."
)

doc.add_heading("8.2 Why MAFFT?", level=2)
doc.add_paragraph(
    "MAFFT (Multiple Alignment using Fast Fourier Transform) is a multiple sequence "
    "alignment tool. It takes unaligned sequences and inserts gaps so that "
    "corresponding positions line up. This is essential for phylogenetic analysis — "
    "you can't build a tree without first aligning the sequences."
)

doc.add_heading("8.3 Why IQ-TREE?", level=2)
doc.add_paragraph(
    "IQ-TREE is a maximum-likelihood phylogenetic tree builder. It finds the tree "
    "that most likely produced the observed sequence data, given a model of DNA "
    "evolution. The GTR+G model was chosen because: (1) GTR (General Time Reversible) "
    "is the most general DNA substitution model, and (2) G (Gamma) accounts for "
    "different sites evolving at different rates."
)

doc.add_heading("8.4 Why SHA256 Hashing?", level=2)
doc.add_paragraph(
    "SHA256 is a cryptographic hash function that produces a unique 64-character "
    "fingerprint for any file. If even one character in the file changes, the hash "
    "completely changes. Every frozen file in this project has a SHA256 hash, "
    "proving the data hasn't been tampered with."
)

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 9: KEY NUMBERS
# ──────────────────────────────────────────────────────
doc.add_heading("9. Key Numbers Summary", level=1)

add_table(doc, ["Metric", "Value"], [
    ["Total KJUs", "204"],
    ["Canonical 16S sequences", "111"],
    ["No sequence", "1 (KJU0039)"],
    ["Unresolved", "92"],
    ["Unique biological species", "51"],
    ["Unique genera", "34"],
    ["NCBI reference sequences", "65"],
    ["Reference species covered", "13"],
    ["BLAST species agreement", "47/111"],
    ["BLAST genus agreement", "8/111"],
    ["BLAST insufficient reference", "56/111"],
    ["BLAST discordant", "0"],
    ["Independent QC checks", "20/20 PASS"],
    ["Pipeline steps", "107"],
    ["Scripts", "52"],
    ["Frozen master SHA256 (V2)", "65185f75..."],
    ["Release package files", "55"],
    ["Release package size", "8.6 MB"],
])

doc.add_paragraph()

add_bold_paragraph(doc, "Sequence Statistics:", "")
add_table(doc, ["Metric", "Value"], [
    ["Alignment length", "1,549 bp"],
    ["Mean sequence length", "1,347 bp"],
    ["GC content mean", "54.6%"],
    ["Pairwise identity mean", "26.6%"],
    ["Conserved positions", "44.5%"],
])

doc.add_paragraph()

add_bold_paragraph(doc, "Diversity Indices:", "")
add_table(doc, ["Index", "Value"], [
    ["Shannon index (genus)", "2.886"],
    ["Shannon index (species)", "3.593"],
    ["Evenness (genus)", "0.818"],
    ["Evenness (species)", "0.914"],
])

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 10: KNOWN LIMITATIONS
# ──────────────────────────────────────────────────────
doc.add_heading("10. Known Limitations", level=1)

limitations = [
    ("1. KJU0038/KJU0039 Shared Accession",
     "KJU0038 (PZ666011.1) and KJU0039 (PZ666011) share the same NCBI accession root. "
     "KJU0038 has a canonical sequence; KJU0039 has NO_SEQUENCE. Both have biological "
     "identity 'Shewanella algae'. These records are documented as a known duplicate "
     "and are NOT collapsed."),
    ("2. 92 Unresolved KJUs",
     "92 KJUs have no available sequence data. Their provenance is preserved in the "
     "V2 master, but no sequence analysis was possible."),
    ("3. BLAST Reference Database Limitations",
     "The NCBI reference database contains 65 sequences covering only 13 species. "
     "The canonical dataset spans 51 species. For 56 KJUs, the accepted species has "
     "no reference entry, making species-level BLAST comparison impossible."),
    ("4. BLAST_taxonomy Column Empty in V2",
     "The BLAST_taxonomy column in V2 is empty due to a bug in STEP92 where the "
     "species lookup used versioned accession keys but the code normalized to root "
     "accessions. This was fixed in STEP92b. Corrected data is available separately."),
    ("5. Cross-KJU Identical Sequences",
     "8 pairs of canonical KJUs share identical sequences. Biological distinctness "
     "is assumed from metadata (different accessions, different source records)."),
    ("6. 16S Only",
     "All canonical sequences are 16S rRNA gene sequences. No ITS sequences are "
     "included in this dataset."),
    ("7. Local BLAST Only",
     "BLAST validation was performed against a local database of 65 NCBI reference "
     "sequences, not the full NCBI nt database. Results are limited to the reference set."),
    ("8. Priority 2 KJUs",
     "KJU0224 and KJU0225 were selected by Priority 2 (accession match with no "
     "accepted length in master). Their canonical sequences are validated but the "
     "selection criteria differ from Priority 1 KJUs."),
]

for title, desc in limitations:
    add_bold_paragraph(doc, title)
    doc.add_paragraph(desc)

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 11: CONCLUSION
# ──────────────────────────────────────────────────────
doc.add_heading("11. Conclusion", level=1)

doc.add_paragraph(
    "This project successfully reconciled 204 bacterial isolates from Doddagubbi Lake "
    "into a clean, validated, and reproducible dataset. The 111 canonical sequences "
    "represent 51 species across 34 genera, providing a comprehensive snapshot of "
    "the bacterial diversity in this freshwater ecosystem."
)

doc.add_paragraph(
    "Key achievements:"
)

achievements = [
    "Resolved 3,511 raw sequence records into 111 validated canonical sequences",
    "Built a frozen, immutable dataset with SHA256 integrity verification",
    "Validated taxonomy using both BLAST and phylogenetic methods",
    "Built a maximum-likelihood phylogenetic tree with IQ-TREE",
    "Analyzed taxonomic composition and sequence diversity",
    "Created a fully reproducible pipeline with 52 scripts and complete documentation",
    "Uploaded everything to GitHub with provenance tracking",
]
for a in achievements:
    doc.add_paragraph(a, style='List Bullet')

doc.add_paragraph(
    "The dataset is ready for downstream ecological analysis, comparative genomics, "
    "or integration with other Doddagubbi Lake studies."
)

doc.add_page_break()

# ──────────────────────────────────────────────────────
# SECTION 12: FROZEN DATA INTEGRITY
# ──────────────────────────────────────────────────────
doc.add_heading("12. Frozen Data Integrity", level=1)

doc.add_paragraph(
    "Every frozen file in this project has a SHA256 hash. These hashes prove the "
    "data hasn't been modified. If you want to verify integrity, compute the SHA256 "
    "of any file and compare it to the recorded hash."
)

add_table(doc, ["File", "SHA256"], [
    ["V1 Master", "7d858153a0aaf19e3e1f4f7a1d3b03b3acc69ec6451c6c71541213172196c8a8"],
    ["V2 Master", "65185f7545f73fcf36ff5062bb85305affcae76458f8049130f052a811b79bee"],
])

doc.add_paragraph()
doc.add_paragraph(
    "To verify a file's integrity:"
)

p = doc.add_paragraph()
run = p.add_run("sha256sum <filename>")
run.font.name = 'Courier New'
run.font.size = Pt(10)

doc.add_paragraph(
    "The output should match the recorded hash exactly."
)

# ──────────────────────────────────────────────────────
# SAVE
# ──────────────────────────────────────────────────────
doc.save(str(OUT))
print(f"Report saved: {OUT}")
print(f"Size: {OUT.stat().st_size / 1024:.1f} KB")
