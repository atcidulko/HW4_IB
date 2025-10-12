# BIO_UTILS: A lightweight educational bioinformatics toolkit written for Python course of **Institute of Bioinformatics**.  
It provides three main utilities:

1. **DNA/RNA sequence processing** (`run_dna_rna_tools`)  
2. **FASTQ record filtering** (`filter_fastq`)  
3. **File parsing utilities** (`bio_files_processor`)

This project was created exclusively for educational purposes and may not be copied or used for any commercial purpose.  
Its primary purpose is to illustrate to esteemed teachers and their assistants the concept of  
"one closed ear" — information provided to me does not leave my mind because it has no way out.

---
## Author
Aglaia Tcidulko  
Tg: @vesfir  
GitHub: https://github.com/atcidulko

--- 
## Overview

This project provides simple tools for working with DNA and RNA sequences, filtering FASTQ reads, and parsing biological files (FASTA, FASTQ, CSV).  
It includes functions for sequence manipulation (reverse, complement, transcription), GC content calculation, quality-based read filtering, and lightweight file handling.

All code has been improved based on feedback from the previous homework:  
	- clearer variable names  
	- proper use of docstrings instead of inline comments  
	- simplified and more readable logic  
	- adherence to PEP8 style guidelines  

### Features
- Reverse, complement, and transcription for DNA & RNA  
- GC content calculation (%)  
- Phred+33 quality decoding and mean score computation  
- Read filtering by GC%, length, and mean quality thresholds  
- FASTA and FASTQ parsing helpers  
- CSV writing and summarization tools  

---
## Project Structure: A Treasure Map

bio_utils/
│
├── utils/
│   ├── dna_rna_tools.py      # Essential DNA/RNA utilities
│   └── fastq_tools.py        # GC%, quality, and filtering helpers
│
├── bio_files_processor.py    # File parsers and format converters
├── main.py                   # Entry point script
├── README.md                 # This very manifesto
└── .gitignore                # Undeletable MacOS system file

-----
## Installation
git clone https://github.com/atcidulko/bio_utils.git  
cd bio_utils  
python3 main.py  

To fork the repository, replace <atcidulko> with your GitHub username.

---
## Usage examples

### DNA/RNA tools
from utils.dna_rna_tools import run_dna_rna_tools

print(run_dna_rna_tools("ATGC", "reverse"))
print(run_dna_rna_tools("ATGC", proc="transcribe"))
print(run_dna_rna_tools("AUGC", proc="complement"))

### FASTQ filtering
from utils.fastq_tools import filter_fastq

seqs_example = {
    "seq1": ("ATGC", "IIII"),
    "seq2": ("GGCC", "####"),
    "seq3": ("ATGCGT", "IIIIII"),
    "seq4": ("ATATAT", "IIIIII")
}

filtered = filter_fastq(
    seqs_example,
    gc_bounds=(40, 60),
    length_bounds=(5, 10),
    quality_threshold=30
)

print("Kept reads:", list(filtered.keys()))

### File processing
from bio_files_processor import parse_fasta, parse_fastq

fasta_records = parse_fasta("example_fasta.fasta")
fastq_records = parse_fastq("example_fastq.fastq")

print("FASTA records:", len(fasta_records))
print("FASTQ records:", len(fastq_records))

---
## Example output

Filtered sequences:  
seq3: seq=ATGCGT, qual=IIIIII  

---
## **Thank you for reading**

