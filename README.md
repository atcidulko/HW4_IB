# BIO_UTILS: A lightweight educational bioinformatics toolkit written for Python course of **Institute of Bioinformatics**.

It provides three main utilities:

1. **Biological sequence processing** (classes `DNASequence`, `RNASequence`, `AminoAcidSequence`)
2. **FASTQ record filtering** (`filter_fastq`)
3. **File parsing utilities** (`convert_multiline_fasta_to_oneline`, `parse_blast_output`)

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

This project has been refactored from a set of standalone functions into a proper object-oriented toolkit.
It now uses classes to represent biological sequences, with shared logic defined in abstract base classes.
FASTQ filtering has been reimplemented using **Biopython**.

All code has been improved based on feedback from previous homeworks:
- clearer variable names
- proper use of docstrings instead of inline comments
- simplified and more readable logic
- adherence to PEP8 style guidelines
- full OOP refactoring: abstraction, inheritance, polymorphism

### Features
- Abstract base class `BiologicalSequence` defining a shared interface
- `DNASequence`: reverse, complement, reverse_complement, transcribe → RNASequence
- `RNASequence`: reverse, complement, reverse_complement
- `AminoAcidSequence`: alphabet validation, molecular weight estimation
- FASTQ filtering by GC%, length, and mean Phred quality — powered by Biopython
- FASTA format converter and BLAST output parser

---
## Project Structure: A Treasure Map
```
hw4_ib/
│
├── main.py                # All classes, filter_fastq, and file utilities
├── requirements.txt       # biopython
├── README.md              # This very manifesto
├── example.fastq          # Test data
└── .gitignore             # Undeletable MacOS system file
```

---
## Installation
```bash
git clone https://github.com/atcidulko/HW4_IB.git
cd HW4_IB
pip install -r requirements.txt
```

---
## Usage examples

### Biological sequences
```python
from main import DNASequence, RNASequence, AminoAcidSequence

dna = DNASequence("ATGC")
print(dna.complement())        # DNASequence(TACG)
print(dna.reverse())           # DNASequence(CGTA)
print(dna.transcribe())        # RNASequence(AUGC)

rna = RNASequence("AUGC")
print(rna.reverse_complement()) # RNASequence(GCAU)

protein = AminoAcidSequence("ACDEF")
print(protein.molecular_weight())  # approximate weight in Da
```

### FASTQ filtering
```python
from main import filter_fastq

filter_fastq(
    input_path="example.fastq",
    output_path="filtered.fastq",
    gc_bounds=(40, 60),
    length_bounds=(50, 300),
    quality_threshold=20
)
```

### File utilities
```python
from main import convert_multiline_fasta_to_oneline, parse_blast_output

convert_multiline_fasta_to_oneline("input.fasta", "output.fasta")
parse_blast_output("blast_results.txt", "best_hits.txt")
```

---
## **Thank you for reading**
