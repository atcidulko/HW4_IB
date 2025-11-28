"""
Main module for running DNA/RNA tools and FASTQ filtering.

According to HW4/HW3 requirements, this file contains ONLY:
- imports
- run_dna_rna_tools
- filter_fastq (imported from utils)
No example code, no debugging code.
"""

from utils.dna_rna_tools import (
    transcribe_dna_to_rna,
    reverse,
    complement,
    reverse_complement,
    run_dna_rna_tools,
)

from utils.fastq_tools import filter_fastq


def main():
    """
    Entry point for the toolkit.
    At this stage, the main function only explains how to use the tools.
    It does NOT contain example run code (as required in HW4).
    """
    print("DNA/RNA tools and FASTQ filter imported successfully.")
    print("Use run_dna_rna_tools() or filter_fastq() in your workflow.")
    print("This main file intentionally contains no demo code.")


if __name__ == "__main__":
    main()
