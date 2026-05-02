"""
BIO_UTILS: A lightweight educational bioinformatics toolkit.

Provides biological sequence classes, FASTQ filtering, and file utilities.
Can be run from the command line to filter FASTQ files.

Usage:
    python main.py --input example.fastq --output filtered.fastq \\
        --gc-min 40 --gc-max 60 --len-min 50 --len-max 300 --quality 20
"""

import argparse
import logging
from abc import ABC, abstractmethod
from typing import Union, Tuple

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# ============================================================
# Logging setup
# ============================================================

def setup_logger(log_file: str = "bio_utils.log") -> logging.Logger:
    """Configure and return a logger that writes to a file and stdout.

    Args:
        log_file: path to the log file.

    Returns:
        Configured Logger instance.
    """
    logger = logging.getLogger("bio_utils")
    logger.setLevel(logging.DEBUG)

    # Remove existing handlers so we can switch to a new log file
    logger.handlers.clear()

    if True:
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        stream_handler.setFormatter(formatter)

        logger.addHandler(file_handler)
        logger.addHandler(stream_handler)

    return logger


logger = setup_logger()


# ============================================================
# Base class
# ============================================================

class BiologicalSequence(ABC):

    def __init__(self, sequence: str) -> None:
        self._seq = sequence
        if not self.check_alphabet():
            raise ValueError(
                f"Invalid alphabet for {type(self).__name__}: {sequence}"
            )

    def __len__(self) -> int:
        return len(self._seq)

    def __getitem__(self, index: Union[int, slice]) -> str:
        return self._seq[index]

    def __str__(self) -> str:
        return f"{type(self).__name__}({self._seq})"

    def __repr__(self) -> str:
        return self.__str__()

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


# ============================================================
# Nucleic acids
# ============================================================

class NucleicAcidSequence(BiologicalSequence):
    COMPLEMENT_MAP: str = None
    VALID_ALPHABET: set = None

    def check_alphabet(self) -> bool:
        if self.VALID_ALPHABET is None:
            raise NotImplementedError("VALID_ALPHABET must be defined in subclass")
        return set(self._seq.upper()).issubset(self.VALID_ALPHABET)

    def complement(self) -> "NucleicAcidSequence":
        if self.COMPLEMENT_MAP is None:
            raise NotImplementedError("COMPLEMENT_MAP must be defined in subclass")
        return type(self)(self._seq.translate(self.COMPLEMENT_MAP))

    def reverse(self) -> "NucleicAcidSequence":
        return type(self)(self._seq[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    COMPLEMENT_MAP = str.maketrans("ATGCatgc", "TACGtacg")
    VALID_ALPHABET = {"A", "T", "G", "C"}

    def transcribe(self) -> "RNASequence":
        rna_seq = self._seq.replace("T", "U").replace("t", "u")
        return RNASequence(rna_seq)


class RNASequence(NucleicAcidSequence):
    COMPLEMENT_MAP = str.maketrans("AUGCaugc", "UACGuacg")
    VALID_ALPHABET = {"A", "U", "G", "C"}


# ============================================================
# Amino acid sequence
# ============================================================

class AminoAcidSequence(BiologicalSequence):
    VALID_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy*-")

    def check_alphabet(self) -> bool:
        return set(self._seq).issubset(self.VALID_ALPHABET)

    def molecular_weight(self) -> float:
        """Approximate molecular weight in Da."""
        weights = {
            "A": 89, "R": 174, "N": 132, "D": 133, "C": 121,
            "E": 147, "Q": 146, "G": 75, "H": 155, "I": 131,
            "L": 131, "K": 146, "M": 149, "F": 165, "P": 115,
            "S": 105, "T": 119, "W": 204, "Y": 181, "V": 117,
        }
        return sum(weights.get(aa.upper(), 0) for aa in self._seq)


# ============================================================
# FastQ filtration via Biopython
# ============================================================

def filter_fastq(
    input_path: str,
    output_path: str,
    gc_bounds: Union[float, Tuple[float, float]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: float = 0,
    log_file: str = "bio_utils.log",
) -> None:
    """Filter a FastQ file by GC-content, length, and average Phred quality.

    Args:
        input_path: path to input .fastq file.
        output_path: path to output .fastq file with passing reads.
        gc_bounds: GC-content range in %. A single number means upper bound only.
        length_bounds: read length range. A single number means upper bound only.
        quality_threshold: minimum average Phred quality score.
        log_file: path to the log file.
    """
    global logger
    logger = setup_logger(log_file)
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, float(gc_bounds))
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, int(length_bounds))

    logger.info(
        "Starting FASTQ filtering: input=%s, output=%s, "
        "gc_bounds=%s, length_bounds=%s, quality_threshold=%s",
        input_path, output_path, gc_bounds, length_bounds, quality_threshold,
    )

    passed = []
    total = 0

    try:
        for record in SeqIO.parse(input_path, "fastq"):
            total += 1
            gc = gc_fraction(record.seq) * 100
            length = len(record)
            avg_quality = sum(record.letter_annotations["phred_quality"]) / length

            if (
                gc_bounds[0] <= gc <= gc_bounds[1]
                and length_bounds[0] <= length <= length_bounds[1]
                and avg_quality >= quality_threshold
            ):
                passed.append(record)

        SeqIO.write(passed, output_path, "fastq")

    except FileNotFoundError:
        logger.error("Input file not found: %s", input_path)
        raise

    except Exception as e:
        logger.error("Unexpected error during filtering: %s", e)
        raise

    logger.info(
        "Filtering complete: %d / %d reads passed. Written to %s",
        len(passed), total, output_path,
    )


# ============================================================
# Functions from bio_files_processor
# ============================================================

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    if output_fasta is None:
        output_fasta = "oneline_" + input_fasta

    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        seq_id = None
        seq_lines = []
        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id is not None:
                    outfile.write(f"{seq_id}\n{''.join(seq_lines)}\n")
                seq_id = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id is not None:
            outfile.write(f"{seq_id}\n{''.join(seq_lines)}")


def parse_blast_output(input_file: str, output_file: str) -> None:
    results = set()

    with open(input_file) as f:
        in_block = False
        for line in f:
            if "Sequences producing significant alignments:" in line:
                in_block = True
                continue
            if in_block:
                stripped = line.strip()
                if not stripped:
                    in_block = False
                    continue
                lower = stripped.lower()
                if lower.startswith("description") or (
                    "score" in lower and "e value" in lower
                ):
                    continue
                results.add(stripped.split()[0])
                in_block = False

    with open(output_file, "w") as out:
        for r in sorted(results):
            out.write(f"{r}\n")


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for FASTQ filtering.

    Returns:
        Namespace with parsed arguments.
    """
    parser = argparse.ArgumentParser(
        prog="bio_utils",
        description="Filter a FASTQ file by GC content, read length, and quality.",
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to the input FASTQ file.",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path to the output FASTQ file with passing reads.",
    )
    parser.add_argument(
        "--gc-min",
        type=float,
        default=0.0,
        help="Minimum GC content in %% (default: 0).",
    )
    parser.add_argument(
        "--gc-max",
        type=float,
        default=100.0,
        help="Maximum GC content in %% (default: 100).",
    )
    parser.add_argument(
        "--len-min",
        type=int,
        default=0,
        help="Minimum read length (default: 0).",
    )
    parser.add_argument(
        "--len-max",
        type=int,
        default=2**32,
        help="Maximum read length (default: no limit).",
    )
    parser.add_argument(
        "--quality",
        type=float,
        default=0.0,
        help="Minimum average Phred quality score (default: 0).",
    )
    parser.add_argument(
        "--log-file",
        default="bio_utils.log",
        help="Path to the log file (default: bio_utils.log).",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Re-setup logger with the user-specified log file
    global logger
    logger = setup_logger(args.log_file)

    filter_fastq(
        input_path=args.input,
        output_path=args.output,
        gc_bounds=(args.gc_min, args.gc_max),
        length_bounds=(args.len_min, args.len_max),
        quality_threshold=args.quality,
    )


if __name__ == "__main__":
    main()
