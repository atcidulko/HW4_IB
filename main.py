from abc import ABC, abstractmethod
from typing import Union, Tuple

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


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
# Aminoacid sequence
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
) -> None:
    """
    Filter FastQ file by GC-content, length, and average Phred quality.

    Args:
        input_path: path to input .fastq file
        output_path: path to output .fastq file
        gc_bounds: GC-content range in % 
        length_bounds: read length range
        quality_threshold: minimum average Phred quality score
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, float(gc_bounds))
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, int(length_bounds))

    passed = []

    for record in SeqIO.parse(input_path, "fastq"):
        gc = gc_fraction(record.seq) * 100
        length = len(record)
        avg_quality = sum(record.letter_annotations["phred_quality"]) / length

        if (gc_bounds[0] <= gc <= gc_bounds[1]
                and length_bounds[0] <= length <= length_bounds[1]
                and avg_quality >= quality_threshold):
            passed.append(record)

    SeqIO.write(passed, output_path, "fastq")


# ============================================================
# Functions from bio_files_processor (they can stay)
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
