"""
Tests for bio_utils: filter_fastq, sequence classes, logging, and CLI.

Run with:
    pytest tests/test_bio_utils.py -v
"""

import logging
import os
import sys
import tempfile

import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Make sure main.py is importable from the repo root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from main import (
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    filter_fastq,
    setup_logger,
)

EXAMPLE_FASTQ = os.path.join(os.path.dirname(__file__), "..", "example.fastq")


# ============================================================
# Helpers
# ============================================================

def make_fastq(records: list[tuple], path: str) -> None:
    """Write a minimal FASTQ file from a list of (name, seq, quals) tuples."""
    bio_records = []
    for name, seq, quals in records:
        r = SeqRecord(Seq(seq), id=name, description="")
        r.letter_annotations["phred_quality"] = quals
        bio_records.append(r)
    SeqIO.write(bio_records, path, "fastq")


def count_reads(path: str) -> int:
    return sum(1 for _ in SeqIO.parse(path, "fastq"))


# ============================================================
# Tests: DNASequence
# ============================================================

class TestDNASequence:

    def test_complement(self):
        dna = DNASequence("ATGC")
        assert str(dna.complement()) == "DNASequence(TACG)"

    def test_reverse_complement(self):
        dna = DNASequence("ATGC")
        assert str(dna.reverse_complement()) == "DNASequence(GCAT)"

    def test_transcribe(self):
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        assert isinstance(rna, RNASequence)
        assert str(rna) == "RNASequence(AUGC)"

    def test_invalid_alphabet_raises(self):
        """Invalid characters must raise ValueError."""
        with pytest.raises(ValueError):
            DNASequence("ATGX")


# ============================================================
# Tests: AminoAcidSequence
# ============================================================

class TestAminoAcidSequence:

    def test_molecular_weight(self):
        aa = AminoAcidSequence("AG")
        # A=89, G=75
        assert aa.molecular_weight() == 164

    def test_invalid_alphabet_raises(self):
        """Invalid amino acid characters must raise ValueError."""
        with pytest.raises(ValueError):
            AminoAcidSequence("ACBZ")


# ============================================================
# Tests: filter_fastq — filtering logic
# ============================================================

class TestFilterFastq:

    def test_all_reads_pass_with_default_params(self, tmp_path):
        """With default (permissive) thresholds all reads should pass."""
        out = str(tmp_path / "out.fastq")
        filter_fastq(EXAMPLE_FASTQ, out)
        assert count_reads(out) == count_reads(EXAMPLE_FASTQ)

    def test_strict_quality_filters_all_out(self, tmp_path):
        """Quality threshold above maximum possible Phred should filter everything."""
        out = str(tmp_path / "out.fastq")
        filter_fastq(EXAMPLE_FASTQ, out, quality_threshold=100)
        assert count_reads(out) == 0

    def test_gc_filter(self, tmp_path):
        """Only reads within gc_bounds should pass."""
        inp = str(tmp_path / "in.fastq")
        out = str(tmp_path / "out.fastq")
        # seq GGGG → GC=100%, seq AAAA → GC=0%
        make_fastq(
            [
                ("high_gc", "GGGG", [40, 40, 40, 40]),
                ("low_gc",  "AAAA", [40, 40, 40, 40]),
            ],
            inp,
        )
        filter_fastq(inp, out, gc_bounds=(50, 100))
        ids = [r.id for r in SeqIO.parse(out, "fastq")]
        assert ids == ["high_gc"]

    def test_length_filter(self, tmp_path):
        """Only reads within length_bounds should pass."""
        inp = str(tmp_path / "in.fastq")
        out = str(tmp_path / "out.fastq")
        make_fastq(
            [
                ("short", "AT",     [40, 40]),
                ("long",  "ATGCAT", [40] * 6),
            ],
            inp,
        )
        filter_fastq(inp, out, length_bounds=(5, 100))
        ids = [r.id for r in SeqIO.parse(out, "fastq")]
        assert ids == ["long"]


# ============================================================
# Tests: filter_fastq — error handling
# ============================================================

class TestFilterFastqErrors:

    def test_missing_input_raises(self, tmp_path):
        """FileNotFoundError must be raised for a non-existent input file."""
        out = str(tmp_path / "out.fastq")
        with pytest.raises(FileNotFoundError):
            filter_fastq("nonexistent_file.fastq", out)


# ============================================================
# Tests: file I/O — output file is created and valid
# ============================================================

class TestFilterFastqFileIO:

    def test_output_file_is_created(self, tmp_path):
        """Output file must exist after filtering."""
        out = str(tmp_path / "result.fastq")
        filter_fastq(EXAMPLE_FASTQ, out)
        assert os.path.exists(out)

    def test_output_is_valid_fastq(self, tmp_path):
        """Output must be a parseable FASTQ file."""
        out = str(tmp_path / "result.fastq")
        filter_fastq(EXAMPLE_FASTQ, out)
        # SeqIO.parse raises if the file is malformed
        records = list(SeqIO.parse(out, "fastq"))
        assert isinstance(records, list)


# ============================================================
# Tests: logging
# ============================================================

class TestLogging:

    def test_log_file_is_created(self, tmp_path):
        """Logger must create the log file."""
        log_path = str(tmp_path / "test.log")
        setup_logger(log_path)
        out = str(tmp_path / "out.fastq")
        filter_fastq(EXAMPLE_FASTQ, out, log_file=log_path)
        assert os.path.exists(log_path)

    def test_log_file_contains_info(self, tmp_path):
        """Log file must contain at least one INFO message after filtering."""
        log_path = str(tmp_path / "test.log")
        out = str(tmp_path / "out.fastq")
        filter_fastq(EXAMPLE_FASTQ, out, log_file=log_path)
        content = open(log_path).read()
        assert "INFO" in content

    def test_log_file_contains_error_on_bad_input(self, tmp_path):
        """Log file must contain an ERROR message when input file is missing."""
        log_path = str(tmp_path / "test.log")
        out = str(tmp_path / "out.fastq")
        with pytest.raises(FileNotFoundError):
            filter_fastq("no_such_file.fastq", out, log_file=log_path)
        content = open(log_path).read()
        assert "ERROR" in content
