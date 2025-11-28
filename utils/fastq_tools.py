from typing import Dict, Tuple


SequenceName = str
Sequence = str
QualityString = str
SeqRecord = Tuple[Sequence, QualityString]
SeqDict = Dict[SequenceName, SeqRecord]


def gc_content(sequence: str) -> float:
    """
    Calculate GC percentage of a sequence.

    Args:
        sequence: Nucleotide sequence (DNA).

    Returns:
        GC-content in percent.
    """
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return 100.0 * gc_count / len(sequence)


def avg_quality(quality_string: str) -> float:
    """
    Calculate average Phred33 quality score of a sequence.

    Args:
        quality_string: ASCII-encoded quality string (Phred+33).

    Returns:
        Average quality score.
    """
    if not quality_string:
        return 0.0
    scores = [ord(char) - 33 for char in quality_string]
    return sum(scores) / len(scores)


def _normalize_bounds(
    bounds: int | float | Tuple[float, float],
    max_default: float,
) -> Tuple[float, float]:
    """
    Normalize bounds argument.

    If `bounds` is a single number → interpret as (0, bounds).
    If `bounds` is a tuple of two numbers → return as (min, max).

    Args:
        bounds: single number or (min, max) pair.
        max_default: default upper bound if input is incorrect.

    Returns:
        (min_bound, max_bound)
    """
    if isinstance(bounds, (int, float)):
        return 0.0, float(bounds)

    if isinstance(bounds, tuple) and len(bounds) == 2:
        lower, upper = bounds
        return float(lower), float(upper)

    return 0.0, float(max_default)


def filter_fastq(
    seqs: SeqDict,
    gc_bounds: float | Tuple[float, float] = (0, 100),
    length_bounds: int | Tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0,
) -> SeqDict:
    """
    Filter FASTQ records by length, GC-content and average quality.

    Args:
        seqs:
            Dictionary with read names as keys and (sequence, quality) as values.
        gc_bounds:
            GC-content interval in percent. Examples:
            - (20, 80) → keep reads with GC between 20 and 80%;
            - 44.4 → keep reads with GC <= 44.4%.
        length_bounds:
            Length interval for filtering. Analogous to gc_bounds:
            - (50, 150);
            - 100 → keep reads with len <= 100.
        quality_threshold:
            Minimum average Phred33 quality score.

    Returns:
        Dictionary with reads that pass all filters.
    """
    gc_min, gc_max = _normalize_bounds(gc_bounds, max_default=100.0)

    length_min_float, length_max_float = _normalize_bounds(
        length_bounds,
        max_default=float(2**32),
    )
    length_min = int(length_min_float)
    length_max = int(length_max_float)

    filtered_sequences: SeqDict = {}

    for read_name, (sequence, quality_string) in seqs.items():
        sequence_length = len(sequence)

        if not (length_min <= sequence_length <= length_max):
            continue

        gc_value = gc_content(sequence)
        if not (gc_min <= gc_value <= gc_max):
            continue

        average_quality = avg_quality(quality_string)
        if average_quality < quality_threshold:
            continue

        filtered_sequences[read_name] = (sequence, quality_string)

    return filtered_sequences
