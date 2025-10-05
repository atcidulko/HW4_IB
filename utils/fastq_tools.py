def filter_fastq(
    seqs: dict,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0
) -> dict:
    """
    Filters FASTQ sequences based on GC content, length, and quality.

    Args:
        seqs (dict): dictionary with sequence names as keys and (sequence, quality) as values
        gc_bounds (tuple): (min_gc, max_gc) percentage bounds for GC content
        length_bounds (tuple): (min_len, max_len) bounds for seq length
        quality_threshold (int): min average Phred33 quality score

    Returns:
        dict: sequences that pass all filters, with names as keys and (sequence, quality) as values
    """

    def gc_content(seq: str) -> float:
        """Calculate the GC percentage of a sequence."""
        return 100 * (seq.count("G") + seq.count("C")) / len(seq)

    def avg_quality(qual: str) -> float:
        """Calculate average Phred33 quality score of a sequence."""
        return sum(ord(c) - 33 for c in qual) / len(qual)

    filtered_sequences = {}

    for name, (seq, qual) in seqs.items():
        if not seq:  # skip empty seq
            continue

        # Check seq length
        if not (length_bounds[0] <= len(seq) <= length_bounds[1]):
            continue

        # Check GC
        gc = gc_content(seq)
        if not (gc_bounds[0] <= gc <= gc_bounds[1]):
            continue

        # Check average quality
        avg_q = avg_quality(qual)
        if avg_q < quality_threshold:
            continue

        # Passed all filters
        filtered_sequences[name] = (seq, qual)

    return filtered_sequences



