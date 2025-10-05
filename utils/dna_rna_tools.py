def is_nucleic_acid(seq: str) -> bool:
    """
    Check if the sequence is a valid DNA or RNA sequence.
    Only A, T, G, C, U are allowed, and T and U cannot appear together.
    """
    if not seq:
        return False

    seq_upper = seq.upper()
    valid_nucleotides = "ATGCU"

    for nucleotide in seq_upper:
        if nucleotide not in valid_nucleotides:
            return False

    if "T" in seq_upper and "U" in seq_upper:
        return False

    return True


def transcribe(seq: str) -> str | None:
    """
    Transcribe DNA to RNA (A->U, T->A, G->C, C->G).
    Returns None for invalid sequences or if input is already RNA.
    """
    if not is_nucleic_acid(seq) or "U" in seq.upper():
        return None

    table = {"A": "U", "T": "A", "G": "C", "C": "G"}
    result = ""

    for nucleotide in seq:
        new_char = table.get(nucleotide.upper(), nucleotide)
        if nucleotide.islower():
            new_char = new_char.lower()
        result += new_char

    return result


def reverse(seq: str) -> str | None:
    """Return the reversed sequence. Returns None for invalid input."""
    if not is_nucleic_acid(seq):
        return None
    return seq[::-1]


def complement(seq: str) -> str | None:
    """
    Return the complementary sequence.
    Uses DNA table if sequence contains T, RNA table if sequence contains U.
    Returns None for invalid input.
    """
    if not is_nucleic_acid(seq):
        return None

    seq_upper = seq.upper()
    dna_table = {"A": "T", "T": "A", "G": "C", "C": "G"}
    rna_table = {"A": "U", "U": "A", "G": "C", "C": "G"}

    if "T" in seq_upper:
        table = dna_table
    else:
        table = rna_table

    result = ""
    for nucleotide in seq:
        new_char = table[nucleotide.upper()]
        if nucleotide.islower():
            new_char = new_char.lower()
        result += new_char

    return result


def reverse_complement(seq: str) -> str | None:
    """Return the reverse complement of a sequence."""
    comp_result = complement(seq)
    if comp_result is None:
        return None
    return reverse(comp_result)


def run_dna_rna_tools(*args) -> str | list | None:
    """
    Apply a procedure to one or more DNA/RNA sequences.
    
    The last argument is the procedure name.
    Returns a single result for one sequence or a list for multiple sequences.
    """
    if not args:
        return None

    *sequences, procedure = args
    results = []

    func_map = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    if procedure not in func_map:
        return None

    for seq in sequences:
        results.append(func_map[procedure](seq))

    return results[0] if len(results) == 1 else results

