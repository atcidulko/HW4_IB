def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """
    Convert a multi-line FASTA file so that each sequence appears on a single line.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str, optional): Path for the output FASTA file.
            If not provided, the output file will be named "oneline_" + input filename.

    This function reads a FASTA file where sequences can be split across multiple lines,
    and writes a new FASTA file where each sequence is joined into a single line.
    """
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


def parse_blast_output(input_file, output_file):
    """
    Parse a BLAST output file and extract best matches.

    Logic:
    - Find the line: 'Sequences producing significant alignments:'
    - Then find the first non-empty and non-technical line below it
      (not starting with 'Description' and not the 'Score    E value' header).
    - Take the first token in that line (sequence ID).
    - Collect all such IDs from all blocks.
    - Write unique IDs sorted alphabetically into output_file.
    """
    results = set()

    with open(input_file) as f:
        in_block = False  # are we inside a 'Sequences producing significant alignments:' section?

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

                hit_id = stripped.split()[0]
                results.add(hit_id)

                in_block = False

    with open(output_file, "w") as out:
        for r in sorted(results):
            out.write(f"{r}\n")

