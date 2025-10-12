
def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """
    Convert a multi-line FASTA file so that each seq appears on a single line.

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
        # write the last sequence
        if seq_id is not None:
            outfile.write(f"{seq_id}\n{''.join(seq_lines)}\n")


def parse_blast_output(input_file, output_file):
    """
    Read a BLAST file and find the best matches

    Args:
        input_file (str): Path to the BLAST output text file.
        output_file (str): Path to save the parsed results.

    This function reads a file with BLAST results. It finds the first line under "Sequences that match well" for each query. 
    It saves all unique matches in alphabetical order in a new file.

    """
    results = set()
    with open(input_file) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        if "Sequences producing significant alignments:" in line:
            i += 1
            while i < len(lines) and lines[i].strip():
                # grab first word/column as description
                desc = lines[i].split()[0]
                results.add(desc)
                i += 1
        else:
            i += 1

    with open(output_file, "w") as out:
        for r in sorted(results):
            out.write(f"{r}\n")

