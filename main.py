
import os
from collections import Counter

def analyze_statuses(fastq_file: str) -> Counter:
    """
    Analyze the BH status of reads in a FASTQ file.

    Each read header is expected to contain a status at the end,
    e.g., BH:ok, BH:failed, BH:changed.

    Args:
        fastq_file (str): Path to the FASTQ file.

    Returns:
        Counter: A dictionary-like object with counts of each status.
    """
    counter = Counter()
    with open(fastq_file, "r") as f:
        for i, line in enumerate(f):
            if i % 4 == 0:  # header line in FASTQ
                parts = line.strip().split()
                if len(parts) > 1:
                    status = parts[-1]
                    counter[status] += 1
    return counter

def main():
    """
    Main function for running the FASTQ BH status analysis.

    - Checks that the FASTQ file exists.
    - Counts the BH statuses.
    - Prints the results to the console.
    """
    fastq_path = "data/sample.fastq"  # <- change this to your file

    if not os.path.exists(fastq_path):
        print(f"File {fastq_path} not found!")
        return

    status_counts = analyze_statuses(fastq_path)

    print("BH statuses in the FASTQ file:")
    for status, count in status_counts.items():
        print(f"{status}: {count}")

    # Add other main.py functionality here if needed
    # For example, sequence processing, filtering, etc.

if __name__ == "__main__":
    main()

