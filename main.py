# main.py
from utils.fastq_utils import read_fastq, write_fastq, filter_fastq

# Parameters
input_fastq = "example_data/test.fastq"
output_fastq = "filtered_test.fastq"

# Filtering criteria
length_bounds = (5, 1000)       # min and max length
gc_bounds = (40, 60)            # min and max GC content percentage
quality_threshold = 30           # minimum average Phred33 quality

# 1. Read FASTQ
seqs = read_fastq(input_fastq)

# 2. Filter sequences
filtered = filter_fastq(
    seqs,
    gc_bounds=gc_bounds,
    length_bounds=length_bounds,
    quality_threshold=quality_threshold
)

# 3. Write filtered FASTQ
write_fastq(filtered, output_fastq)

