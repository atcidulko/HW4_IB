from utils.fastq_tools import filter_fastq

if __name__ == "__main__":
    seqs_example = {
        "seq1": ("ATGC", "IIII"),       # len 4, GC=50%, avg_q=40
        "seq2": ("GGCC", "####"),       # len 4, GC=100%, avg_q=2
        "seq3": ("ATGCGT", "IIIIII"),   # len 6, GC=50%, avg_q=40
        "seq4": ("ATATAT", "IIIIII")    # len 6, GC=0%, avg_q=40
    }

    filtered = filter_fastq(
        seqs_example,
        length_bounds=(5, 10),
        gc_bounds=(40, 60),
        quality_threshold=30
    )

    print("Filtered sequences:")
    for name, (seq, qual) in filtered.items():
        print(f"{name}: seq={seq}, qual={qual}")
