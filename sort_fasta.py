from Bio import SeqIO
import argparse

def sort_fasta_by_length(input_file, output_file):
    # Read the sequences from the input FASTA file
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Sort the sequences based on length
    sorted_sequences = sorted(sequences, key=lambda x: len(x.seq))

    # Write the sorted sequences to the output FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(sorted_sequences, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sort a FASTA file based on sequence length")
    parser.add_argument("input_fasta", help="proteins.fa or dna.fa")
    args = parser.parse_args()

    sort_fasta_by_length(args.input_fasta, args.input_fasta + ".sorted")
