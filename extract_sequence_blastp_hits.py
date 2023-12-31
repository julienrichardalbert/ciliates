import sys
import argparse
from Bio import SeqIO

def extract_sequences(blast_results_file, db, output_file):
    # Step 1: Read BLASTP results and extract sequence IDs
    with open(blast_results_file, "r") as blast_results:
        # Create a set to store the sequence IDs from BLAST results
        sequence_ids = {line.split()[1] for line in blast_results if not line.startswith("#")}

    # Step 2: Read protein database sequences and extract hits
    # Create an index of sequences in the protein database
    sequence_index = SeqIO.index(db, "fasta")
    # Create a dictionary to store the sequences of hits from the database
    sequences = {}

    # Loop through the sequence IDs and extract the corresponding sequences
    for seq_id in sequence_ids:
        if seq_id in sequence_index:
            sequences[seq_id] = str(sequence_index[seq_id].seq)

    # Step 3: Write the extracted sequences to a new FASTA file
    with open(output_file, "w") as output_fasta:
        for sequence_id, sequence in sequences.items():
            # Write sequence ID and sequence data to the output FASTA file
            header_line = f">{sequence_id}\n"
            sequence_data = f"{sequence}\n"
            output_fasta.write(header_line)
            output_fasta.write(sequence_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract amino acid sequence from a blastp output")
    parser.add_argument("blast_results_file", help="proteinA.fa.blastp")
    parser.add_argument("db", help="database.fa")
    args = parser.parse_args()

    output_file = args.blast_results_file + ".fa"
    extract_sequences(args.blast_results_file, args.db, args.blast_results_file + ".fa")
