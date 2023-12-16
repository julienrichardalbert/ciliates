import sys
import argparse
from Bio import SeqIO
import numpy as np

def extract_sequences(blast_results_file, db, output_file):
    print('Extracting sequences')
    # Step 1: Read BLAST results and extract sequence IDs
    with open(blast_results_file, "r") as blast_results:
        # Create a set to store the sequence IDs from BLAST results
        # BLAST results were in format: "-outfmt 6 qseqid sseqid ..."
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


def filter_long_sequences(input_file, output_file, input_multiplier):
    # Read sequences from the input FASTA file
    sequences = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

    # Calculate sequence lengths
    sequence_lengths = [len(seq) for seq in sequences.values()]

    # Calculate median, third quartile, and interquartile range
    median_length = np.median(sequence_lengths)
    third_quartile = np.percentile(sequence_lengths, 75)
    first_quartile = np.percentile(sequence_lengths, 25)
    iqr = third_quartile - first_quartile
    print(f'Blastp hits to filter: {len(sequence_lengths)}') # get the length of the lengths, genius
    print(f'median length: {median_length}')
    print(f'iqr: {iqr}')

    # Define length threshold (third quartile + 2 * IQR)
    length_threshold_minimum = 50
    length_threshold_iqr = third_quartile + input_multiplier * iqr
    length_threshold = length_threshold_iqr + length_threshold_minimum if (length_threshold_iqr - median_length < 50) else length_threshold_iqr
    print(f'Setting minimum length difference: {length_threshold_minimum}')
    print(f'Length threshold as a function of iqr: {length_threshold_iqr}')
    print(f'Length threshold used: {length_threshold}')

    output_log = input_file + '.size.log'
    # Print the number of input and output entries
    with open(output_log, 'w') as ofile:
        ofile.write(f'median length: {median_length}\n')
        ofile.write(f'iqr: {iqr}\n')
        ofile.write(f'Setting minimum length difference: {length_threshold_minimum}\n')
        ofile.write(f'Length threshold as a function of iqr: {length_threshold_iqr}\n')
        ofile.write(f'Length threshold used: {length_threshold}\n')

    # Write filtered sequences to the output FASTA file
    with open(output_file, "w") as output_fasta:
        for sequence_id, sequence in sequences.items():
            # Filter out sequences longer than the threshold
            if len(sequence) < length_threshold:
                # Write sequence ID and sequence data to the output FASTA file
                header_line = f">{sequence_id}\n"
                sequence_data = f"{sequence.seq}\n"
                output_fasta.write(header_line)
                output_fasta.write(sequence_data)
            else:
                # Log filtered sequences to stdout
                diff = int(len(sequence) - length_threshold)
                print(f"Sequence {sequence_id} is {diff}bp longer than the threshold and has been filtered out.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract amino acid sequence from a blastp output")
    parser.add_argument("blast_results_file", help="proteinA.fa.blastp")
    parser.add_argument("db", help="database.fa")
    parser.add_argument("-m", "--multiplier", help="Multiplication factor used for calculating sequence length thresholds", default='3')
    args = parser.parse_args()


    output_file = args.blast_results_file + ".fa"
    extract_sequences(args.blast_results_file, args.db, output_file)

    # write over file
    filter_long_sequences(output_file, output_file + ".filtLen", float(args.multiplier))
