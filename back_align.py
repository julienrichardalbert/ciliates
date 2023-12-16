# do pairwise alignment of the trimmed protein with the original

import argparse
import os
import subprocess
from Bio import SeqIO


def faSomeRecordPy(query_file, target_file, output_file):
    # Read the query sequence
    with open(query_file, 'r') as qfile:
        query_record = SeqIO.read(qfile, 'fasta')
        #print(query_record)

    # Search for the query sequence in the target file
    with open(target_file, 'r') as tfile, open(output_file, 'w') as ofile:
        target_sequences = SeqIO.to_dict(SeqIO.parse(tfile, 'fasta'))

        if query_record.id in target_sequences:
            # Write the matching sequence to the output file
            SeqIO.write(target_sequences[query_record.id], ofile, 'fasta')
            print(f"Sequence found and written to {output_file}")
        else:
            print("Sequence not found in the target file.")

def back_align(original, trimmed):
    output = trimmed + '.pairAln'
    command = ['needle', trimmed, original, output, '-gapopen', '0', '-gapextend', '0', '-aformat', 'pair']
    subprocess.run(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pairwise alignment of a trimmed protein sequence with its original.")
    parser.add_argument("input_trim", help="Input trimmed protein. Can be a multi FASTA.")
    parser.add_argument("input_ori", help="Input original protein. Must be a single-sequence FASTA.")
#   maybe add an optional argument to add a Database of known protein domains?
#   parser.add_argument("-db", "--database", help="Database of known protein domains", default='/path/to/db')
    args = parser.parse_args()


    faSomeRecordPy(args.input_ori, args.input_trim, args.input_trim + '.processed')
    back_align(args.input_ori, args.input_trim + '.processed')
