# annotate tree using sequence length

import argparse
from Bio import SeqIO


def annotate_seq_size(input_file, output_file):
    print('Annotating using sequence lengths')
    sequences = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    with open(output_file, "w") as output_table:
        output_table.write('#name\tsize\n') # header
        for seq_id, sequence in sequences.items(): # for each sequence name and sequence
                output_table.write(f'{seq_id}\t{len(sequence.seq)}\n') # write the name and the sequence length

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate input alignment on sequence length")
    parser.add_argument("input_seqs", help="proteinA.fa.blastp")
    args = parser.parse_args()

    output_file = args.input_seqs + ".size.annotation"
    annotate_seq_size(args.input_seqs, output_file)
