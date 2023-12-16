import argparse
from Bio import SeqIO

# Function to extract and write the longest isoform
def extract_longest_isoform(input_file, output_file):
    records = list(SeqIO.parse(input_file, "fasta"))

    # Dictionary to store the longest isoform for each unique prefix
    longest_isoforms = {}

    for record in records:
        sequence = str(record.seq)
        seq_id = record.id

        # Extract the unique prefix
        # >9986.NC_067389_1_cds_XP_051677352_1_36497
        # >9986.NC_067389
        prefix = seq_id.split('_')[0] + seq_id.split('_')[1]

        if prefix in longest_isoforms:
            # Check if the current sequence is longer
            if len(sequence) > len(longest_isoforms[prefix].seq):
                longest_isoforms[prefix] = record
        else:
            longest_isoforms[prefix] = record

    # Write the longest isoforms to the output file
    with open(output_file, 'w') as out_file:
        SeqIO.write(list(longest_isoforms.values()), out_file, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract the longest isoform from a FASTA of CDSs")
    parser.add_argument("fasta_file", help="CDS_geneA.fa")
    args = parser.parse_args()
    output_file = args.fasta_file + ".longestIso"
    extract_longest_isoform(args.fasta_file, output_file)
