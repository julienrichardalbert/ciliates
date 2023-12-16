# Outputs ohnologs of single input protein/transcript/gene(s)

import argparse
import re
import subprocess
from extract_sequence_blast_hits import extract_sequences
from fasta_align_tree import perform_alignment
from fasta_align_tree import create_phylogenetic_tree
from fasta_align_tree import save_phylogenetic_tree

def search_and_output(input_file, ohno_srch, cds_srch, output_file):
    # Read the query sequence
    with open(input_file, 'r') as qfile:
        query_header = qfile.readline().strip()

    # Remove '>' from the start of the header
    query_header = query_header.lstrip('>')

    # trying to make it ignore T/P/G is impossible
    query_pattern = query_header

    # Read the search file
    match = []
    with open(ohno_srch, 'r') as sfile:
        for line in sfile:
            columns = line.strip().split('\t')
            # Check if the query pattern matches the first column
            if re.match(query_pattern, columns[0]):
                # Output the query and the rest of the columns as a list
                match = columns
                #print(len(match))

    # Print whether matches were found or not
    if len(match) > 1:
        print(f"{len(match)-1} ohnologues found!")

        # In order to work with FigTree, the first column must match the cds_srch nomenclature
        # In order to extract fasta seqs, the second column must match teh cds_srch nomenclature
        with open(output_file, "w") as output:
            output.write('seq\tseq2\tWGD\tfullseq\n')
            for hit in match:
                name, wgd = hit.split('-')
                output.write(name + '\t' + name + '\t' + wgd + '\t' + hit + '\n')

    elif len(match) == 1:
        print(f"An entry was found for query: {query_header} but it has no annotated Ohnolog. It's a Nonolog!")
    else:
        print(f"No matches found for query: {query_header}. Try replacing P/G in the header name with T!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="From a single-entry FASTA, use Sherlock annotations to define Ohnologs.")
    parser.add_argument("input_file", help="Fasta file with the query sequence in format TaxID.Species_strain_T00001X.")
    parser.add_argument("ohno_srch", help="WGD annotations from Sherlock.")
    parser.add_argument("cds_srch", help="Fasta containing all coding transcripts.")
    args = parser.parse_args()

    search_and_output(args.input_file, args.ohno_srch, args.cds_srch, args.input_file + '.ohnologs.txt')

    extract_sequences(blast_results_file=args.input_file + '.ohnologs.txt', db=args.cds_srch, output_file=args.input_file + '.ohnologs.fa')

    alignment_file = perform_alignment(args.input_file + '.ohnologs.fa', aligner="muscle")
    phylogenetic_tree = create_phylogenetic_tree(alignment_file)
    tree_file = args.input_file + '.ohnologs.fa' + '.aligned.nwk'
    save_phylogenetic_tree(phylogenetic_tree, tree_file)
