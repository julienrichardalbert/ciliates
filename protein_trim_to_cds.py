# "trim" multiple sequence alignment with TrimAl
# output triplets amenable to ete3 evolTree
#
import argparse
import os
import re
import subprocess
from Bio import SeqIO
from fasta_align_tree import create_phylogenetic_tree
from fasta_align_tree import save_phylogenetic_tree

def run_trimal(input_file, db, newdb, nogaps, automated1, noallgaps):

    print(f"Making smaller cds db for speed")
    with open(input_file, 'r') as names_in, open('temp.txt', 'w') as names_out:
        for line in names_in:
            if line.startswith('>'):
                gene_name = re.sub('>', '', line.strip())
                names_out.write(gene_name + '\n')
    command = ['faSomeRecords', db, 'temp.txt', newdb]
    subprocess.call(command)
    os.remove('temp.txt')

    options = []
    # If no optional arguments are provided, set -automated1 to be True
    if not any([nogaps, automated1, noallgaps]):
        automated1 = True

    if nogaps:
        options.append("-nogaps")
    if automated1:
        options.append("-automated1")
    if noallgaps:
        options.append("-noallgaps")
    print(f"Trimming alignment file with option: {options}")

    option_str = "_".join(options).replace("-", "") # there's probably a cleaner way
#   trimal_output = f"{input_file}.trimal.{option_str}" # DONT DO THIS
    trimal_output_cds = f"{input_file}.trimal.{option_str}.cds"
    html_trimal_output = f"{input_file}.trimal.{option_str}.html"

    print("Trimming protein sequence AND converting to cds")
    command = ["trimal",
        "-in", input_file,
        "-out", trimal_output_cds,
        "-backtrans", input_file + '.cds.db', "-ignorestopcodon",
        "-htmlout", html_trimal_output
    ] + options
    subprocess.run(command)
    return trimal_output_cds


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper script for trimal command.")
    parser.add_argument("input_file", help="Input file")
    parser.add_argument("db", help="CDS database")
    parser.add_argument("-nogaps", action="store_true", help="Enable -nogaps option")
    parser.add_argument("-automated1", action="store_true", help="Enable -automated1 option")
    parser.add_argument("-noallgaps", action="store_true", help="Enable -noallgaps option")

    args = parser.parse_args()
    trimal_output_cds = run_trimal(args.input_file, args.db, args.input_file + '.cds.db', args.nogaps, args.automated1, args.noallgaps)
    # make a new tree using fixed alignment, no?
    phylogenetic_tree = create_phylogenetic_tree(alignment_file=trimal_output_cds)
    tree_file = trimal_output_cds + '.nwk'
    save_phylogenetic_tree(phylogenetic_tree, tree_file)
