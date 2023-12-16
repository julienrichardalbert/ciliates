# "trim" multiple sequence alignment with TrimAl
# output triplets amenable to ete3 evolTree
#
import argparse
import os
import subprocess
from Bio import SeqIO
from fasta_align_tree import create_phylogenetic_tree
from fasta_align_tree import save_phylogenetic_tree

def run_trimal(input_file, nogaps, automated1, noallgaps):
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
#    print(option_str)
    trimal_output = f"{input_file}.trimal.{option_str}"
#    print(trimal_output)
    html_trimal_output = f"{input_file}.trimal.{option_str}.html"

    command = ["trimal", "-in", input_file, "-out", trimal_output, "-htmlout", html_trimal_output] + options
    subprocess.run(command)
    return trimal_output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper script for trimal command.")
    parser.add_argument("input_file", help="Input file")
    parser.add_argument("-nogaps", action="store_true", help="Enable -nogaps option")
    parser.add_argument("-automated1", action="store_true", help="Enable -automated1 option")
    parser.add_argument("-noallgaps", action="store_true", help="Enable -noallgaps option")

    args = parser.parse_args()
    trimal_output = run_trimal(args.input_file, args.nogaps, args.automated1, args.noallgaps)
    # make a new tree using fixed alignment, no?
    phylogenetic_tree = create_phylogenetic_tree(alignment_file=trimal_output)
    tree_file = trimal_output + '.nwk'
    save_phylogenetic_tree(phylogenetic_tree, tree_file)
