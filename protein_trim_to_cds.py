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
from back_align import faSomeRecordPy
from log_progress import log


def run_trimal(input_file, db, newdb, args):
    log(f"Making smaller cds db for speed")
    faSomeRecordPy(
        query_file=input_file,
        target_file=db,
        output_file=newdb
    )

    option_str = str(args.trimal_strat)
#   trimal_output = f"{input_file}.trimal.{option_str}" # DONT DO THIS EVER
    trimal_output_cds = f"{input_file}.trimal.{option_str}.cds"
    html_trimal_output = f"{input_file}.trimal.{option_str}.html"
    options = str('-' + option_str)

    log(f'Trimming protein sequence and converting to cds using option: {args.trimal_strat}')
    command = ["trimal",
        "-in", input_file,
        "-out", trimal_output_cds,
        "-backtrans", input_file + '.cds.db', "-ignorestopcodon",
        "-htmlout", html_trimal_output,
        options]
    log(' '.join(map(str, command)))
    subprocess.run(command)
    return trimal_output_cds


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Wrapper script for trimal command.")
    parser.add_argument("input_file", help="Input file")
    parser.add_argument("db", help="CDS database")
    parser.add_argument('-trimal_strat', '--trimal_strat', action='store_true', default='automated1', help='TrimAl strategy. Please choose from: automated1, nogaps, noallgaps')

    args = parser.parse_args()
    trimal_output_cds = run_trimal(args.input_file, args.db, args.input_file + '.cds.db', args.nogaps, args.automated1, args.noallgaps)
    # make a new tree using fixed alignment, no?
    phylogenetic_tree = create_phylogenetic_tree(alignment_file=trimal_output_cds)
    tree_file = trimal_output_cds + '.nwk'
    save_phylogenetic_tree(phylogenetic_tree, tree_file)
