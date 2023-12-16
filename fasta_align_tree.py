# From a starting FASTA file, perform MAFFT alignment and create a phylogenetic tree in newick format
# Requires muscle and/or mafft to be installed and in the PATH
#
# EDIT 10 NOV 2023 The function skips MAFFT alignment if it detects the file is pre-aligned (contains "trimal" in filename, for now)
# EDIT 30 NOV 2023 Replace MAFFT alignment with Muscle alignment. Include MAFFT as optional parameter
# EDIT 30 NOV 2023 Break function into two: alignment, tree building

import argparse
import os
import subprocess
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo


def perform_alignment(fasta_file, aligner):
    # Check if the input file contains "trimal"
    skip_alignment = "trimal" in fasta_file.lower()

    # Read the FASTA sequences
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    # Perform multiple sequence alignment using the specified aligner (if not skipped)
    if not skip_alignment:
        print(f"Multiple sequence alignment using {aligner}")
        alignment_file = fasta_file + ".aligned"

        if aligner.lower() == "mafft":
            align_command = f"mafft --quiet {fasta_file} > {alignment_file}"
        elif aligner.lower() == "muscle":
            align_command = f"muscle -align {fasta_file} -output {alignment_file}"
        else:
            raise ValueError("Invalid aligner choice. Use 'mafft' or 'muscle'.")

        subprocess.run(align_command, shell=True)
    else:
        print(f"Skipping {aligner} alignment for pre-aligned file")


    return alignment_file if not skip_alignment else fasta_file

def create_phylogenetic_tree(alignment_file):
    # Calculate the distance matrix
    print("Calculating the distance matrix using BioPython")
    alignment = AlignIO.read(alignment_file, "fasta")
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    # Step 4: Construct the phylogenetic tree
    print("Creating phylogenetic tree using BioPython")
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    return tree

def save_phylogenetic_tree(tree, output_file):
    # Step 5: Save the tree to a file
    Phylo.write(tree, output_file, "newick")
    print("Phylogenetic tree created and saved to", output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Multiple sequence alignment followed by phylogenetic tree construction from a FASTA.")
    parser.add_argument("fasta_file", help="proteins.fa, CDS.fa, or aligned file (e.g., proteins.trimal.fa)")
    parser.add_argument("-a", "--aligner", choices=["mafft", "muscle"], default="muscle", help="Alignment method: 'mafft' or 'muscle'")
    args = parser.parse_args()

    alignment_file = perform_alignment(args.fasta_file, args.aligner)
    phylogenetic_tree = create_phylogenetic_tree(alignment_file)

    if "trimal" in args.fasta_file.lower():
        tree_file = args.fasta_file + ".nwk"
    else:
        tree_file = args.fasta_file + ".aligned.nwk"

    save_phylogenetic_tree(phylogenetic_tree, tree_file)
