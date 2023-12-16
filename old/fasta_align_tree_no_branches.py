import argparse
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import subprocess

def create_phylogenetic_tree(fasta_file):
    # Step 1: Read the protein sequences from the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    # Step 2: Perform multiple sequence alignment using MAFFT
    print("Multiple sequence alignment using MAFFT on file", fasta_file)
    alignment_file = fasta_file + ".aligned"
    mafft_cline = MafftCommandline(input=fasta_file)
    stdout, stderr = mafft_cline()
    with open(alignment_file, "w") as handle:
        handle.write(stdout)

    # Step 3: Calculate the distance matrix
    print("Calculating the distance matrix using BioPython on file", alignment_file)
    alignment = AlignIO.read(alignment_file, "fasta")
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    # Step 4: Construct the phylogenetic tree
    print("Creating phylogenetic tree using BioPython")
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Simplify the tree by setting branch lengths to a constant value (e.g., 1.0)
    for clade in tree.find_clades(order="postorder"):
        clade.branch_length = 1.0

    # Save the simplified tree to a file
    simplified_tree_file = fasta_file + ".spTree.nwk"
    Phylo.write(tree, simplified_tree_file, "newick")

    # get Ugly

    # Define the sed command as a list of arguments
    sed_command = ["sed", "-I", ".ori", "-e", "s/:1\.00000//g", "-e", "s/Inner.//g", simplified_tree_file]
    # Run the sed command
    subprocess.call(sed_command)

    print("Simplified phylogenetic tree created and saved to", simplified_tree_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a phylogenetic tree from a FASTA file.")
    parser.add_argument("fasta_file", help="proteinA.fa")
    args = parser.parse_args()

    create_phylogenetic_tree(args.fasta_file)
