import argparse
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

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

    # Step 5: Save the tree to a file
    tree_file = fasta_file + ".aligned.nwk"
    Phylo.write(tree, tree_file, "newick")

    print("Phylogenetic tree created and saved to", tree_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a phylogenetic tree from a FASTA file.")
    parser.add_argument("fasta_file", help="proteinA.fa")
    args = parser.parse_args()

    create_phylogenetic_tree(args.fasta_file)
