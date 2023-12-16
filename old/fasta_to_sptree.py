# create a species relatedness tree based on the header sequences of a fasta file
# NOTE: fasta file headers must be in the form "NCBITaxaID.ProteinID" (e.g. 9606.Q53XC5)

from ete3 import NCBITaxa, Tree
from Bio import SeqIO
import subprocess


# Initialize the NCBI Taxonomy database
ncbi = NCBITaxa()

# Function to extract taxIDs from a FASTA file and create an unrooted taxonomic tree
def create_unrooted_taxonomic_tree_from_fasta(fasta_file, output_tree_file):
    taxid_counts = {}
    taxIDs = []

    with open(fasta_file, 'r') as file:
        current_taxid = None
        for line in file:
            if line.startswith('>'):
                parts = line.strip('>').split('.')
                if len(parts) > 1:
                    current_taxid = int(parts[0])
                    taxIDs.append(current_taxid)
                    if current_taxid in taxid_counts:
                        taxid_counts[current_taxid] += 1
                    else:
                        taxid_counts[current_taxid] = 1

    if not taxIDs:
        print("No taxIDs found in the input file.")
        return None


    tree = ncbi.get_topology(taxIDs)
#    print(tree.get_ascii(attributes=["sci_name", "rank"]))

    # Save the unrooted tree in Newick format
    tree.write(outfile=output_tree_file)

#    sed_command = ["sed", "-I", ".ori", "-e", "s/:1//g", "-e", "s/)1/)/g", output_tree_file] # for DGINN spTree
#    subprocess.call(sed_command) # for DGINN spTree
    print(taxid_counts)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create an unrooted taxonomic tree from taxIDs in a FASTA file and save in Newick format")
    parser.add_argument("fasta_file", help="Input FASTA file with taxIDs in headers")
    args = parser.parse_args()
    output_tree_file = args.fasta_file + ".spTree"
    # Create the unrooted taxonomic tree and get taxID counts
    result = create_unrooted_taxonomic_tree_from_fasta(args.fasta_file, output_tree_file)
