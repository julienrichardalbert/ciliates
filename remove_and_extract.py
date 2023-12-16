# REMOVE AND EXTRACT COLOURED TAXA FROM TREE
#
# Given a tree, with undesired taxa coloured, this script will remove
# those taxa from the original fasta file (outputs *.fasta.cleaned).
# This script is compatible with phylogenies made using RAxML, FastTree,
# and IQ-Tree.
#
# You can also extract sequences with a certain colour using the fourth argument.
#
# Usage: python remove_and_extract.py [fasta_file] [coloured_tree] [colour_code_to_remove] [colour_code_to_extract]
#
# Example 1: Remove red taxa and extract blue taxa
# python remove_and_extract.py example.fasta coloured.tre ff0000 0000ff
#
# Example 2: Just remove red taxa
# python remove_and_extract.py example.fasta coloured.tre ff0000
#
# Example 3: Just extract blue taxa
# python remove_and_extract.py example.fasta coloured.tre NA ff0000

import sys
import argparse
from fasta_align_tree import perform_alignment
from fasta_align_tree import create_phylogenetic_tree
from fasta_align_tree import save_phylogenetic_tree

colour = sys.argv[3]

# do you want to remove sequences?
if colour != 'NA':
    output = open(sys.argv[1] + '.cleaned', 'w')

# do you want to extract sequences?
extract = 'no'
if len(sys.argv) > 4:
    extract = 'yes'
    extract_colour_input = sys.argv[4]

    if extract_colour_input == 'ff0000' or extract_colour_input == 'red':
        extract_colour_name = 'red'
        extract_colour = 'ff0000'

    elif extract_colour_input == '00ff00' or extract_colour_input == 'green':
        extract_colour_name = 'green'
        extract_colour = '00ff00'

    elif extract_colour_input == '0000ff' or extract_colour_input == 'blue':
        extract_colour_name = 'blue'
        extract_colour = '0000ff'

    ext_out_name = sys.argv[1] + '.' + extract_colour_name + '.extract'
    ext_out = open(ext_out_name,'w')

# load in the fasta file
fasta_file = open(sys.argv[1], 'r').read().split('>')
fasta = []
for seq in fasta_file[1:]:
    fasta.append('>'+seq.split('\n')[0])
    sequence = ''
    for part in seq.split('\n')[1:]:
        sequence = sequence + part.strip()
    fasta.append(sequence)

# load in coloured tree file
tree = open(sys.argv[2], 'r').readlines()

# get lines with the coloured taxa info
coloured_taxa = []
extract_taxa = []

for lines in tree:
    if colour != 'NA':
        if colour in lines:
            coloured_taxa.append(lines.split("[")[0].strip("\t").strip("'").replace('@','_'))
    if extract == 'yes':
        if extract_colour in lines:
            extract_taxa.append(lines.split("[")[0].strip("\t").strip("'").replace('@','_'))

# remove and extract coloured taxa from the fasta file:
if colour != 'NA':
    n = 0
    while n < len(fasta):
        if fasta[n].startswith('>') and fasta[n].split(' ')[0].strip('>').strip('\n').replace('@','_') not in coloured_taxa:
            output.write(fasta[n].strip('\n') + '\n' + fasta[n+1].strip('\n') + '\n')
            n += 2
        elif fasta[n].startswith('>') and fasta[n].split(' ')[0].strip('>').strip('\n').replace('@','_') in extract_taxa:
            ext_out.write(fasta[n].strip('\n') + '\n' + fasta[n+1].strip('\n') + '\n')
            n += 2
        else:
            n += 2

if len(sys.argv) > 4:
    n = 0
    while n < len(fasta):
        if fasta[n].startswith('>') and fasta[n].split(' ')[0].strip('>').strip('\n').replace('@','_') in extract_taxa:
            ext_out.write(fasta[n].strip('\n') + '\n' + fasta[n+1].strip('\n') + '\n')
            n += 2
        else:
            n += 2

if colour != 'NA':
    output.close()
if len(sys.argv) > 4:
    ext_out.close()


# automatically align output sequences using Muscle and create phylo tree
alignment_file = perform_alignment(ext_out_name, aligner="muscle")
phylogenetic_tree = create_phylogenetic_tree(alignment_file)

if "trimal" in ext_out_name.lower(): # This shoud never happen
    print("WHY DO YOUR SEQUENCES YOU WANT TO EXTRACT HAVE TRIMAL IN THE NAME??? HMM? THINK ABOUT IT!")
    tree_file = ext_out_name + ".nwk"
else:
    tree_file = ext_out_name + ".aligned.nwk"

save_phylogenetic_tree(phylogenetic_tree, tree_file)
