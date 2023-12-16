# Use Sherlock WGD annotations to annotate a tree.
# Input aligned fasta, like the other annotation scripts.
# Searches for the query gene name in the first column of a WGD annotation file
# Performs a reciprocal search of ohnologs in the query file to make the annotation

import argparse
from fasta_align_tree import create_phylogenetic_tree

def parse_metadata(metadata_file):
    print("Reading Ohnolog annotation file")
    metadata_dict = {}
    with open(metadata_file, 'r') as f:
        for line in f:
            entries = line.strip().split()
            gene_name = entries[0][:-5] # remove -WGDX
            wgd_info = entries[1:]
            metadata_dict[gene_name] = wgd_info
    #print(metadata_dict)
    return metadata_dict

''' e.g.
'65130.PTRED_209_T71800001295740123': ['65130.PTRED_209_T71800001293890467-WGD1', '65130.PTRED_209_T71800001294460098-WGD2', '65130.PTRED_209_T71800001293760101-WGD2']
'''

def process_query(query_file, metadata_dict, output_annotation):
    print("Annotating tree using Ohnolog annotations")
    # first make a list of all the header names in the query fasta file
    all_query_names = []
    with open(query_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                query_name = line.strip()
                gene_name = query_name[1:]
                all_query_names.append(gene_name)
    #print(all_query_names)

    with open(output_annotation, "w") as output:
    # now search for the gene name in the WGD annotation file
        output.write("seq\tohnoA\twgdA\tohnowgdA\tohnoB\twgdB\tohnowgdB\tohnoC\twgdC\tohnowgdC\tohnoD\twgdD\tohnowgdD\tohnoE\twgdE\tohnowgdE\tohnoF\twgdF\tohnowgdF\n")
        for gene_name in all_query_names:
            output_line_tmp = []
            if gene_name in metadata_dict:
                ohnologs = metadata_dict[gene_name]
                #print(gene_name, ohnologs)
                if ohnologs:
                    for entry_wgd in ohnologs:
                        #print(entry_wgd)
                        entry = entry_wgd.split('-')[0]  # don't consider -WGDX in the search
                        # print(entry)
                        if entry in all_query_names:  # there is a match. Look whether there is more than one ohnolog
                            print("Within-tree Ohnologue found! Ohno!")
                            output_line_tmp.append(entry + '\t' + entry_wgd[-4:] + '\t' + entry_wgd)
                    if output_line_tmp:
                        output_line_tmp = '\t'.join(output_line_tmp)
                    else:
                        output_line_tmp = "No_ohnolog_in_tree"
                else: # gene name is in annotation and has no ohnologs
                    output_line_tmp = "No_WGD_annotated"
            else: # gene name is not in annotation file
                output_line_tmp = "No_gene_in_annotation"
            output_line = gene_name +'\t' + output_line_tmp
            #print(output_line)
            output.write(output_line + '\n')
        print(f"Ohnologue annotations written to file {output_annotation}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process aligned fasta sequences.")
    parser.add_argument("query_file", help="FASTA alignment.")
    parser.add_argument("metadata_file", help="Ohnolog annotation file.")
    args = parser.parse_args()

    metadata_dict = parse_metadata(args.metadata_file)
    process_query(args.query_file, metadata_dict, args.query_file + ".inTreeOhnologs.annotation")
