# ANNOTATE A PHYLOGENY WITH TAXONOMY INFORMATION
#
# Create an annotation file for FigTree from taxa ids in the headers of a fasta
#
# NOTE: fasta file headers must be in the form "NCBITaxaID.ProteinID" (e.g. 9606.Q53XC5)
#

import argparse
from ete3 import NCBITaxa

def annotate_taxonomy(input_file, output_file):
    print("Annotating using TaxID")
    ncbi = NCBITaxa()

    lines = open(input_file, 'r').readlines()
    otu_d = {}

    for i in lines:
        if i.startswith('>'):
            otu_d[i.strip('>').strip()] = [i.split('.')[0].strip('>').strip()]

    with open(output_file, 'w') as output:
        output.write('otu\ttaxaid\tspecies\tsupergroup\tdomain\n')

        for i in otu_d:
            taxa = i.split('.')[0]
            try:
                otu_d[i].append(ncbi.get_taxid_translator([int(taxa)])[int(taxa)])
            except:
                otu_d[i].append('NA')
            try:
                otu_d[i].append(str(ncbi.get_taxid_translator([ncbi.get_lineage(int(taxa))[3]])[ncbi.get_lineage(int(taxa))[3]]))
            except:
                otu_d[i].append('NA')
            try:
                domain = 'NA'
                if 2759 in ncbi.get_lineage(int(taxa)):
                    domain = 'Eukaryota'
                elif 10239 in ncbi.get_lineage(int(taxa)):
                    domain = 'Virus'
                elif 2 in ncbi.get_lineage(int(taxa)):
                    domain = 'Bacteria'
                elif 2157 in ncbi.get_lineage(int(taxa)):
                    domain = 'Archaea'
                else:
                    domain = 'Virus'
                otu_d[i].append(domain)
            except:
                otu_d[i].append('NA')

            output.write(i + '\t' + '\t'.join(otu_d[i]) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate a phylogeny with taxonomy information.")
    parser.add_argument("input_file", help="Path to the input file")
    args = parser.parse_args()
    output_file = args.input_file +  '.taxonomy.annotation'
    annotate_taxonomy(args.input_file, output_file)
