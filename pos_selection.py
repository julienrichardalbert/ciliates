import argparse
import os
from config_loader import load_config
from log_progress import setup_logger, log

from Bio import SeqIO
from evol_chunks import genes_to_trees, trees_to_graphs
from ete3 import NCBITaxa
ncbi = NCBITaxa()

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
from ete3 import NCBITaxa
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a single instance of a combination of steps, which perform the steps: blastp | align | build tree | annotate | trim | cds | evolution")
    defaults = load_config()
    parser.add_argument('-s', '--start_file', required=True, help='Either a protein/gene name (e.g. 5888.PTET_51_T0310084) when paired with -genestart or a tree (culled tree ending in .nwk) when paired with -treestart.')
    parser.add_argument('-db_prot', '--db_prot',  default=defaults['db_prot'], help='full path to protein database.')
    parser.add_argument('-db_cds', '--db_cds', default=defaults['db_cds'], help='full path to cds database.')
    parser.add_argument('-db_ohno', '--db_ohno', default=defaults['db_ohno'], help='full path to ohnolog annotation file.')
    parser.add_argument('-db_uniprot', '--db_uniprot', default=defaults['db_uniprot'], help='full path to ohnolog annotation file.')

    parser.add_argument('-genestart', '--genestart', action='store_true', default=False, help='Keep all blastp results, make a tree, curate it and stop.')
    parser.add_argument('-treestart', '--treestart', action='store_true', default=False, help='Start from a manually curated tree, run the rest of the pipeline.')
    parser.add_argument('-fullpipe', '--fullpipe', action='store_true', default=False, help='Only keep top blastp hit per species and continue to the end.')

    parser.add_argument('-trimal_strat', '--trimal_strat', action='store_true', default=defaults['trimal_strat'], help='TrimAl strategy. Please choose from: automated1, nogaps, noallgaps')
    parser.add_argument('-lenmultiplier', '--lenmultiplier', type=float, default=defaults['lenmultiplier'], help='Base multiplier on which sequences are filtered. It is itself multiplied by the iqr of sizes.')
    parser.add_argument('-blaste', '--blaste', type=float, default=defaults['blaste'], help='BLASTp evalue to threshold results.')
    parser.add_argument('-evol', '--evolution', action='store_true', default=False, help='Run the evolution analysis script.')
    parser.add_argument('-refname', '--refname', action='store_true', default=defaults['refname'], help='String to match reference species name e.g. PTET_51.')
    parser.add_argument('-allmodels', '--allmodels', action='store_true', default=False, help='Run all evolutionary models.')

    args = parser.parse_args()
    ori_dir = os.getcwd()
    setup_logger(log_file=args.start_file)

    log('''
    ╭───────────────────────────────────╮
    │       Initiating a new run        │
    ╰───────────────────────────────────╯
    ''')
    log('Setting defaults')
    for key, value in defaults.items():
        log(f'{key}: {value}')
    log('')
    log(f'Input file: {args.start_file}')

    # start at the beginning
    if args.genestart:
        with open(args.start_file, 'r') as file:
            gene_record = SeqIO.read(file, 'fasta')
            gene_name = gene_record.id
            genes_to_trees(gene_name, ori_dir, args)

    # start in the middle and input a curated tree
    elif args.treestart:
        trees_to_graphs(alignment=args.start_file.replace('.nwk',  ''), args=args)

    # start with a gene, take the top hits for each species and continue to the end
    elif args.fullpipe:
        with open(args.start_file, 'r') as file:
            gene_record = SeqIO.read(file, 'fasta')
            gene_name = gene_record.id
            alignment_file = genes_to_trees(gene_name, ori_dir, args)
            trees_to_graphs(alignment=alignment_file, args=args)
    else:
        print('Please select either -genestart or -treestart\nType -h for more help.')  # DON'T TRY TO MAKE THE PIPELINE DO A FULL RUN. ALWAYS CURATE THE INITIAL ALIGNMENT!
    os.chdir(ori_dir)
    log('''>> Done! <<''')
