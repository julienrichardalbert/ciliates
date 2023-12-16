import os
import re
import argparse
import subprocess
import shutil
import logging
import datetime
import pandas as pd
from Bio import SeqIO

from config_loader import load_config
from log_progress import setup_logger, log
from blastp import run_blastp, add_header_to_blastp
from filter_blast_hits import filter_best_hits
from extract_sequence_blast_hits import extract_sequences
from extract_sequence_blast_hits import filter_long_sequences
from fasta_align_tree import perform_alignment, create_phylogenetic_tree, save_phylogenetic_tree
from annotate_using_ohnolog import parse_metadata, process_query
from annotate_using_blast import run_blastp_uniprot, annotate_swiss
from annotate_using_taxid import annotate_taxonomy
from annotate_using_size import annotate_seq_size
from protein_trim_to_cds import run_trimal
from evolution import calculate_dS, run_models, get_pvals, evol_graphs, modify_leaf_names_reroot
from back_align import faSomeRecordPy

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
from ete3 import NCBITaxa
import matplotlib.pyplot as plt

ncbi = NCBITaxa()

def genes_to_trees(gene_name, ori_dir, args):
    logging.info(f"Processing {gene_name}...")
    os.makedirs(gene_name, exist_ok=True)
    with open(gene_name + '.txt', "w") as one_only_needs_a_name:
        one_only_needs_a_name.write('>'+gene_name)
    shutil.copy(gene_name + '.txt', gene_name)
    os.chdir(gene_name)

    evalue = str(args.blaste)

    faSomeRecordPy(query_file=gene_name + '.txt', target_file=args.db_prot, output_file=gene_name + '.fa')
    run_blastp(db=args.db_prot + '.db', query=gene_name + '.fa', evalue=evalue)
    os.rename(gene_name + '.fa.blastp', gene_name + '.fa.blastp.e' + evalue) # this should be recoded into blastp.py but my other pipeline would suffer.
    add_header_to_blastp(gene_name + '.fa' + '.blastp.e' + evalue)

    logging.info(f"Making a tree from {len(gene_name + '.fa' + '.blastp.e' + evalue)} blastp hits")

    extract_sequences(
        blast_results_file=gene_name + '.fa' + '.blastp.e' + evalue,
        db=args.db_prot,
        output_file=gene_name + '.fa' + '.blastp.e' + evalue + '.fa'
    )

    filter_long_sequences(
        input_file=gene_name + '.fa' + '.blastp.e' + evalue + '.fa',
        output_file=gene_name + '.fa' + '.blastp.e' + evalue + '.fa' + '.filtLen',
        input_multiplier=3)

    hits = gene_name + '.fa' + '.blastp.e' + evalue + '.fa' + '.filtLen'
    alignment_file = perform_alignment(hits, aligner='muscle')
    phylogenetic_tree = create_phylogenetic_tree(alignment_file)
    save_phylogenetic_tree(phylogenetic_tree, output_file=hits + '.aligned' + '.nwk')
    annotate_taxonomy(hits, hits + '.taxonomy.annotation')
    ohnolog_dict = parse_metadata(metadata_file=args.db_ohno)
    process_query(hits, ohnolog_dict, hits + '.inTreeOhnologs.annotation')
    run_blastp_uniprot('/Users/jra/Dropbox/ciliates/old/uniprot/uniprot_sprot.fasta.db', hits,
        hits + '.uniprot.annotation')
    annotate_swiss(hits, hits + '.uniprot.annotation')
    annotate_seq_size(hits, hits + '.size.annotation')

    # combine annotation files into one
    annotation_files = [
        hits + '.taxonomy.annotation',
        hits + '.inTreeOhnologs.annotation',
        hits + '.blast.annotation',
        hits + '.size.annotation',
    ]

    # Read each file into a DataFrame and merge them on the shared first column
    annot_dfs = [pd.read_csv(annot, delimiter='\t') for annot in annotation_files]

    # Combine annotation files into one DataFrame
    result = pd.concat(annot_dfs, axis=1, join='inner')
    result = result.fillna("NA")
    output_annotations = hits + '.combined.annotation.txt'
    result.to_csv(output_annotations, sep='\t', index=False)



def trees_to_graphs(alignment, input_arguments):

    # use TrimAl to trim the multifasta
    trimal_output_cds = run_trimal(
        input_file=alignment,
        db=args.db_cds,
        newdb=alignment + '.cds.db',
        automated1=True,
        nogaps=False,
        noallgaps=False
    )

    # from run_trimal function: trimal_output_cds = f"{input_file}.trimal.{option_str}.cds"

    phylogenetic_tree = create_phylogenetic_tree(alignment_file=alignment + '.trimal.automated1.cds')
    tree_file = alignment + '.trimal.automated1.cds.nwk'
    save_phylogenetic_tree(phylogenetic_tree, tree_file)

    # Measure positive selection
    if args.evolution:
        logging.info("Measuring positive selection")
        tree = EvolTree(newick=tree_file, format=1)
        tree.link_to_alignment(alignment + '.trimal.automated1.cds')
        tree.workdir = os.getcwd()  # set working directory where models will be saved

        logging.info(tree)  # because its fun to do
        calculate_dS(alignment + '.trimal.automated1.cds', tree, '5888.PTET_51')  # partial name of the reference leaf

        # this has to be this way, unless I load the config.ini into the evolution script.
        try:
            preference_dictionary = defaults['preference_dictionary']
        except KeyError:
            print("preference_dictionary not found in defaults")
            preference_dictionary = {'tet': 1,
                                    'bia': 2, 'dec': 2, 'dod': 2, 'jen': 3, 'nov': 2, 'oct': 2, 'pen': 2, 'pri': 2, 'qua': 2,
                                    'sex': 3, 'son': 3, 'tre': 2,
                                    'cau': 4}

        if args.allmodels:
            models = ['M1', 'M2', 'M7', 'M8', 'fb']
            run_models(alignment + '.trimal.automated1.cds', tree, models)
            get_pvals(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', 'p2')  # alt, neg
            get_pvals(alignment + '.trimal.automated1.cds', tree, 'M8', 'M7', 'p10')
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '')  # alt, neg, suffix
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M8', 'M7', '')
            modify_leaf_names_reroot(tree, preference_dictionary)  # change leaf name nomenclature from TaxID.gene to spe and reroot
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '_sp')
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M8', 'M7', '_sp')
            os.chdir(ori_dir)
            tree = ''
        else:
            models = ['M1', 'M2']
            run_models(alignment + '.trimal.automated1.cds', tree, models)
            get_pvals(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', 'p2')  # alt, neg
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '')  # alt, neg, suffix
            modify_leaf_names_reroot(tree, preference_dictionary)  # change leaf name nomenclature from TaxID.gene to spe and reroot
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '_sp')
            tree = ''
    else:
        logging.info("Skipping evolutionary analysis.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a single instance of a combination of steps, which perform the steps: blastp | align | build tree | annotate | trim | cds | evolution")
    defaults = load_config()
    parser.add_argument('-s', '--start_file', required=True, help='Either a protein/gene name (e.g. 5888.PTET_51_T0310084) when paired with -genestart or a tree (culled tree ending in .nwk) when paired with -treestart.')
    parser.add_argument('-db_prot', '--db_prot',  default=defaults['db_prot'], help='full path to protein database.')
    parser.add_argument('-db_cds', '--db_cds', default=defaults['db_cds'], help='full path to cds database.')
    parser.add_argument('-db_ohno', '--db_ohno', default=defaults['db_ohno'], help='full path to ohnolog annotation file.')

    parser.add_argument('-genestart', '--genestart', action='store_true', default=False, help='Keep all blastp results, make a tree, curate it and stop.')
    parser.add_argument('-treestart', '--treestart', action='store_true', default=False, help='Start from a manually curated tree, run the rest of the pipeline.')

    parser.add_argument('-blaste', '--blaste', type=float, default=defaults['blaste'], help='BLASTp evalue to threshold results.')
    parser.add_argument('-evol', '--evolution', action='store_true', default=False, help='Run the evolution analysis script.')
    parser.add_argument('-allmodels', '--allmodels', action='store_true', default=False, help='Run all evolutionary models.')

    args = parser.parse_args()
    ori_dir = os.getcwd()
    setup_logger(log_file=args.start_file)

    log('Setting defaults')
    for key, value in defaults.items():
        log(f'{key}: {value}')
    log('Starting pipeline')

    # start at the beginning
    if args.genestart:
        print(f'Input file: {args.start_file}')
        with open(args.start_file, 'r') as file:
            gene_record = SeqIO.read(file, 'fasta')
            gene_name = gene_record.id
            hits = genes_to_trees(gene_name, ori_dir, args)

    # start in the middle and input a curated tree
    elif args.treestart:
        print(f'Input file: {args.start_file}')
        trees_to_graphs(alignment=args.start_file.replace('.nwk',  ''), input_arguments=args)

    else:
        print('Please select either -genestart or -treestart\nType -h for more help.')  # DON'T TRY TO MAKE THE PIPELINE DO A FULL RUN. ALWAYS CURATE THE INITIAL ALIGNMENT!
    os.chdir(ori_dir)
    log('Script done')
