# run the whole thing
# you must start with a list of single-copy genes
# how those are defined is probably really important, I'm just taking it easy and extracting all genes that don't have ohnologs.
# for now this works on a FASTA file... would be nice to simply input gene names
# TODO:  make high-level functions that can be called by yet another script. Rename this one something like "function amasser"
#
# REQUIREMENTS
# faSomeRecords
#

import os
import re
import argparse
import subprocess
import shutil

from Bio import SeqIO
from blastp import run_blastp, add_header_to_blastp
from filter_blast_hits import filter_best_hits
from extract_sequence_blast_hits import extract_sequences
from extract_sequence_blast_hits import filter_long_sequences

from fasta_align_tree import perform_alignment, create_phylogenetic_tree, save_phylogenetic_tree
from annotate_using_ohnolog import parse_metadata, process_query
from annotate_using_blast import run_blastp_uniprot, annotate_swiss
from annotate_using_taxid import annotate_taxonomy
# from annotate_using_taxid  # no functions but don't forget to run it
# A homologue identifying script here if I ever code it.
from protein_trim_to_cds import run_trimal
from evolution import file_exists, organize_your_life, get_species_name, calculate_dS, run_models, get_pvals, evol_graphs, modify_leaf_names_reroot, evol_graphs
from back_align  import faSomeRecordPy, back_align

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
from ete3 import NCBITaxa
import matplotlib.pyplot as plt
ncbi = NCBITaxa()



def parser_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="blastp | align | build tree | annotate | trim | cds | evolution")
    parser.add_argument("fasta_file", help="protein(s).fa")
    parser.add_argument("db_prot", help="protein database")
    parser.add_argument("db_cds", help="cds database")
    parser.add_argument("ohnologs", help="ohnolog annotation file.txt")
    parser.add_argument("-evol", "--evolution", action="store_true", default=False, help="Run the evolution analysis script.")
    parser.add_argument("-run1", "--keep_blastp", action="store_true", default=False, help="Keep all blastp results and stop.")
    parser.add_argument("-run2", "--runrest", action="store_true", default=False, help="From a manually curated tree, run the rest of the pipeline.")
    parser.add_argument("-runall", "--runall", action="store_true", default=False, help="Run all evolutionary models.")

    return parser.parse_args()

args = parser_arguments()
ori_dir = os.getcwd()

with open(args.fasta_file, 'r') as file:
    gene_records = SeqIO.parse(file, 'fasta')
    for gene_record in gene_records:
        gene_name = gene_record.id
        print(f"Processing {gene_name}...")

        os.makedirs(gene_name, exist_ok=True)
        with open(gene_name + '.txt', "w") as one_only_needs_a_name:
            one_only_needs_a_name.write('>'+gene_name)
        shutil.copy(gene_name + '.txt', gene_name)
        os.chdir(gene_name)

        faSomeRecordPy(query_file = gene_name + '.txt', target_file = args.db_prot, output_file = gene_name + '.fa')
        run_blastp(db=args.db_prot + '.db', query=gene_name + '.fa', evalue=str(0.1))
        add_header_to_blastp(gene_name + '.fa' + '.blastp')

        if args.keep_blastp:
            print("Keeping all blastp hits")
            extract_sequences (
                blast_results_file = gene_name + '.fa' + '.blastp',
                db = args.db_prot,
                output_file = gene_name + '.fa' + '.blastp' + '.fa'
            )
            filter_long_sequences (
                input_file  = gene_name + '.fa' + '.blastp' + '.fa',
                output_file = gene_name + '.fa' + '.blastp' + '.fa' + '.filtLen',
                input_multiplier = 3)
            hits = gene_name + '.fa' + '.blastp' + '.fa' + '.filtLen'
            alignment_file = perform_alignment(hits, aligner='muscle')
            phylogenetic_tree = create_phylogenetic_tree(alignment_file)
            save_phylogenetic_tree (phylogenetic_tree, output_file=hits + '.aligned' + '.nwk')
            # annotate the tree
            annotate_taxonomy(hits, hits + '.taxonomy.annotation')
            ohnolog_dict = parse_metadata(metadata_file=args.ohnologs)
            process_query(hits, ohnolog_dict, hits + '.inTreeOhnologs.annotation')
            run_blastp_uniprot("/Users/jra/Dropbox/ciliates/old/uniprot/uniprot_sprot.fasta.db", hits,  hits + ".uniprot.annotation")
            annotate_swiss(hits, hits + ".uniprot.annotation")
            continue # break out of the loop since we don't want to run the whole pipeline on ALL blast hits
                     # i.e. use this time to cull the tree and use that as input to a new pipeline run.
        # Proceed with filtering blastp results if -keep is not set

        if not filter_best_hits(
            input_file=gene_name + '.fa' + '.blastp',
            output_file=gene_name + '.fa' + '.blastp' + '.top'
        ):
            print(f"Skipping {gene_name} as it only returned 1 hit.")
            os.chdir(ori_dir)  # move back to the original directory
            continue

        extract_sequences (
            blast_results_file = gene_name + '.fa' + '.blastp' + '.top',
            db = args.db_prot,
            output_file = gene_name + '.fa' + '.blastp' + '.top' + '.fa'
        )

        filter_long_sequences (
            input_file = gene_name + '.fa' + '.blastp' + '.top' + '.fa',
            output_file = gene_name + '.fa' + '.blastp' + '.top' + '.fa' + '.filtLen',
            input_multiplier = 3)

        hits = gene_name + '.fa' + '.blastp' + '.top' + '.fa' + '.filtLen'

        alignment_file = perform_alignment(hits, aligner='muscle')
        phylogenetic_tree = create_phylogenetic_tree(alignment_file)

        save_phylogenetic_tree (
            phylogenetic_tree,
            output_file=hits + '.aligned' + '.nwk'
        )

        # annotate the tree
        annotate_taxonomy(hits, hits + '.taxonomy.annotation')
        ohnolog_dict = parse_metadata(metadata_file=args.ohnologs)
        process_query(hits, ohnolog_dict, hits + '.inTreeOhnologs.annotation')

        alignment = hits + '.aligned' # make it nicer in here

        trimal_output_cds = run_trimal (
            input_file=alignment,
            db=args.db_cds,
            newdb=alignment+'.cds.db',
            automated1=True,
            nogaps=False,
            noallgaps=False
        )

        phylogenetic_tree = create_phylogenetic_tree(alignment_file=trimal_output_cds)
        tree_file = trimal_output_cds + '.nwk'
        save_phylogenetic_tree(phylogenetic_tree, tree_file)

        # Measure positive selection
        if args.evolution:
            print("Measuring positive selection")
            trimal_output = alignment + '.trimal.automated1.cds' # this should be recoded in the trimal script
            tree = EvolTree(newick = tree_file, format=1)
            tree.link_to_alignment(trimal_output)
            current_directory = os.getcwd()
            tree.workdir = current_directory # set working directory where models will be saved

            print(tree) # because its fun to do
            calculate_dS(trimal_output, tree, '5888.PTET_51') # partial name of the reference leaf

            preference_dictionary = {'tet':1,
                'bia':2, 'dec':2, 'dod':2, 'jen':2, 'nov':2, 'oct':2,'pen':2, 'pri':2, 'qua':2, 'sex':2, 'son':2, 'tre': 2,
                                    'cau':3}

            if args.runall:
                models = ['M1', 'M2', 'M7', 'M8']
                run_models(trimal_output, tree, models)
                get_pvals(trimal_output, tree, 'M2', 'M1', 'p2') # alt, neg
                get_pvals(trimal_output, tree, 'M8', 'M7', 'p10')
                evol_graphs(trimal_output, tree, 'M2', 'M1', '') # alt, neg, suffix
                evol_graphs(trimal_output, tree, 'M8', 'M7', '')
                modify_leaf_names_reroot(tree, preference_dictionary) # change leaf name nomenclature from TaxID.gene to spe and reroot
                evol_graphs(trimal_output, tree, 'M2', 'M1', '_sp')
                evol_graphs(trimal_output, tree, 'M8', 'M7', '_sp')
                os.chdir(ori_dir)
                tree = ''
            else:
                models = ['M1', 'M2']
                run_models(trimal_output, tree, models)
                get_pvals(trimal_output, tree, 'M2', 'M1', 'p2') # alt, neg
                evol_graphs(trimal_output, tree, 'M2', 'M1', '') # alt, neg, suffix
                modify_leaf_names_reroot(tree, preference_dictionary) # change leaf name nomenclature from TaxID.gene to spe and reroot
                evol_graphs(trimal_output, tree, 'M2', 'M1', '_sp')
                os.chdir(ori_dir)
                tree = ''
        else:
            print("Skipping evolutionary analysis.")
            os.chdir(ori_dir)
