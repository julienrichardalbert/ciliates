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
from translate import translate_cds
from evolution import calculate_dS, run_models, get_pvals_branch, get_pvals, evol_graphs, modify_leaf_names_reroot
from back_align import faSomeRecordPy, back_align

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
from ete3 import NCBITaxa
import matplotlib.pyplot as plt

ncbi = NCBITaxa()

def genes_to_trees(gene_name, ori_dir, args):
    log(f"Processing {gene_name}...")
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



    if args.fullpipe: # skip manual tree curation
        seqs_to_align = gene_name + '.fa' + '.blastp.e' + evalue + '.top'
        filter_best_hits(gene_name + '.fa' + '.blastp.e' + evalue, seqs_to_align, args)
    else:
        seqs_to_align = gene_name + '.fa' + '.blastp.e' + evalue

    extract_sequences(
        blast_results_file=seqs_to_align,
        db=args.db_prot,
        output_file=seqs_to_align + '.fa'
    )

    filter_long_sequences(
        input_file=seqs_to_align + '.fa',
        output_file=seqs_to_align + '.fa' + '.filtLen',
        input_multiplier=args.lenmultiplier)

    hits = seqs_to_align + '.fa' + '.filtLen'
    alignment_file = perform_alignment(hits, aligner='muscle')
    phylogenetic_tree = create_phylogenetic_tree(alignment_file)
    save_phylogenetic_tree(phylogenetic_tree, output_file=hits + '.aligned' + '.nwk')
    annotate_taxonomy(hits, hits + '.taxonomy.annotation')
    ohnolog_dict = parse_metadata(metadata_file=args.db_ohno)
    process_query(hits, ohnolog_dict, hits + '.inTreeOhnologs.annotation')
    run_blastp_uniprot(args.db_uniprot, hits, hits + '.uniprot.annotation')
    annotate_swiss(hits, hits + '.uniprot.annotation')
    annotate_seq_size(hits, hits + '.size.annotation')

    # combine annotation files into one
    annotation_files = [
        hits + '.taxonomy.annotation',
        hits + '.inTreeOhnologs.annotation',
        hits + '.blast.annotation',
        hits + '.size.annotation',
    ]

    # This doesn't actually work. Damn. Keep making it to shame yourself.
    annot_dfs = [pd.read_csv(annot, delimiter='\t') for annot in annotation_files]

    # Combine annotation files into one DataFrame
    result = pd.concat(annot_dfs, axis=1, join='inner')
    result = result.fillna("NA")
    output_annotations = hits + '.combined.annotation.txt'
    result.to_csv(output_annotations, sep='\t', index=False)

    return alignment_file # in order to input into function 2/2



def trees_to_graphs(alignment, args):

    # use TrimAl to trim the multifasta
    trimal_output_cds = run_trimal(
        input_file=alignment,
        db=args.db_cds,
        newdb=alignment + '.cds.db',
        args=args
    )

    translate_cds(
        input_file  = alignment + '.trimal.automated1.cds',
        output_file = alignment + '.trimal.automated1.aa'
    )
    #THIS ISNT DOING WHAT YOU THINK IT dOES. IT HSOULD RETURN JUST ONE AA SEQ!
    faSomeRecordPy(
        query_file  = args.start_file,
        target_file = alignment + '.trimal.automated1.aa',
        output_file = alignment + '.trimal.automated1.aa.starter'
    )
    back_align(
        original = alignment,
        trimmed  = alignment + '.trimal.automated1.aa.starter',
        output   = alignment + '.trimal.automated1.aa.starter.pairAln'
    )

    # from run_trimal function: trimal_output_cds = f"{input_file}.trimal.{option_str}.cds"

    phylogenetic_tree = create_phylogenetic_tree(alignment_file=alignment + '.trimal.automated1.cds')
    tree_file = alignment + '.trimal.automated1.cds.nwk'
    save_phylogenetic_tree(phylogenetic_tree, tree_file)

    # Measure positive selection
    if args.evolution or args.fullpipe:
        log("Measuring positive selection")
        tree = EvolTree(newick=tree_file, format=1)
        tree.link_to_alignment(alignment + '.trimal.automated1.cds')
        tree.workdir = os.getcwd()  # set working directory where models will be saved

        log(tree)  # because its fun to do
        calculate_dS(alignment + '.trimal.automated1.cds', tree, args.refname)  # partial name of the reference leaf

        # this has to be this way, unless I load the config.ini into the evolution script.
        try:
            preference_dictionary = load_config()['preference_dictionary']
        except KeyError:
            print("preference_dictionary not found in defaults")
            preference_dictionary = {'tet': 1,
                                    'bia': 2, 'dec': 2, 'dod': 2, 'jen': 3, 'nov': 2, 'oct': 2, 'pen': 2, 'pri': 2, 'qua': 2,
                                    'sex': 3, 'son': 3, 'tre': 2,
                                    'cau': 4}

        if args.allmodels:
            models = ['M1', 'M2', 'M7', 'M8', 'b_free', 'M0']
            run_models(alignment + '.trimal.automated1.cds', tree, models)
            get_pvals_branch(alignment + '.trimal.automated1.cds', tree, 'b_free', 'M0')
            get_pvals(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', 'p2')  # alt, neg
            get_pvals(alignment + '.trimal.automated1.cds', tree, 'M8', 'M7', 'p10')
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '')  # alt, neg, suffix
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M8', 'M7', '')

            log('Changing leaf names to species name')
            modify_leaf_names_reroot(tree, preference_dictionary)  # change leaf name nomenclature from TaxID.gene to spe and reroot
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '_sp')
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M8', 'M7', '_sp')
            tree = ''
        else:
            models = ['M1', 'M2', 'b_free', 'M0']
            run_models(alignment + '.trimal.automated1.cds', tree, models)
            get_pvals_branch(alignment + '.trimal.automated1.cds', tree, 'b_free', 'M0' )
            get_pvals(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', 'p2')  # alt, neg
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '')  # alt, neg, suffix
            modify_leaf_names_reroot(tree, preference_dictionary)  # change leaf name nomenclature from TaxID.gene to spe and reroot
            evol_graphs(alignment + '.trimal.automated1.cds', tree, 'M2', 'M1', '_sp')
            tree = ''
    else:
        log("Skipping evolutionary analysis.")
