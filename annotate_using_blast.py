# ANNOTATE A TREE WITH BLAST HITS
# Adapted from https://github.com/nickatirwin/Figtree-Analysis/blob/master/blast_annotation.py
# Annotate sequences within a tree using protein annotations from SWISS-PROT

'''
set up database :
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz # download swiss-prot and unzip
gunzip uniprot_sprot.fasta.gz
sed -i 's/ /_/g' uniprot_sprot.fasta # replace spaces with underscore so that the full annotation is included in the blast output
makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot.fasta.db -dbtype prot # make a blast database
set location of blast database at the bottom of this script:
run_blastp_uniprot("/Users/jra/Dropbox/ciliates/uniprot/uniprot_sprot.fasta.db", args.fasta_file)
'''

import sys
import argparse
import subprocess

def run_blastp_uniprot(db, query, output_file):
    print("Annotating tree using BLAST")
    print("Running blastp...")
    cline_blastp = [
        "blastp",
        "-db", db,
        "-query", query,
        "-out", output_file,
        "-outfmt", "6",
        "-evalue",  "1e-5"
    ]
    subprocess.run(cline_blastp)

def annotate_swiss(fasta_file, blastp_results):
    print("Parsing results and annotating tree")
    # 1. load in fasta file
    fasta = open(fasta_file,'r').read()
    seq_d = {}
    for line in fasta.split('>')[1:]:
        seq_d[line.split('\n')[0].strip().strip('>')] = []

    # 2. load in blast file
    blast = open(blastp_results,'r').readlines()
    blast_d = {}
    evalue_d = {}
    for line in blast:
        if line.split('\t')[0].strip() in blast_d:
            pass
        else:
            blast_d[line.split('\t')[0].strip()] = line.split('\t')[1].split('_',1)[1].split('_',1)[1].split('OS=')[0].strip('_')
            evalue_d[line.split('\t')[0].strip()] = line.split('\t')[-2].strip()

    # 3. Output annotation
    out = open(fasta_file + '.blast.annotation', 'w')
    out.write('seq\tseq2\tprotein\tevalue\n')
    for s in list(seq_d.keys()):
        try:
            out.write(s + '\t' + s + '\t' + blast_d[s] + '\t' + evalue_d[s] + '\n')
        except:
            out.write(s + '\t' + s + '\tno_annotation\tNA\n')
    out.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate protein sequences using SWISS-PROT annotations")

    defaults = load_config()
    parser.add_argument("fasta_file", help="proteinA.fa")
    parser.add_argument('-db_uniprot', '--db_uniprot',  default=defaults['db_uniprot'], help='full path to uniprot database.')
    args = parser.parse_args()

    run_blastp_uniprot(args.db_uniprot, args.fasta_file,  args.fasta_file + ".uniprot.out")
    annotate_swiss(args.fasta_file, args.fasta_file + ".uniprot.out")
