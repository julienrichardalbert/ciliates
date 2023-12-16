# do pairwise alignment of the trimmed protein with the original

import argparse
import os
import subprocess
import re
from Bio import SeqIO
from log_progress import log


def faSomeRecordPy(query_file, target_file, output_file):
    with open(query_file, 'r') as names_in, open('temp.txt', 'w') as ofile:
        for line in names_in:
            if line.startswith('>'):
                gene_name = re.sub('>', '', line.strip())
                ofile.write(gene_name + '\n')
    command = ['faSomeRecords', target_file, 'temp.txt', output_file]
    subprocess.call(command)
    os.remove('temp.txt')



def back_align(original, trimmed, output):
    log('Pairwise alignment of trimmed protein sequence to the original')
    command = ['needle', trimmed, original, output, '-gapopen', '0', '-gapextend', '0', '-aformat', 'pair']
    subprocess.run(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pairwise alignment of a trimmed protein sequence with its original.")
    parser.add_argument("input_trim", help="Input trimmed protein. Can be a multi FASTA.")
    parser.add_argument("input_ori", help="Input original protein. Must be a single-sequence FASTA.")
#   maybe add an optional argument to add a Database of known protein domains?
#   parser.add_argument("-db", "--database", help="Database of known protein domains", default='/path/to/db')
    args = parser.parse_args()


    faSomeRecordPy(args.input_ori, args.input_trim, args.input_trim + '.processed')
    back_align(args.input_ori, args.input_trim + '.processed', args.input_trim + '.processed' + '.pairAln')
