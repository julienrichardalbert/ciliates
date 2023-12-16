# translate CDS into protein
# kinda from https://biopython.org/docs/1.75/api/Bio.Seq.html

import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from log_progress import log


'''
codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'Q', 'TAG':'Q',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
'''


def translate_cds(input_file, output_file):
    log(f'Translating {input_file} using translate table #6')
    with open(input_file, 'r') as file:
        gene_records = SeqIO.parse(file, 'fasta')

        with open(output_file, "w") as ofile:
            for gene_record in gene_records:
                gene_name = gene_record.id
                cds_sequence = gene_record.seq
                prot_sequence = cds_sequence.translate(table='6') # table 6 for ciliates
                #print(gene_name)
                #print(prot_sequence)
                ofile.write('>' + gene_name + '\n')
                ofile.write(str(prot_sequence) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate DNA into amino acids like in the good ol' days!")
    parser.add_argument('-s', '--start_file', required=True, help='Input CDS FASTA.')
    parser.add_argument('-table', '--table', type=float, default=6, help='NCBI codon table.')

    args = parser.parse_args()
    translate_cds(args.start_file, args.start_file + '.aa')
