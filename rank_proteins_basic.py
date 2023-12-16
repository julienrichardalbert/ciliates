import argparse
import pandas as pd
from Bio import SeqIO

def calculate_basic_score(sequence):
    basic_score = 0
    counts = {'R': 0, 'K': 0, 'H': 0, 'D': 0, 'E': 0}

    for aa in sequence:
        if aa in ['R', 'K', 'H']:
            basic_score += 1
            counts[aa] += 1
        elif aa in ['D', 'E']:
            basic_score -= 1
            counts[aa] += 1

    return basic_score, counts

def process_fasta(input_file):
    result_table = []

    with open(input_file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_name = record.id
            sequence = str(record.seq)

            basic_score, aa_counts = calculate_basic_score(sequence)

            # Append results to the table
            result_table.append({
                'protein_name': protein_name,
                'basic_score': basic_score,
                'R_count': aa_counts['R'],
                'K_count': aa_counts['K'],
                'H_count': aa_counts['H'],
                'D_count': aa_counts['D'],
                'E_count': aa_counts['E'],
                'protein_len': len(sequence),
                'basic_percent': round(100*basic_score/len(sequence), 2)
            })

    return result_table

def write_results_to_csv(result_table, output_file):
    df = pd.DataFrame(result_table)
    df.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rank proteins on their basic amino acid content.")
    parser.add_argument('-s', '--start_file', required=True, help='Protein FASTA file.')
    args = parser.parse_args()
    result_table = process_fasta(args.start_file)
    write_results_to_csv(result_table, args.start_file + '.basic.txt')
