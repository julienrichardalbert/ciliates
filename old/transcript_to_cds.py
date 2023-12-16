from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import six_frame_translations
from Bio.Data import CodonTable

# Define your custom codon table
ciliate_codon_table = {
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
        'TGC':'C', 'TGT':'C', 'TGA':'.', 'TGG':'W',
        '---':'-', 'nnn':'x', 'NNN':'X'
}



# Create a Biopython CodonTable object
custom_table_obj = CodonTable.CodonTable(forward_table=ciliate_codon_table)


def find_longest_orf(transcript_sequence, codon_table):
    seq = Seq(transcript_sequence)
    translations = []

    for frame in range(3):
        # Translate the sequence using the custom codon table
        translated_seq = seq[frame:].translate(codon_table)

        # Ensure the translated sequence is a multiple of three
        translated_seq = translated_seq[:len(translated_seq)//3 * 3]

        translations.append(str(translated_seq))

    # Find the longest translated sequence
    longest_orf = max(translations, key=len)
    return longest_orf

def process_fasta(input_file, output_file, codon_table):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            header = record.description
            sequence = str(record.seq)
            longest_orf = find_longest_orf(sequence, codon_table)

            # Write to the output file
            outfile.write(f">{header}\n{longest_orf}\n")

# Example usage:
input_fasta = "/Users/jra/Dropbox/ciliates/new_pipeline/input.fasta"
output_fasta = "/Users/jra/Dropbox/ciliates/new_pipeline/output.fasta"

# Use the custom codon table directly
process_fasta(input_fasta, output_fasta, custom_table_obj)
