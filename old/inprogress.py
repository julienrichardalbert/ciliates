import subprocess
from Bio import SeqIO
import re
import argparse

def run_tblastn(input_protein, reference_genome, output_coding_sequence, output_table):
    # Step 1: Run tblastn to search for the input protein in the reference genome
    blast_command = f"tblastn -query {input_protein} -db {reference_genome} -out {output_table} -db_gencode 6 -outfmt '6 qseqid sseqid evalue bitscore sstart send sseq' -max_target_seqs 1"
    subprocess.run(blast_command, shell=True, check=True)

    # Step 2: Parse the BLAST results to get the subject sequence (best match)
    with open(output_table, "r") as blast_results:
        line = blast_results.readline()
        _, subject_sequence, _, _, start, end, subject_alignment = line.strip().split("\t")

    # Step 3: Extract the coding sequence from the reference genome
    with open(reference_genome, "r") as genome_file:
        coding_sequence = ""
        in_coding_region = False

        for record in SeqIO.parse(genome_file, "fasta"):
            if record.id == subject_sequence:
                in_coding_region = True
                coding_sequence += str(record.seq)

        coding_sequence = coding_sequence[int(start) - 1 : int(end)]

    # Step 4: Write the coding sequence to an output file
    with open(output_coding_sequence, "w") as output_file:
        output_file.write(f">{subject_sequence}\n{coding_sequence}")

    print("Coding sequence extracted and saved in", output_coding_sequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find coding sequence of a given protein in a reference genome using tblastn.")
    parser.add_argument("input_protein", help="Input protein sequence in FASTA format")
    parser.add_argument("reference_genome", help="Reference genome assembly in FASTA format")
    parser.add_argument("output_coding_sequence", help="Output file for coding sequence in FASTA format")
    parser.add_argument("output_table", help="Output file for BLAST results table")

    args = parser.parse_args()

    run_tblastn(args.input_protein, args.reference_genome, args.output_coding_sequence, args.output_table)
