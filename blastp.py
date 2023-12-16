# from a single protein FASTA, go fishing using blastp and return hits in .txt format
# UPDATE 30 NOV 2023: extract the sequences after performing blastp
# UPDATE 05 DEC 2023: added option to include the -evalue in blastp (default=10)

import argparse
import subprocess
from extract_sequence_blast_hits import extract_sequences
from extract_sequence_blast_hits import filter_long_sequences

def run_blastp(db, query, evalue):
    print(f'Running blastp with query {query}')
    output = query + ".blastp"

    cline_blastp = [
        "blastp",
        "-db", db,
        "-query", query,
        "-out", output,
        "-evalue", evalue,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    ]
    subprocess.run(cline_blastp)

def add_header_to_blastp(output):
    header = "#qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
    with open(output, "r+") as file:
        content = file.read()
        file.seek(0, 0)
        file.write(header + content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run blastp and add header to the output. Do not include .db in database name")
    parser.add_argument("query", help="proteinA.fa")
    parser.add_argument("db", help="database.fa")
    parser.add_argument("-e", "--evalue", help="BLASTp evalue. Hits < evalue results will be kept.", default='10')
    args = parser.parse_args()

    run_blastp(args.db + ".db", args.query, args.evalue)
    add_header_to_blastp(args.query + ".blastp")

    # run extract_sequence_blast_hits.py
    extract_output = args.query + ".blastp" + ".fa"
    extract_sequences(args.query + ".blastp", args.db, extract_output)
    filter_long_sequences(extract_output, extract_output + ".filtLen", input_multiplier=3)
