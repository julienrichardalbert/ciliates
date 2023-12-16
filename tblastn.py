# Search a translated nucleotide database using a protein query
import argparse
import subprocess

def run_tblastn(db, query):
    output = query + ".tblastn"

    cline_tblastn = [
        "tblastn",
        "-db", db,
        "-query", query,
        "-out", output,
        "-outfmt", "6 qseqid sseqid pident sstart send length mismatch gapopen qstart qend sstart send evalue bitscore"
    ]
    subprocess.run(cline_tblastn)

def add_header_to_tblastn(output):
    header = "#qseqid\tsseqid\tpident\tsstart\tsend\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
    with open(output, "r+") as file:
        content = file.read()
        file.seek(0, 0)
        file.write(header + content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run tblastn and add header to the output. Do not include .db in database name")
    parser.add_argument("query", help="protein.fa")
    parser.add_argument("db", help="database.fa")
    args = parser.parse_args()

    run_tblastn(args.db + ".db", args.query)
    add_header_to_tblastn(args.query + ".tblastn")
