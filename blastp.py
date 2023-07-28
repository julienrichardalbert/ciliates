import argparse
import subprocess

def run_blastp(db, query):
    output = query + ".blastp"

    cline_blastp = [
        "blastp",
        "-db", db,
        "-query", query,
        "-out", output,
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
    args = parser.parse_args()

    run_blastp(args.db + ".db", args.query)
    add_header_to_blastp(args.query + ".blastp")
