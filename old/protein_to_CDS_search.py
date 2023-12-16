import re
import argparse
from Bio import Entrez, SeqIO

Entrez.email = "jrichardalbert@gmail.com"

def extract_info_from_protein_name(protein_name):
    # Use regular expressions to extract accession, taxid, and gene name
    accession_match = re.search(r'\|([^|]+)\|', protein_name)
    taxid_match = re.search(r'OX=(\d+)', protein_name)
    gene_name_match = re.search(r'GN=([^ ]+)', protein_name)

    if accession_match and taxid_match and gene_name_match:
        accession = accession_match.group(1)
        taxid = taxid_match.group(1)
        gene_name = gene_name_match.group(1)
        return accession, taxid, gene_name
    else:
        return None, None, None

def fetch_and_save_cds_sequence(protein_name, accession, taxid, gene_name, output_file):
    try:
        # Search for the CDS using the derived accession, taxid, and gene name
        search_term = f"{accession} OR {gene_name} AND txid{taxid}[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
        record = Entrez.read(handle)

        if "IdList" in record and record["IdList"]:
            # Get the CDS record
            print(f"Searching for: {protein_name}")
            cds_id = record["IdList"][0]
            cds_handle = Entrez.efetch(db="nucleotide", id=cds_id, rettype="gb", retmode="text")
            cds_record = SeqIO.read(cds_handle, "genbank")
            description = cds_record.description
            print(f"Found match - CDS ID: {cds_id}: {description}" )

            # Write the original protein sequence and CDS sequence to the output file
            with open(output_file, "a") as outfile:
                outfile.write(f"{protein_name}\n{cds_record.seq}\n")

        else:
            print(f"No CDS found for {protein_name}")

    except Exception as e:
        print(f"Error processing {protein_name}")
        print(e)

def main():
    parser = argparse.ArgumentParser(description="From proteins get CDS from NCBI.")
    parser.add_argument("input_file", help="File containing protein sequences in FASTA format.")
    parser.add_argument("output_file", help="File to save CDS sequences.")
    args = parser.parse_args()

    with open(args.input_file, "r") as infile:
        for line in infile:
            if line.startswith(">"): # it doesn't have to, as the regex will return nothing if this processes a sequence line, but good to do anyway
                protein_name = line.strip()
                accession, taxid, gene_name = extract_info_from_protein_name(protein_name)
                if accession and taxid and gene_name:
                    fetch_and_save_cds_sequence(protein_name, accession, taxid, gene_name, args.output_file)

if __name__ == "__main__":
    main()
