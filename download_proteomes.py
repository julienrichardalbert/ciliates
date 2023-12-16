import argparse
import requests
from datetime import datetime
import shutil
from ete3 import NCBITaxa
ncbi = NCBITaxa()


def get_species_name(taxid):
    latin = ncbi.get_taxid_translator([int(taxid)])[int(taxid)].replace(" ", "_")
#    print(latin)
    today = datetime.now().strftime("%Y-%m-%d")
    output_file = f"{latin}_{taxid}_uniprotkb_{today}.fasta"
    return output_file


def download_proteome(taxid, output):
    print(f"Downloading {taxid}...")
    url = f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28%28taxonomy_id%3A{taxid}%29%29'
    with requests.get(url, stream=True) as request:
        request.raw.decode_content = True # https://github.com/psf/requests/issues/2155
        with open(output, 'wb') as f:
            shutil.copyfileobj(request.raw, f) # John Zwinck answer to https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
    return output


def main():
    parser = argparse.ArgumentParser(description="Download proteome data from a list of TaxIDs.")
    parser.add_argument("input_file", help="Tab-separated file with taxids in the first column, or a single TaxID integer.")
    args = parser.parse_args()

    # if a taxid is input
    if args.input_file.isdigit():
        taxid = args.input_file
        output_file = get_species_name(taxid)
        download_proteome(taxid, output_file)

    # if you input a file with taxids in column 1
    else:
        with open(args.input_file, 'r') as infile:
            for line in infile:
                taxid = line.strip().split("\t")[0]
                output_file = get_species_name(taxid)
                download_proteome(taxid, output_file)

if __name__ == "__main__":
    main()
