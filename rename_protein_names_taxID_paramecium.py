'''
We must first ensure the inital file is in format:
	>speciesName_proteinName
In my case, I had some proteomes with no separator, and others with period separators.
So, I used sed to find and replace:
for x in *fa; do sed 's/\./_/g' $x > tmp && mv tmp $x; done
unfortunately, sed -i doesn't work on JRA's MBP in 2023...
and
19240  [2023-07-25 17:54:53] sed 's/TBOREA/TBOREA_/g' T.borealis.protein.fa > tmp && mv tmp T.borealis.protein.fa
19244  [2023-07-25 17:55:15] sed 's/PMMNP/PMMNP_/g' multimicronucleatum_MO3c4_annotation_v1.protein.fa > tmp && mv tmp multimicronucleatum_MO3c4_annotation_v1.protein.fa

EDIT 27 OCTOBER 2023: keep the original name of the species in the protein file. That way, you can pull out strains in the future.
'''

import argparse

def read_metadata(metadata_file):
    metadata = {}
    with open(metadata_file, "r") as file:
        for line in file:
            name, taxid = line.strip().split()
            metadata[name] = taxid
    return metadata

def replace_protein_names(proteome_file, metadata_file):
    metadata = read_metadata(metadata_file)

    with open(proteome_file, "r") as file:
        proteome_data = file.read()

    for old_ID, taxid in metadata.items():
        proteome_data = proteome_data.replace(f"{old_ID}_", f"{taxid}.{old_ID}_")

    output_file = proteome_file.replace(".fa", "_taxID.fa")
    with open(output_file, "w") as file:
        file.write(proteome_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace protein names in a proteome file using metadata.")
    parser.add_argument("proteome_file", help="BIG proteome file.")
    parser.add_argument("metadata_file", help="2-column metadata file.")
    args = parser.parse_args()

    replace_protein_names(args.proteome_file, args.metadata_file)
