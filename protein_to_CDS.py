# GET PROTEIN CODING SEQUENCES #


import argparse

# Parse fasta files because sequences can span many lines and make Julien's life difficult
# Returns a >name:sequence dictionary
def parse_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as file:
        current_sequence = ""
        current_header = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = current_sequence
                current_header = line[1:]
                current_sequence = ""
            else:
                current_sequence += line
        if current_header:
            sequences[current_header] = current_sequence
    return sequences


# Search the transcript database file for the input protein name
# In the Paramecium DB files, protein names contain P and transcripts T, so replace those first
# This doesn't screw up species names if you convert species names to TaxID first
def find_matching_transcript(protein_names, transcript_sequences):
    matching_transcripts = {}
    protein_names_with_t = [name.replace("P", "T") for name in protein_names]

    for transcript_name, transcript_sequence in transcript_sequences.items():
        if transcript_name in protein_names_with_t:
            matching_transcripts[transcript_name] = transcript_sequence

    return matching_transcripts

# Save the file
def write_output_fasta(output_file, matching_transcripts):
    with open(output_file, 'w') as file:
        for transcript_name, transcript_sequence in matching_transcripts.items():
            file.write(f">{transcript_name}\n{transcript_sequence}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract coding sequence from a protein fasta")
    parser.add_argument("protein_query_file", help="proteinA.fa")
    parser.add_argument("transcript_database_file", help="ciliates19.transcript_taxID.fa")
    args = parser.parse_args()

    output_file = args.protein_query_file + ".transcripts.fa"

    protein_names = parse_fasta(args.protein_query_file).keys()
    transcript_sequences = parse_fasta(args.transcript_database_file) # this could be sped up by only parsing the transcript db for matching transcripts, but whatever
    matching_transcripts = find_matching_transcript(protein_names, transcript_sequences)

    print(f"Number of protein sequences provided: {len(protein_names)}")
    print(f"Number of matching transcripts found: {len(matching_transcripts)}")

    write_output_fasta(output_file, matching_transcripts)
