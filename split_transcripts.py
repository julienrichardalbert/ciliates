import argparse
import os

def get_transcript_name(transcript_header):
    # Split the header by spaces and return the first part
    return transcript_header.split()[0].split('.')[-1]

def split_transcripts(input_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    transcripts = []
    current_transcript = {"name": "", "sequence": ""}
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_transcript["name"]:
                transcripts.append(current_transcript)
                current_transcript = {"name": "", "sequence": ""}
            current_transcript["name"] = get_transcript_name(line)
        else:
            current_transcript["sequence"] += line
    
    if current_transcript["name"]:
        transcripts.append(current_transcript)
    
    return transcripts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split the input file into individual transcripts")
    parser.add_argument("input_file", help="Input file containing transcript annotations")
    args = parser.parse_args()

    # Create a folder with the input file's base name
    output_folder = os.path.splitext(os.path.basename(args.input_file))[0]
    os.makedirs(output_folder, exist_ok=True)

    transcripts = split_transcripts(args.input_file)
    
    for i, transcript in enumerate(transcripts, start=1):
        # Use the transcript name directly (already cleaned)
        output_file = os.path.join(output_folder, f"{transcript['name']}.txt")
        with open(output_file, "w") as out:
            out.write(">" + transcript["name"] + "\n")
            out.write(transcript["sequence"])
        print(f"Transcript {i} saved to {output_file}")
