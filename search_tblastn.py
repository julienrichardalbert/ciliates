# Starting with a protein sequence, find all genomic loci that encode it

import os
import argparse
import subprocess
import pandas as pd
import pybedtools

def run_tblastn(input_file, db, output_file, eval):
    cline_tblastn = [
        'tblastn',
        '-db', db,
        '-query', input_file,
        '-out', output_file,
        '-evalue', eval,
        '-db_gencode', '6',
        '-outfmt', '6 qseqid sseqid pident sstart send length mismatch gapopen qstart qend sstart send evalue bitscore'
    ]
    subprocess.run(cline_tblastn)


def add_header_to_tblastn(output):
    header = "#qseqid\tsseqid\tpident\tsstart\tsend\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
    with open(output, "r+") as file:
        content = file.read()
        file.seek(0, 0)
        file.write(header + content)

def tblastn_to_fasta_no_introns(input_file, db, output_file):
    df = pd.read_csv(input_file, sep='\t')
    df['species'] = df['sseqid'].apply(lambda x: x.rsplit('_', 1)[0])

    # group by species and get the row with the lowest 'evalue'
    result_df = df.sort_values(['species', 'evalue'], ascending=[True, True]).groupby('species').head(1)

    if not result_df.empty:
        #print(result_df)

        # .iloc[0] gets the value without the index
        out_chr = str(result_df['sseqid'].iloc[0])
        out_start = int(result_df['sstart'].iloc[0])
        out_end = int(result_df['send'].iloc[0])

        if out_start < out_end:
            out_strand = '+'
            out_start = out_start - 1 # BED is 0-based. Ensures CDS starts with ATG
        else:
            out_strand = '-'
            out_start, out_end = out_end, out_start  # swap values
            out_start = out_start - 1 # CDS ends with TGA



        bed_data = {
            'chromosome': [out_chr],
            'start': [out_start],
            'end': [out_end],
            'name': ['.'],
            'score': ['0'],
            'strand': [out_strand],
        }
        df = pd.DataFrame(bed_data)

        # Create a BedTool from the DataFrame
        bed_df = pybedtools.BedTool.from_dataframe(df)
        seq = bed_df.sequence(fi=db, s=True, fo=output_file)
    else:
        print(f'No tBLASTn hits for genome: {db}')




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run a single instance of a combination of steps, which perform the steps: blastp | align | build tree | annotate | trim | cds | evolution')
    parser.add_argument('-s', '--start_file', required=True, help='Protein FASTA file.')
    parser.add_argument('-db_genomes', '--db_genomes',  default='/Users/jra/Dropbox/ciliates/references/genomes/', help='full path to protein database(s).')
    parser.add_argument('-blaste', '--blaste', type=float, default=0.1, help='BLASTp evalue to threshold results.')
    args = parser.parse_args()


    # Loop over reference genomes in the directory
    for database in os.listdir(args.db_genomes):
        if database.endswith('.fa'):
            database = database + '.db'
            print(f'BLAST against {database}')
            # Construct the full path to the database file
            database_path = os.path.join(args.db_genomes, database)
            output_file = f"{database.replace('.fa.db', '')}.tblastn.txt"
            run_tblastn(args.start_file, database_path, output_file, str(args.blaste))
            add_header_to_tblastn(output_file)
            tblastn_to_fasta_no_introns(output_file, database_path.replace('.fa.db', '.fa'), output_file + '.fa')
            # combine all files, make small db
            # translate to protein
            # trimal, backtranslate using small db
            # evolution
