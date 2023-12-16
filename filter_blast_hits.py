# filter a .blast output file for the top hit of each species
# Species names must be in format: TaxID.Species_name_can_be_long_TranscriptID

import argparse
import pandas as pd
from log_progress import log


def filter_best_hits(input_file, output_file, args):

    log('Filtering blastp hits for the top per species')
    df = pd.read_csv(input_file, sep='\t')

    # Extract species from 'sseqid' column
    df['species'] = df['sseqid'].apply(lambda x: x.rsplit('_', 1)[0])

    # Group by species and get the row with the highest 'pident' and lowest 'evalue'
    result_df = df.sort_values(['species', 'evalue'], ascending=[True, True]).groupby('species').head(1)

    # Save the result to the output file
    result_df.to_csv(output_file, sep='\t', index=False)

    log(f'{input_file}: {len(df)}')
    log(f'{output_file}: {len(result_df)}')

    return len(result_df) > 1 # returns True if there is more than 1 hit.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter top blast hits, returning one per species.")
    parser.add_argument("hits", help="Input file containing blast hits data.")
    args = parser.parse_args()
    output = args.hits + '.top'


    filter_best_hits(args.hits, output, args)
