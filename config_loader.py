# Load default paths to databases and thresholds.
# The config.ini file MUST be in the same directory as this script.

import configparser
import os
import ast

def load_config():
    # Get the absolute path of the directory containing the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    config_file_path = os.path.join(script_dir, 'config.ini')
#    print(f'Using config file: {config_file_path}')

    config = configparser.ConfigParser()
    config.read(config_file_path)

    defaults = {
        'db_prot': config.get('paths', 'db_prot'),
        'db_cds': config.get('paths', 'db_cds'),
        'db_ohno': config.get('paths', 'db_ohno'),
        'db_uniprot': config.get('paths', 'db_uniprot'),
        'blaste': float(config.get('defaults', 'blaste')),
        'preference_dictionary': ast.literal_eval(config.get('defaults', 'preference_dictionary')),
        'refname': config.get('defaults', 'refname'),
        'lenmultiplier': config.get('defaults', 'lenmultiplier'),
        'trimal_strat': config.get('defaults', 'trimal_strate'),
    }

#    print('Setting defaults')
#    for key, value in defaults.items():
#        print(f'{key}: {value}')

    return defaults

if __name__ == "__main__":
    config = load_config()
