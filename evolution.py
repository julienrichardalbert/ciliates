# From  an input alignemnt and tree, get the amino acid positions under positive selection
#
# following my other script nomenclature, files must follow nomenclature
# alignment: sequences.extract.trimal.aligned
# tree:      sequences.extract.trimal.aligned.nwk

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
from ete3 import NCBITaxa
import matplotlib.pyplot as plt
import os
import shutil
import argparse
from log_progress import log

ncbi = NCBITaxa()

def file_exists(filename):
    current_directory = os.getcwd()
    if os.path.exists(filename):
        print(f"The file {filename} exists in the current working directory {current_directory}.")
    else:
        print(f"The file {filename} does not exist in the current working directory {current_directory}.")

def organize_your_life(input_alignment):
    # get the full path of the input file
    full_input_file_path = os.path.abspath(input_alignment)
    full_input_newick_path = full_input_file_path + '.nwk'

    directory, alignment_file_name = os.path.split(full_input_file_path)
    directory, newick = os.path.split(full_input_newick_path)

    path = os.path.join(directory, f"{alignment_file_name}.evolTree")
    os.makedirs(path, exist_ok=True)

    # copy input files the new directory
    shutil.copy(full_input_file_path, path)
    shutil.copy(full_input_newick_path, path)

    # move to the new directory
    os.chdir(path)
    print(f"EvolTree running in new folder: {path}")

    # load the tree
    tree = EvolTree(newick = full_input_newick_path, format=1)
    # associate alignment to tree
    tree.link_to_alignment(input_alignment)
    tree.workdir = path # set working directory where models will be saved
    alnPrefix = alignment_file_name # might need this in the future, like for saving files?
    return tree, alnPrefix

def organize_your_life2(input_alignment):
    # get the full path of the input file
    tree = input_alignment + '.nwk'
    # load the tree
    tree = EvolTree(newick = tree, format=1)
    # associate alignment to tree
    tree.link_to_alignment(input_alignment)
    current_directory = os.getcwd()
    tree.workdir = current_directory # set working directory where models will be saved
    alnPrefix = input_alignment # might need this in the future, like for saving files?
    return tree, alnPrefix

# Get species name from TaxID
def get_species_name(taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        species_name = names[taxid]
        return species_name
    except:
        return f"Unknown_species_{taxid}"

# calculate dS values relative to a reference leaf
def calculate_dS(prefix, tree, reference):
    partial_leaf_name = reference
    output_txt = prefix + '.ref.dS.txt'
    output_fig  = prefix + '.ref.dS.pdf'
    # Initialize variables to store the reference leaf and its full name
    reference_leaf = None
    reference_leaf_name = None

    # Search for the leaf with a name that starts with the specified partial name
    for leaf in tree:
        if leaf.name.startswith(partial_leaf_name):
            reference_leaf = leaf
            reference_leaf_name = leaf.name

    if reference_leaf:
        # Calculate dS values for leaves relative to the reference leaf
        dS_values = {}

        for leaf in tree:
            # Extract TaxID from the leaf name
            taxid = int(leaf.name.split('.')[0])

            # Get the scientific name associated with the TaxID
            species_name = ncbi.get_taxid_translator([taxid]).get(taxid, leaf.name)

            # Calculate the dS value for the leaf relative to the reference leaf
            dS = reference_leaf.get_distance(leaf)
            dS_values[leaf.name] = {'Species': species_name, 'dS': dS}

        with open(output_txt, "w") as f:
            f.write('#gene\tspecies\tdS\n')
            for leaf_name, data in dS_values.items():
                #print(f"{leaf_name}\t{data['Species']}\t{data['dS']:.4f}")
                f.write(f"{leaf_name}\t{data['Species']}\t{data['dS']:.4f}\n")
            log(f"Ds values have been written to {output_txt}")

        # Create a bar chart
        geneName = str(prefix.split("_")[0])
        leaf_names = list(dS_values.keys())
        species_names = [data['Species'] for data in dS_values.values()]
        dS_values_list = [data['dS'] for data in dS_values.values()]

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(species_names, dS_values_list, color='black')
        ax.set_xlabel('dS Values')
        ax.set_title('dS Values relative to reference leaf: ' + partial_leaf_name + ' for gene ' + geneName)

        # Set x-axis limit to 1.1 if dS values are less than 1.1
        max_dS = max(dS_values_list)
        if max_dS < 1.1:
            ax.set_xlim(0, 1.1)

        # Add a red dotted line at x=1
        ax.axvline(x=1, color='red', linestyle='--')

        # Save the plot
        plt.savefig(output_fig)

    else:
        log(f"No leaf found in the tree with a name starting with '{partial_leaf_name}'.")


# run the evolutionary models (longest step)
# intelligently checking whether a model exists was a massive waste of time because ete3 is a fucking fuck and won't "find" the models after it is relaunched
def run_models(prefix, tree, models):
    #def folder_exists(folderPath):
    #    return os.path.exists(folderPath) and os.path.isdir(folderPath)
    #model = ["M1", "M2", "M7", "M8", "M0", "fb"]
    #for model in models:
    #    if not folder_exists(model):
    #        tree.run_model(model)
    #    else:
    #        print(f"The model {model} exists. Skipping")
    # TABARNAK WHY THE FUCK CANT THIS STUPID FUCK ETE3 FIND THE FUCKING MODELS IN THE FUCKING WORKING DIRECTORY I FUCKING SPECIFIED
    log(f'Computing model set: {models}')
    for model in models:
        log(f'Computing model: {model}')
        tree.run_model(model)

    output_models = prefix + '.models.txt'
    with open(output_models, "w") as output:
        for model in models:
            #print(model)
            get_model = tree.get_evol_model(model)
            #print(get_model)
            output.write(str(get_model)+'\n')

def get_pvals(prefix, tree, alt, neg, pnum):

    altModel = tree.get_evol_model(alt)
    pval = tree.get_most_likely(alt, neg)
    log(f'{alt} vs {neg} pval: {pval}')
    output_file = prefix + '.pvals.txt'
    with open(output_file, "a") as output:
        output.write('#altModel\tnegModel\tpval\n')
        output.write(alt + '\t' + neg  +  '\t' + str(pval) + '\n')

    if pval < 0.05:
        output_file2 = prefix + '.sigAAs.txt'
        with open(output_file2, "a") as output2:
            output2.write('#alt_model\tnegmodel\taa\tposition\tprobability\tomega\n')
            log(f'{alt} model wins.')
            for s in range(len(altModel.sites['BEB']['aa'])):
                    if altModel.sites['BEB'][pnum][s] > 0.95:
                        pos_acid  = altModel.sites['BEB']['aa'][s]
                        pos_site  = str(s+1)
                        pos_prob  = str(altModel.sites['BEB'][pnum][s])
                        pos_omega = str(altModel.sites['BEB']['w'][s])

                        output_line2 = alt + '\t' + neg + '\t' + pos_acid + '\t' + pos_site + '\t' + pos_prob + '\t' + pos_omega + '\n'
                        log(output_line2)
                        output2.write(output_line2)
    else:
        log(f'{neg} model is not rejected')

def evol_graphs(prefix, tree, alt, neg, suffix):
    # Change the colours
    col_line = {'NS' : '#777777',
           'RX' : '#777777', 'RX+': '#777777',
           'CN' : '#777777', 'CN+': '#777777',
           'PS' : '#777777', 'PS+': '#777777'}
    col_bar = {'NS' : '#BCBCBC',
           'RX' : '#5D63AB', 'RX+': '#5D63AB',
           'CN' : '#659A62', 'CN+': '#659A62',
           'PS' : '#F4C95D', 'PS+': '#F4C95D'}

    # Print the first tree without naming branches
    alt_model = tree.get_evol_model(alt)
    neg_model = tree.get_evol_model(neg)

    alt_model.set_histface(up=True, colors=col_line, errors=True, kind='curve', ylim=[0,20], hlines = [1], hlines_col=['black'])
    tree.render(prefix +  '.' + alt + suffix + '.line.pdf', histfaces=[alt])
    alt_model.set_histface(up=True, colors=col_line, errors=True, kind='curve', ylim=[0,10], hlines = [1], hlines_col=['black'])
    tree.render(prefix +  '.' + alt + suffix + '.line2.pdf', histfaces=[alt])

    alt_model.set_histface(up=True, colors=col_bar, errors=True, kind='stick', ylim=[0,4], hlines = [1], hlines_col=['black'])
    tree.render(prefix +  '.' + alt + suffix + '.bar.pdf',  histfaces=[alt])
    alt_model.set_histface(up=True, colors=col_line, errors=False, kind='curve', ylim=[0,4])
    tree.render(prefix +  '.' + alt + suffix + '.barl.pdf', histfaces=[alt])

def modify_leaf_names_reroot(tree, prefDic):
    # Modify the tree leaf names
    # ONLY DO THIS AS THE LAST STEP
    for leaf in tree.iter_leaves():
        taxid = int(leaf.name.split('.')[0])
        species_name = get_species_name(taxid)
        leaf.name = species_name.split(' ')[1] # discard genus, setting outgroup takes only first 3 characters
    root_point = tree.get_farthest_oldest_leaf(prefDic)
    tree.set_outgroup(root_point)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Measure positive selection using ete3.")
    parser.add_argument("alignment_file", help="Coding sequence multiple sequence alignment")
    parser.add_argument("-r", "--reference_leaf", help="Reference leaf name", default='5888.PTET_51')
    args = parser.parse_args()

    tree, alnPrefix = organize_your_life(args.alignment_file)
    print(tree) # because its fun to do
    calculate_dS(alnPrefix, tree, args.reference_leaf) # partial name of the reference leaf
    models = ['M1', 'M2', 'M7', 'M8']
    run_models(alnPrefix, tree, models)
    get_pvals(alnPrefix, tree, 'M2', 'M1', 'p2') # alt, neg
    get_pvals(alnPrefix, tree, 'M8', 'M7', 'p10')
    evol_graphs(alnPrefix, tree, 'M2', 'M1', '') # alt, neg, suffix
    evol_graphs(alnPrefix, tree, 'M8', 'M7', '')
    preference_dictionary = {'tet':1,
        'bia':2, 'dec':2, 'dod':2, 'jen':2, 'nov':2, 'oct':2,'pen':2, 'pri':2, 'qua':2, 'sex':2, 'son':2, 'tre': 2,
                            'cau':3}
    modify_leaf_names_reroot(tree, preference_dictionary) # change leaf name nomenclature from TaxID.gene to spe and reroot
    evol_graphs(alnPrefix, tree, 'M2', 'M1', '_sp')
    evol_graphs(alnPrefix, tree, 'M8', 'M7', '_sp')
