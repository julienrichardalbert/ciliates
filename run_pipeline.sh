# notes on phylogenetics of Paramecium cenH3

# INSTALLATION
'''
I dont remember exactly how I got ete3 itself to work, but I managed to first install it without any external apps.

Maybe like this?:
CONDA_SUBDIR=osx-64 conda create -n myenv_x86 python=3.9
conda activate myenv_x86
conda config --env --set subdir osx-64

no... I had to specifically install python v3.6, so this could not have been the correct command!

Then, I used
ete3 upgrade-external-tools
to install the external tools. Running a command like ete3 build resulted in errors that pointed towards where the external apps were being looked for.

It was a matter of finding the external app:
which clustalo
and making a link to the conda env directory.



'''

# creating databases of transcriptomes
# see old pipeline doc

# index database using NCBI (only run once, even though its fast)
makeblastdb -in 17_paramecia_transcriptomes_taxID.fa -out 17_paramecia_transcriptomes_taxID.fa.db -dbtype nucl

# start with a protein sequence. The nomenclature is very important
# >TaxID.Species_strain_transcriptID
# you probably have to convert Protein (P) or Gene (G) to Transcript (T)
# the script is not smart enough to do this automatically
'''
>5888.PTET_51_T0460249
MGKKETKDKKPTQFKKIYELLSKYTQVIIVG
...
'''

# it would be great to be able to start with a list of proteins that all belong to the same gene family. That way, we can compare positive and negative selection in a focused context.
#

# identify annotated ohnologs and output a tree
python ~/github/ciliates/find_ohnologues.py P0460249_riboprotein.fa ../references/wgd_annotations/17_ciliates_wgd.txt ../references/transcriptomes/17_paramecia_transcriptomes_taxID.fa

# if no ohnologs found, skip to the tblastn stage


# trim the alignment and compare trees
# THIS STEP IS NOT REALLY NECESSARY. SKIP.
python ~/github/ciliates/fix_alignments.py P0460249_riboprotein_ohnologs.fa.aligned

'''
Now we know we are dealing with a gene that has 3 ohnologs; WGD2 (x2) and WGD1
'''

# identify potential homologues using tblastn
python ~/github/ciliates/tblastn.py P0460249_riboprotein.fa ../references/transcriptomes/17_paramecia_transcriptomes_taxID.fa

# extracted the amino coding sequence of hits
python ~/github/ciliates/extract_sequence_blast_hits.py P0460249_riboprotein.fa.tblastn ../references/transcriptomes/17_paramecia_transcriptomes_taxID.fa

# align coding sequence and create a newick tree
python ~/github/ciliates/fasta_align_tree.py P0460249_riboprotein.fa.tblastn.fa

# trim alignment and compare trees (I'm still trying to get a hang of how to make good trees, duplicating effort)
python ~/github/ciliates/fix_alignments.py P0460249_riboprotein.fa.tblastn.fa.aligned

'''
Comparing trees, there is a big difference using trimAl to make trees better.
Is this reflected in the html output? YES!
For example, tredecaurelia has this super long isoform in the alignment that blows out non-trimmed trees, and this is clear in the html output
Maybe I should ALWAYS trim...
'''

'''
from the trimmed output, there are clearly 4 clusters. There are 4 ohnologs. Test ohnolog annotation script on trimmed alignment.
'''

# Annotate the tree
python ~/github/ciliates/annotate_using_wgd.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1 ../references/wgd_annotations/17_ciliates_wgd.txt
python ~/github/ciliates/Figtree-Analysis/ncbi_taxid_to_annotation.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1
python ~/github/ciliates/annotate_using_blastx.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1

'''
It is super helpful to load the original ohnologue annotations from the start into this tree!!! (fullseq)
'''

# open tree in FigTree
# Selection Mode: Taxa
# Import annotations
# save the tree (File > Save As...)
# extract only the specific blue cluster
# JRA added the tree creation step to this
python ~/github/ciliates/remove_and_extract.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1 P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.coloured.nwk NA 0000ff

# annotate the new tree
python ~/github/ciliates/annotate_using_taxid.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract
python ~/github/ciliates/annotate_using_blastx.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract
python ~/github/ciliates/annotate_using_wgd.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract ../references/wgd_annotations/17_ciliates_wgd.txt

# and then open the newick file in FigTree and repeat
python ~/github/ciliates/remove_and_extract.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1 P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract.coloured.nwk NA 0000ff
python ~/github/ciliates/annotate_using_taxid.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract
python ~/github/ciliates/annotate_using_blastx.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract
python ~/github/ciliates/annotate_using_wgd.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract ../references/wgd_annotations/17_ciliates_wgd.txt

#GOOD. No more ohnologues. Begin EvolTree testing!!
python python ~/github/ciliates/evolution.py P0460249_riboprotein.fa.tblastn.fa.trimal.automated1.extract

# Done!




names_tmp=$(cat PTET_51_solo_genes.txt)
for name_tmp in $names_tmp; do
    echo $name_tmp > tmp
    name=$(cut -f1 tmp | cut -f1 -d '-')
    echo $name
done | head
