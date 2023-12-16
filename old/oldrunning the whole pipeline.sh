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


# CREATING DATABASES

Downloaded from https://paramecium.i2bc.paris-saclay.fr/
and http://ciliates.ihb.ac.cn/database/home/#pp
and http://ciliates.ihb.ac.cn/tcgd/database/download/#vo

'''
Thankfully almost all protein annotation files had a corresponding transcript annotation file, with the exception of
pbursaria_Dd1_mac_JV3_6
…and any pbursaria strain for that matter.
Thankfully there was a gff file:
pbursaria_Dd1_mac_JV3_6_annotation_v1.0.gff3
And reference file
pbursaria_Dd1_mac_JV3_6_diploid_v1.0.fa
'''

So I used cufflinks gffread
gffread -w test2.fa -g pbursaria_Dd1_mac_JV3_6_diploid_v1.0.fa pbursaria_Dd1_mac_JV3_6_annotation_v1.0.gff3
sed 's/ .*//g' test2.fa > test3.fa
mv test3.fa pbursaria_Dd1_mac_JV3_6_annotation_v1.0.transcript.fa

then, I used
rename_protein_names_taxID.py

to rename protein names to include TaxIDs (used for downstream analyses).
# make sure each ID was converted properly:
grep -e ">" ciliates19.transcript_taxID.fa | cut -f1 -d '.' | uniq -c

# index the database example
makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot.fasta.db -dbtype prot

# NEW let's try using a nucleotide database instead
makeblastdb -in paramericum_27oct2023.fasta -out paramericum_27oct2023.fasta.db -dbtype nucl


# DATABASE 1 (19 ciliates)
# grep -e ">" ./proteomes/ciliates_taxID.fa | cut -f1 -d '.' | sort | uniq -c
35534 >219701
33793 >224956
13186 >266149
40810 >304693
41085 >304694
44396 >43137
41676 >43138
37098 >44029
17834 >44030
18673 >5885
42615 >5886
40460 >5888
26996 >5911
40261 >65126
36094 >65128
49951 >65129
36179 >65130
31192 >74790
20694 >78437
{5885: 'Paramecium caudatum',
5886: 'Paramecium primaurelia',
5888: 'Paramecium tetraurelia',
5911: 'Tetrahymena thermophila',
43137: 'Paramecium octaurelia',
43138: 'Paramecium pentaurelia',
44029: 'Paramecium jenningsi',
44030: 'Paramecium multimicronucleatum',
65126: 'Paramecium biaurelia',
65128: 'Paramecium sexaurelia',
65129: 'Paramecium sonneborni',
65130: 'Paramecium tredecaurelia',
74790: 'Paramecium bursaria',
78437: 'Taonius borealis',
219701: 'Paramecium novaurelia',
224956: 'Paramecium quadecaurelia',
266149: 'Pseudocohnilembus persalinus',
304693: 'Paramecium decaurelia',
304694: 'Paramecium dodecaurelia'}
>PPERSA	266149
>TBOREA	78437
>TTHERM	5911
>PMMNT	44030
>PBIA_V1_4_1	65126
>PBUR_Dd1_1	74790
>PCAU_43c3d_1	5885
>PDEC_223_1	304693
>PDODEC_274_1	304694
>PJENN_M_1	44029
>PNOV_TE_1	219701
>POCTA_138_1	43137
>PPENT_87_1	43138
>PPRIM_AZ9-3_1	5886
>PQUADEC_NiA_1	224956
>PSEX_AZ8_4_1	65128
>PSON_ATCC_30995_1	65129
>PTET_51_1	5888
>PTRED_209_2	65130

./transcriptomes/ciliates19.transcript_taxID.fa
grep -e ">" ./transcriptomes/ciliates19.transcript_taxID.fa | cut -f1 -d '.' | sort | uniq -c | wc -l
      19





# START WITH CENH3 IN PARAMECIUM
Starting with an aa sequence provided by Sandra:
>CenH3_1_tetraurelia
MANKKTTKENNNQSFQVDNNEKMPSFMHSLFDSDEKSVIAEKSNERSKKS
EKKKVERISAIQVQKARDKLQKRNKPMSKVLQEIRQLQASSVLVCRRAGF
QRFVRQTGIKVSDELGFKEFRYSSKSLECLQTLTEQYMVDLFEDSVQCTF
HAKRVTLMAKDLNLTARIRGIEQPLQEIRNLKLR

# it would be great to be able to start with a list of proteins that all belong to the same gene family. That way, we can compare positive and negative selection in a focused context.
#

# start with a protein with the following nomenclature:


# identify annotated ohnologs and output a tree

# identified potential homologues using blastp
python ~/github/ciliates/blastp.py Ptetraurelia_CenH3_1.fa ~/Dropbox/ciliates/proteomes/ciliates_taxID.fa

# extracted the amino acid sequences of hits
python ~/github/ciliates/extract_sequence_blastp_hits.py Ptetraurelia_CenH3_1.fa.blastp ./proteomes/ciliates_taxID.fa

# align amino acids and create a newick tree
python ~/github/ciliates/fasta_align_tree.py Ptetraurelia_CenH3_1.fa.blastp.fa

# created annotation files for tree: blast
python ~/github/ciliates/blast_annotation.py Ptetraurelia_CenH3_1.fa.blastp.fa

# created annotation files for tree: taxon
python ~/github/ciliates/Figtree-Analysis/ncbi_taxid_to_annotation.py Ptetraurelia_CenH3_1.fa.blastp.fa

# created annotation files for tree: pfam (I did not test this properly)
python ~/github/ciliates/Figtree-Analysis/pfam_annotation.py Ptetraurelia_CenH3_1.fa.blastp.fa.extract Ptetraurelia_CenH3_1.fa.blastp.fa.extract.pfam

# open tree in FigTree
# Selection Mode: Taxa
# Import annotations and select the cluster of cenH3 proteins and change to blue
# save the tree (File > Save As...)
# extract only the specific cluster of coloured cenH3 (blue)
python ~/github/ciliates/Figtree-Analysis/remove_and_extract.py Ptetraurelia_CenH3_1.fa.blastp.fa Ptetraurelia_CenH3_1.fa.blastp.fa.aligned.coloured.nwk NA 0000ff

# repeat the tree construction step to see your work
python ~/github/ciliates/fasta_align_tree.py Ptetraurelia_CenH3_1.fa.blastp.fa.extract
python ~/github/ciliates/blast_annotation.py Ptetraurelia_CenH3_1.fa.blastp.fa.extract
python ~/github/ciliates/Figtree-Analysis/ncbi_taxid_to_annotation.py Ptetraurelia_CenH3_1.fa.blastp.fa.extract
# and then open the newick file in FigTree

# Branch point
#
# Branch 1: if not interested in dN/dS, just keep the aas, align those
#
# Branch 2: if interested in testing evolutionary models, you need the coding sequences to calculate dN/dS and do modeling
python ~/github/ciliates/protein_to_CDS.py Ptetraurelia_CenH3_1.fa.blastp.fa.extract transcriptomes/ciliates19.transcript_taxID.fa

# align the codiing sequences using the amino acid sequences to conserve codons
# I can make this step better
# I had to edit the transcript names because they have to match with protein names and now they differ  by P and T. Should just remove all Ps and Ts!!!
 sed -i .dak 's/T/P/' Ptetraurelia_CenH3_1.fa.blastp.fa.extract.transcripts.fa

# set threshold to 0.0 to force conversion to nucleotides
# http://etetoolkit.org/cookbook/ete_build_mixed_types.ipynb
ete3 build -a Ptetraurelia_CenH3_1.fa.blastp.fa.extract.renamed -n Ptetraurelia_CenH3_1.fa.blastp.fa.extract.transcripts.fa -o mixed_types/ -w standard_fasttree --clearall --nt-switch-threshold 0.0



I NEED HELP WITH ALIGNMENTS. I AM NOT DOING IT CORRECTLY~!!!!

then I fucked around in test_evolutionary_models.ipynb to get ete3 to cooperate. Its been a pain but I think it works now. Will have to code a python command for this step.



manually downloaded accession:
A. thaliana	NM_001035850
O. pumila	AY612788
C. bursa pastoris	AY612789
C. flexuosa	AB299174
A. hirsuta	AB299170
E. sativa	AB299180
R. sativus	AB299184
B. napus	GU166737
B. oleracea	GU166739
V. vinifera	XM_002281037


for x in *fa; do sed 's/\ /_/g' $x > $x".tmp"; done
cat *tmp > rosia_cenH3_10.fa
python ~/github/ciliates/translate.py -i rosia_cenH3_10.fa

ete3 build -a rosia_cenH3_10.fa.aa -n rosia_cenH3.fa -o rosia_10_mixed_types/ -w standard_fasttree --clearall --nt-switch-threshold 0.0






# I am still missing the superposition of known protein domains with the tree
Repeat with PRDM9
Repeat with Eno1
Repeat with canonical H3 or H4
