**************************************************************************************************************
# Protein Embedding Based Alignments (PEbA)
**************************************************************************************************************

Paper can be viewed at https://www.authorea.com/users/623259/articles/646069-protein-embedding-based-alignment.

This project uses Rostlab's ProtT5-XL-UniRef50 encoder (https://huggingface.co/Rostlab/prot_t5_xl_uniref50) to
embed protein sequences and then calculates the cosine similarity between each amino acid's vector to align them 
using traditional sequence alignment algorithms. The cosine similarity scores are used in place of substitution
matrices, such as BLOSUM, with the goal of producing accurate alignments for sequences that share low character
identity (<20-35%), referred to as the "twilight zone" of protein sequence alignment.

PEbA is implemented in peba.py, where it can use either the Needleman-Wunsch 'global' algorithm or the Smith-Waterman
'local' algorithm to align two sequences. It is designed to be run from the command line. The following arguments
are allowed, with -embed1 and -embed2 being optional arguments if the corresponding sequences from -file1 and -file2
have already been been embedded and saved as a numpy array:

    -a,  --align <str>          : Alignment algorithm to use (global or local) (default: local)
    -f1, --file1 <file1.fa>     : Path to first sequence fasta file
    -f2, --file2 <file2.fa>     : Path to second sequence fasta file
    -e1, --embed1 <file1.txt>   : Path to first sequence embedding file
    -e2, --embed2 <file2.txt>   : Path to second sequence embedding file
    -go, --gopen <int/float>    : Gap opening penalty (default: -11)
    -ge, --gext <int/float>     : Gap extension penalty (default: -1)
    -e, --encoder <str>         : Encoder used to generate embeddings (default: ProtT5)
    -o, --output <str>          : Desired output format, either 'msf' or 'fa' (default: msf)
    -s, --savefile <str>        : Path to save alignment file, prints to console if not specified

Embeddings from any model can be used as long as the embeddings are a 1D array the same length as the sequence.
peba.py loads embeddings using numpy.loadtext(), so if you are pre-generating your own embeddings make sure they
are compatible with said numpy function. The -encoder argument is used to specify the model used and written to
the output for reference.

**************************************************************************************************************
# Installing Requirements
**************************************************************************************************************

To install PEbA, clone this repository and install the necessary requirements.

```
git clone https://github.com/mgtools/PEbA.git
cd PEbA
pip install -r base_requirements.txt
```

This project was developed using python 3.10.6. Necessary requirements for running PEbA can be installed from
requirements.txt.

**************************************************************************************************************
# Running PEbA
**************************************************************************************************************

An example of running PEbA with two fasta files and their corresponding embeddings is shown below:

```
python scripts/peba.py -f1 data/example/1j46_A.fa -f2 data/example/1k99_A.fa -e1 data/example/1j46_A.txt -e2 data/example/1k99_A.txt
```

On default settings, this will produce the following output in the console:

```
PileUp



   MSF: 93  Type: P  Method: ProtT5_Sim  Gopen: -11.0  Gext: -1.0

 Name: 1j46_A oo  Len:  93  Start/End:  0,85
 Name: 1k99_A oo  Len:  93  Start/End:  6,91

//



1j46_A      MQDRVKRPMN AFIVWSRDQR RKMALENPRM RNSEISKQLG YQWKMLTEAE 
1k99_A      HPDFPKKPLT PYFRFFMEKR AKYAKLHPEM SNLDLTKILS KKYKELPEKK 

1j46_A      KWPFFQEAQK LQAMHREKYP NYKYRPRRKA KMLPK
1k99_A      KMKYIQDFQR EKQEFERNLA RFREDHPDLI QNAKK
```

To embed sequences, PEbA will use the ProtT5-XL-UniRef50 model and T5encoder from HuggingFace. If -e1 and
-e2 are specified, meaning that the embeddings are already generated, these models will not be used.

**************************************************************************************************************
# Comparing PEbA to Reference Alignments
**************************************************************************************************************

To determine if PEbA alignments produce more accurate alignments than those with subsitution matrices, both 
types of alignments are compared to reference alignments from BAliBASE3 (https://www.lbgi.fr/balibase/), found
in the data/BAliBASE_R1-5 folder. BAliBASE alignments are structure based benchmarks used to compare new methods of
alignment, such as this one. References RV11 and 911 contain sequences with <20% sequence identity, references
RV12 and RV912 with 20-40% sequence identity, and reference RV913 with 40-80% sequence identity.

The 'Sum-of-Pairs' (SP) score is used to compare the performance between the alignments generated by this
project to the reference alignments. This percentage indicates how many residue pairs are shared between two
alignments. compute_score.py computes this metric between two alignments, the first being the reference alignment,
and the second being the alignment being compared. An example of the output is shown below:


SP: 0.551   ref_length: 690   comparison_length: 356   similarity: 0.146


SP score is a proportion between 0 and 1, 0 representing no shared residue pairs between two alignments, and 1
representing that all residue pairs are found in both alignments. The SP score is calculated by dividing the
number of shared residue pairs between the reference and test alignment by the number of pairs in the reference
alignment. ref_length is the length of the reference alignment, and comparison_length is the number of residue
pairs being compared (gaps are not counted). The similarity is the percentage of identical residues between
the two sequences that are matched together i.e. if a pair has the same two residues, it is a match.

F1 score can also be calculated to compare alignments, but SP score is the default.

**************************************************************************************************************
# Data and Methods Compared
**************************************************************************************************************

From each BAliBASE benchmark MSA we extracted each sequence and each pairwise alignment. We took each pair of
sequences and generated alignments using PEbA with ProtT5 embeddings, PEbA with ESM2 embeddings, BLOSUM, DEDAL,
vcMSA, and FATCAT (only on RV11).

**************************************************************************************************************
# Language Models and Other Methods
**************************************************************************************************************

The following guide was used for working with ProtT5:
https://github.com/agemagician/ProtTrans/blob/master/Embedding/PyTorch/Advanced/ProtT5-XL-UniRef50.ipynb
We adapted this code and used it to embed sequences in embed_seqs.py.

The readme for ESM2 can be found here, which contains instructions for using the model:
https://github.com/facebookresearch/esm
We adapted their code and used it to embed sequences in embed_seqs.py.

The readme for DEDAL can be found here, which contains instructions for using the model:
https://github.com/google-research/google-research/tree/master/dedal
We adapted their code and used it to run their model in get_aligns.py.

FATCAT github can be found here:
https://github.com/GodzikLab/FATCAT-dist

vcMSA github can be found here:
https://github.com/clairemcwhite/vcmsa

**************************************************************************************************************
# Results and Figures
**************************************************************************************************************

The 'data/alignments' folder contains all of the pairwise alignments that we generated and used for analysis.
The 'figures' folder contains several scripts used to generate the figures and tables in the paper.

**************************************************************************************************************
# License
**************************************************************************************************************

Licensed under the Academic Free License version 3.0. 
