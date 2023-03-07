**************************************************************************************************************
# Protein Embedding Based Alignments (PEbA)
**************************************************************************************************************

This project uses Rostlab's ProtT5-XL-UniRef50 encoder (https://huggingface.co/Rostlab/prot_t5_xl_uniref50) to
embed protein sequences, and then calculates the cosine similarity between each residue's vector to align them 
using traditional sequence alignment algorithms. The cosine similarity scores are used in place of substitution
matrices, such as BLOSUM, with the goal of producing accurate alignments for sequences that share low character
identity (<20-35%), referred to as the "twilight zone" of protein sequence alignment.

PEbA is implemented in local_PEbA.py, where it uses the Smith-Waterman 'local' algorithm to align two sequences.
The script for substitution matrix based alignments are in local_MATRIX.py, which currently supports all BLOSUM 
matrices and PFASUM60 for testing purposes.

The PEbA script is designed to be run from the command line. The following arguments are allowed, with -embed1
-embed2 being optional arguments if the corresponding sequences from -file1 and -file2 have already been
been embedded and saved as a numpy array:

    -file1 <file1.fa>   : Path to first sequence fasta file
    -file2 <file2.fa>   : Path to second sequence fasta file
    -embed1 <file1.txt> : Path to first sequence embedding file
    -embed2 <file2.txt> : Path to second sequence embedding file
    -gopen <int/float>  : Gap opening penalty
    -gext <int/float>   : Gap extension penalty
    -encoder <str>      : Encoder used to generate embeddings

Embeddings from any model can be used assuming the model outputs a 1D array the same length as the sequence.
The -encoder argument is used to specify the model used and written to the output .msf file for reference.

**************************************************************************************************************
# Comparing PEbA to Reference Alignments
**************************************************************************************************************

To determine if PEbA alignments produce more accurate alignments than those with, subsitution matrices, both 
types of alignments are compared to reference alignments from BAliBASE 4 (https://www.lbgi.fr/balibase/).
BAliBASE alignments are structure based benchmarks used to compare new methods of alignment, such as this one. 
The references of most interest include RV11, sequences with <20% identity, and RV12, sequences with 20-40% 
identity, but other references with higher % identity are also used for comparison.

t_coffee's "aln_compare" function is used to compare the performance between BLOSUM and PEbA alignments to the
reference alignments. An example of aln_compare output:

*****************************************************
seq1       seq2          Sim        ALL           Tot  
PEbA_2        2           7.9    78.9 [100.0]   [  152]

As explained in their documentation (https://tcoffee.org/Projects/tcoffee/documentation/index.html) seq1 is the 
alignment being compared; seq2 is the number of sequences in the alignment; Sim is the average identity of sequences 
found in the first alignment (seq1); ALL is the fraction of columns in the first alignment that are found identically 
aligned in the second (100.0 meaning the percentage of columns being used in the comparison); and Tot is the number of 
paired residues contained in the first alignment. A higher score in the 4th column, under ALL, would indicate a more 
similar alignment to the one being compared. In this example, the 78.9% of columns in the PEbA alignment were found 
identically aligned in the reference BAliBASE alignment, both pairwise alignments of two short protein sequences.