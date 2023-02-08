**************************************************************************************************************
# Protein Embedding Based Alignments
**************************************************************************************************************

This project uses Rostlab's ProtT5-XL-UniRef50 encoder (https://huggingface.co/Rostlab/prot_t5_xl_uniref50) to
embed protein sequences and then calculates the cosine similarity between residues to align them using 
traditional Needleman-Wunsch and Smith-Waterman algorithms. The cosine similarity scores are used in place of 
substitution matrices, such as BLOSUM, with the goal of producing accurate alignments for sequences that share 
low character identity (<20%), referred to as the "twilight zone" of protein sequence alignments.

To determine if protein embedding based alignments (PEbA) produce more accurate alignments than those with
subsitution matrices, both types of alignments are compared to reference alignments from BAliBASE 4
(https://www.lbgi.fr/balibase/). BAliBASE alignments are structure based benchmarks used to compare new methods
of alignment, such as this one. The current references of interest include RV11, sequences with <20% identity,
and RV12, sequences with 20-40% identity.

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