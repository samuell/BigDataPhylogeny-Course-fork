# INTRODUCTION

In this practical session we will learn and work on the most important steps of a phylogenetic reconstruction pipeline. The pipeline we will work with consists of three main steps: 
  1) multiple sequence alignment (MSA), 
  2) alignment trimming, and 
  3) phylogenetic inference. 

The outcome of this methodological pipeline is a phylogenetic tree which we will be able to inspect to assess its robustness and the phylogenetic relationships between the sequences represented.

## 1. Alignment
The alignment is a useful tool to verify that we are comparing homologous positions in our sequences. Once you have some certainty that the sequence set belongs to genes that have a common ancestry 
(i.e., they are homologous, all descend from a given ancestral gene), the alignment is also important because it allows the identification of homologous sites between sequences. A molecular phylogeny 
can be reconstructed by looking at the similarities and differences between the aligned sequences on each site. There are several alignment methods, from basic to progressive algorithms. In this section, 
we will focus on what should be inspected from a multiple sequence alignment to ensure a reliable phylogenetic inference process.
We should be aware of the fact that alignments are rarely perfect. Below we can see two examples of this:


Example A. An example of a shifted fragment. It is caused by the falsely homologous positions in the green columns that are highlighted in red (threonine, T).


### Objectives
1. Manipulate and align homologous sequences
2. Obtain definitive matrices for each gene (gene family) and concatenate them for subsequent phylogenetic analyses

### Programs
• MEGA: Sequence manipulation and alignment
It is freely distributed and you can obtain the latest version (12.1) at: https://www.megasoftware.net/

### Data
• You have the data from the previous OrthoFinder run or you can download them from the materials folder

The data is in FASTA format.

You can download them from the folder "sequences".

