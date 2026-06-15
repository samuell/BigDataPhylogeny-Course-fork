# DTL, ASR and comparison of evolutionary scenarios

## INTRODUCTION
In modern comparative genomics, the challenge is no longer just to obtain sequences, but to understand the complex history behind each gene. Although we often assume that genes follow the path of species (vertical inheritance), phenomena such as duplication, horizontal gene transfer (HGT), and gene loss (DLT) create constant phylogenetic conflicts.

This lab will allow you to move beyond individual sequence analysis to enter the world of large-scale phylogenomic reconciliation. We will use data from previous labs to explore how different evolutionary events shape the gene repertoire of extant species.

### Objectives
The main objective is that you learn to use probabilistic modeling tools to reconstruct the evolutionary history of complex gene families. Specifically, in this session you will achieve the following objectives:

  #### Reconciliation by Likelihood: 
  You will learn to use GeneRax to resolve the conflict between gene trees and species trees using maximum likelihood models.

  #### ASR (Ancestral State Reconstruction): 
  You will infer the ancestral states of gene characters to understand which genetic traits were present in the ancestors of the main mollusk clades.

  #### Comparison of evolutionary scenarios: 
  You will evaluate the robustness of different species topologies by analysing the statistical likelihood of DTL scenarios.

  #### Biological interpretation: 
  You will learn to differentiate between phylogenetic noise (error in tree inference) and real biological signals (transfers and duplications).

## DTL analysis
### Software: 
GeneRax
IQ-TREE (to generate the initial tree), 
Python 3.x with ete3.

### Required data:

- Species tree (S): A rooted Newick file (species_tree.nwk).

- Gene families (G): For each family you need:

    ".fasta": Sequence alignment.

    ".nwk": "Raw" gene tree obtained with:
  ```
  iqtree -s family.fasta -m LG+G
  ```

- Mapping: A "mapping.txt" file that connects sequence names with species names (format: species:gene1,gene2).


```
generax --families families.txt --species-tree S_tree.nwk --prefix output --rec-model UndatedDTL --strategy RESOLVE
```
