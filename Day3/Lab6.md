# Reconciling Gene and Species Trees (DTL), Ancestral Sequence Reconstruction (ASR), and Comparing Evolutionary Scenarios

## INTRODUCTION
In modern comparative genomics, the challenge is no longer just to obtain sequences, but to understand the complex history behind each gene. Although we often assume that genes follow the path of species (vertical inheritance), phenomena such as duplication, horizontal gene transfer (HGT), and gene loss (DLT) create constant phylogenetic conflicts.

This lab will allow you to move beyond individual sequence analysis to enter the world of large-scale phylogenomic reconciliation. We will use data from previous labs to explore how different evolutionary events shape the gene repertoire of extant species.

### 0. Before you start: what is this dataset, really?

The files you will use today look like a small clade of molluscs (clams, snails, octopuses, chitons...) on the same 20-taxon tree used in Chen et al.'s genome-based mollusc phylogeny in Science. That is intentional, but  the actual protein sequences underneath are real, public sequences from twenty different strains of the bacterium *Streptococcus*, simply relabelled with mollusc species names so that the tree topology and branch pattern you'll work with come from real, published, peer-reviewed phylogenetic relationships rather than an arbitrary random tree.

Why bother with this disguise instead of just calling the species "Streptococcus sp. 1, 2, 3..."? Two reasons:
  - First, it keeps the focus on the method (how do you reconcile a gene tree with a species tree? how do you reconstruct an ancestral sequence?) rather than on bacterial taxonomy, which is not the point of today's lab.
  - Second, every species in the tree has a distinct, memorable identity, which makes it much easier to talk about "what happened on the branch leading to the octopuses" than "what happened on branch 14."

None of the biological conclusions you reach today (which mollusc gained a gene, which lineage received a horizontal transfer) are real biology — they describe what happened to relabelled bacterial genes, or in two cases, to gene families we deliberately engineered to contain a duplication, a loss, or a transfer so that you have something to find. Section 1 explains exactly which parts of the data are real and which are constructed, so you always know what you are looking at.

### 1. Objectives
The main objective is that you learn to use probabilistic modeling tools to reconstruct the evolutionary history of complex gene families. Specifically, in this session you will achieve the following objectives:

  #### Reconciliation by Likelihood: 
  You will learn to use GeneRax to resolve the conflict between gene trees and species trees using maximum likelihood models.

  #### ASR (Ancestral State Reconstruction): 
  You will infer the ancestral states of gene characters to understand which genetic traits were present in the ancestors of the main mollusk clades.

  #### Comparison of evolutionary scenarios: 
  You will evaluate the robustness of different species topologies by analysing the statistical likelihood of DTL scenarios.

  #### Biological interpretation: 
  You will learn to differentiate between phylogenetic noise (error in tree inference) and real biological signals (transfers and duplications).

### 2. Data overview

  #### 2.1 The species tree

  "Data/Mollusk_SPP_tree.tre" contains a fixed, rooted, bifurcating tree of 20 species, with *Lingula anatina* (a brachiopod, not a mollusc) as the outgroup, matching the outgroup choice used in the source study. Internal   nodes are pre-labelled N1 to N19 so that you can refer to specific ancestors by name during the exercises (for example, "the ancestor of the cephalopods is N14").

  You can look at the tree by opening it in any tree viewer that reads Newick (FigTree, iTOL, ETE3) as we have done in previous lab sessions.

  #### 2.2 The gene set (real data, single-copy orthologues)

  Data/core_gene_set/ contains six real, single-copy orthologous protein-coding genes, each present in all 20 species exactly once. They were chosen because they are short, well-conserved, easy-to-align genes — exactly   the kind of "boring, reliable" gene family you would use to first sanity-check a pipeline in a real research project, before tackling messier, larger gene families.


## DTL analysis
### Software: 
GeneRax
IQ-TREE (to generate the gene trees), 
Python 3.x with ete3.



- Gene families (G): For each family you need:

    ".fa": Sequence alignment. As starting data we will use the single-copy orthologous families obtained with OrthoFinder (Single_Copy_Orthologue_Sequences).

    ".nwk": "Raw" gene tree 

- Mapping file: A "mapping.txt" file that connects sequence names with species names (format: species:gene1,gene2).

### Practical
#### Gene tree processing
- Decompress gene family sequences using tar:
  ```
  tar -zxvf Single_Copy_Sequences.tar.gz
  ```

- First of all we need to clean headers for further analyes. With that purpose, download the script *"clean_headers.py"* from the "Scripts" folder and run it like this:

  ```
  for f in *.fa; do python clean_headers.py $f; done
  ```
  
- Align sequences:
  ```
  mkdir -p aln
  for f in *.fa; do mafft --auto $f > aln/${f%.fa}_aln.fa; done
  ```

- Trim alignments:
  ```
  mkdir -p trim
  for f in *.fa; do trimal -in $f -out trim/${f%.fa}_trim.fasta -gappyout; done
  ```

- Infer gene trees:

  ```
  for f in *.fasta; do iqtree -s $f -m LG+G -bb 1000 -nt AUTO; done
  ```

```
generax --families families.txt --species-tree S_tree.nwk --prefix output --rec-model UndatedDTL --strategy RESOLVE
```
