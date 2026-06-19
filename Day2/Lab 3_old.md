# INTRODUCTION

In this practical session we will learn and work on the most important steps of a phylogenetic reconstruction pipeline. The pipeline we will work with consists of three main steps: 
  1) multiple sequence alignment (MSA), 
  2) alignment trimming, and 
  3) phylogenetic inference. 

The outcome of this methodological pipeline is a phylogenetic tree which we will be able to inspect to assess its robustness and the phylogenetic relationships between the sequences represented.

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

## 1. Alignment
The alignment is a useful tool to verify that we are comparing homologous positions in our sequences. Once you have some certainty that the sequence set belongs to genes that have a common ancestry 
(i.e., they are homologous, all descend from a given ancestral gene), the alignment is also important because it allows the identification of homologous sites between sequences. A molecular phylogeny 
can be reconstructed by looking at the similarities and differences between the aligned sequences on each site. There are several alignment methods, from basic to progressive algorithms. In this section, 
we will focus on what should be inspected from a multiple sequence alignment to ensure a reliable phylogenetic inference process.
We should be aware of the fact that alignments are rarely perfect. Below we can see two examples of this:


#### Example A. An example of a shifted fragment. It is caused by the falsely homologous positions in the green columns that are highlighted in red (threonine, T).
<img width="932" height="358" alt="image" src="https://github.com/user-attachments/assets/87a4511d-5a0a-40ff-bf16-0e0247d831b5" />

### Example B. Indels (insertions and deletions) can complicate the alignment process. In this example, a single insertion has caused a misalignment of multiple sites in one of the sequences.
<img width="445" height="411" alt="image" src="https://github.com/user-attachments/assets/898deb0d-d3e6-4bde-b2c0-b8da55cca5da" />


## 2. Trimming
Phylogenetic reconstruction is a process that can require several computational resources. Moreover, it is very likely that not every position will be informative. Thus, it is important to remove positions that do not add enough information to the alignment but add noise instead. There are two main methods: 1) checking the proportion of gaps for a given alignment position (e.g., trimAl software); and 2) checking the entropy associated with each position (e.g., BMGE software). In some cases, gappy regions in the alignment can include interesting information such as a protein domain that could have evolved de novo in a small subset of species.

<img width="693" height="236" alt="image" src="https://github.com/user-attachments/assets/28b54703-62e4-4e3c-b6a0-05cf5b7f0afc" />

## 3. Phylogenetic tree inference
Once we have a good alignment (and have trimmed it, if needed), we can infer the phylogenetic relationships between the aligned sequences. This reconstruction consists of grouping sequences by similarity and inferring the branch lengths, which is an estimation of the number of substitutions per site under a given evolutionary model.

It is important to differentiate between a species tree and a gene tree. A species tree aims to represent the divergence process between a set of species, while a gene tree aims to represent the evolutionary history of a given gene family. Gene trees often differ from the species tree due to the occurrence of gene duplications, gene losses or horizontal gene transfer events during the evolution of the gene family. In this particular session, we will focus on reconstructing a gene tree. In future sessions we will discuss species tree reconstruction.

Once the phylogenetic tree has been reconstructed, we can inspect the topology to perform inferences about the evolutionary history of the gene family. For instance, if the sequences that descend from a given node in the tree belong to the same species, the most plausible explanation is that this node of the tree is representing a recent duplication that occurred in the evolutionary path towards this species.

## 4. Tree support and root
Bootstrapping is the most used statistical tool to measure the consistency of the recovered tree. Briefly, this consists of randomly subsampling the alignment sites multiple times (usually >= 100 times) and re-infer a phylogeny from each subsampled alignment. The bootstrap support of every clade in the phylogeny is measured by checking the fraction of bootstrapped trees that recovered each clade.

<img width="940" height="413" alt="image" src="https://github.com/user-attachments/assets/7d6c3513-7258-44be-a65e-ccd5bb63d508" />
(Yang and Rannala, 2012, www.doi.org/10.1038/nrg3186)

The interpretation of a phylogenetic tree often requires us to root it. The root can be set by adding a group of sequences (outgroup sequences) belonging to species that are close enough but do not belong to our group of interest (ingroup); this method is known as outgroup rooting. A less optimal choice is to root the tree at the branch that is in the middle of the most distant path between two tips (midpoint rooting).





## Concatenate alignment


To infer our phylogenomic tree we need to concatenate the trimmed single-gene alignments we generated. There are many tools that you can use for this step (e.g [concat_fasta.pl](https://github.com/santiagosnchez/concat_fasta) or [catsequences](https://github.com/ChrisCreevey/catsequences)). Here, we will use [FASconCAT](https://github.com/PatrickKueck/FASconCAT-G), which will read in all `\*.fas` `\*.phy` or `\*.nex` files in the working directory and concatenate them (in random order).

```
for f in *bmge; do mv $f $f.fas; done
mkdir concatenation/
mv *bmge.fas concatenation/

cd concatenation/
perl ~/Desktop/software/FASconCAT-G_v1.04.pl -l -s
```

Is your concatenated file what you expected? It should contain 23 taxa and 21 genes. You might check the concatenation (`FcC_supermatrix.fas`) and the file containing the coordinates for gene boundaries (`FcC_supermatrix_partition.txt`). Looking good? Then your concatenated dataset is ready to rock!!



## Concatenation: Maximum likelihood


One of the most common approaches in phylogenomics is gene concatenation: the signal from multiple genes is "pooled" together with the aim of increasing resolution power. This method is best when among-gene discordance is low.

We will use [IQTREE](http://www.iqtree.org/), an efficient and accurate software for maximum likelihood analysis. Another great alternative is [RAxML](https://github.com/stamatak/standard-RAxML). The most simple analysis is to treat the concatenated dataset as a single homogeneous entity. We need to provide the number of threads to use (`-nt 1`) input alignment (`-s`), tell IQTREE to select the best-fit evolutionary model with BIC (`-m TEST -merit BIC -msub nuclear`) and ask for branch support measures such as non-parametric bootstrapping and approximate likelihood ratio test (`-bb 1000 -alrt 1000 -bnni`):

```
iqtree2 -s FcC_supermatrix.fas -m TEST -msub nuclear -bb 1000 -alrt 1000 -nt AUTO -bnni -pre unpartitioned
```

A more sophisticated approach would be to perform a partitioned maximum likelihood analysis, where different genes (or other data partitions) are allowed to have different evolutionary models. This should provide a better fit to the data but will increase the number of parameters too. To launch this analysis we need to provide a file containing the coordinates of the partitions (`-p`) and we can ask IQTREE to select the best-fit models for each partition, in this case, according to AICc (more suitable for shorter alignments).

```
iqtree2 -s FcC_supermatrix.fas -p FcC_supermatrix_partition.txt -m TEST -msub nuclear -merit AICc -bb 1000 -alrt 1000 -nt AUTO -bnni -pre partitioned
```

Congratulations!! If everything went well, you should get your maximum likelihood estimation of the vertebrate phylogeny (`.treefile`)! Looking into the file you will see a tree in parenthetical (newick) format. See below how to create a graphical representation of your tree.


## Coalescence analysis


An alternative to concatenation is to use a multispecies coalescent approach. Unlike maximum likelihood, coalescent methods account for incomplete lineage sorting (ILS; an expected outcome of evolving populations). These methods are particularly useful when we expect high levels of ILS, e.g. when speciation events are rapid and leave little time for allele coalescence.

We will use [ASTRAL](https://github.com/smirarab/ASTRAL), a widely used tool that scales up well to phylogenomic datasets. It takes a set of gene trees as input and will generate the coalescent "species tree". ASTRAL assumes that gene trees are estimated without error.

Thus, before running ASTRAL, we will need to estimate individual gene trees. This can be easily done by calling IQTREE in a for loop:

```
for f in *bmge.fas; do iqtree2  -s $f -m TEST -msub nuclear -merit AICc -nt AUTO; done
```

After all gene trees are inferred, we should put them all into a single file:

```
cat *bmge.fas.treefile > my_gene_trees.tre
```

Now running ASTRAL is trivial, providing the input file with the gene trees and the desired output file name:

```
java -jar ~/Desktop/software/Astral/astral.5.7.8.jar -i my_gene_trees.tre -o species_tree_ASTRAL.tre 2> out.log
```

Congratulations!! You just got your coalescent species tree!! Is it different from the concatenated maximum likelihood trees? 


## Tree visualization


Trees are just text files representing relationships with parentheses; did you see that already? But it is more practical to plot them as a graph, for which we can use tools such as [iTOL](https://itol.embl.de), [FigTree](https://github.com/rambaut/figtree/releases), [TreeViewer](https://treeviewer.org/), [iroki](https://www.iroki.net/), or R (e.g. [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html), [phytools](http://www.phytools.org/); see provided R code).

Upload your trees to iTOL. Trees need to be rooted with an outgroup. Click in the branch of *Callorhinchus milii* and the select "Tree Structure/Reroot the tree here". Branch support values can be shown under the "Advanced" menu. The tree can be modified in many other ways, and finally, a graphical tree can be exported. Similar options are available in FigTree.

[Well done!](https://media.giphy.com/media/wux5AMYo8zHgc/giphy.gif)


## Software links

* Orthofinder (https://github.com/davidemms/OrthoFinder)
* PREQUAL (https://github.com/simonwhelan/prequal)
* MAFFT (https://mafft.cbrc.jp/alignment/software/source.html)
* MUSCLE v5 (https://github.com/rcedgar/muscle)
* TrimAL (https://vicfero.github.io/trimal/)
* FASTCONCAT (https://github.com/PatrickKueck/FASconCAT-G)
* IQTREE (http://www.iqtree.org/)
* Phylobayes (https://github.com/bayesiancook/phylobayes/tree/master)
* ASTRAL (https://github.com/smirarab/ASTRAL)
* FIGTree V1.4.4 (https://github.com/rambaut/figtree/releases) 
* TreeViewer (https://treeviewer.org/)
* iTOL (https://itol.embl.de/)
* SeaView (https://doua.prabi.fr/software/seaview_data/seaview5-64.tgz)
  