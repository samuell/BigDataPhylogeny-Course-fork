# Lab 2 Inferring Orthology

OrthoFinder.
 Orthofinder output interpretation.  MAFFT alignment, trimming (Trimal). Alignment visualisation (JalView). 



## Objective and data

We are ready to start our first phylogenomic study.
We want to study the evolution of molluscs, based on this [paper](https://www.science.org/doi/10.1126/science.ads0215)

Due to time and computational constrains, we will work on a subset of this dataset. We will focus on 20 proteomes.

First we connect to the server
!!!ADD CONNECTION AMAZON CLUSTER!!!
(upload fasta files so they can copy in their repo)
```
#Copy the proteomes folder in your own personal folder 
cp BigDataPhylogeny/proteomes ~/.

```

## Inferring ortholog groups

To identify homologs among all the proteins, we used [OrthoFinder](https://github.com/davidemms/OrthoFinder) with default parameters, just providing the folder containing the proteome files.

First we need to activate the conda env with Orthofinder

```
conda activate Orthofinder
```
This is a virtual environment, with different softwares installed that will help us running Orthofinder and all its dependencies.
If you want to know more about conda you can have a look at the [cheat sheet here](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf).

Now we are ready to infer our orthologous genes, it will take a while to run. Lunch it and please keep going with the tutorial, while it runs to completion. We will go back to it at the end of this practical session.
```
orthofinder -a 8 -f /user/work/yp19290/BigData_physalia/Proteomes 
```
It will take 40 minutes to run, let it run, you will go back to the output later.


## Comparing trees

 The main output of a phylogenomic study is a tree, so reading trees correctly is fundamental. 
 
 Here you have two competing topologies with the same species. 

<img Nosenko>

- List at least 3 differences you can see
- Can you identify the paraphyletic group?

- Is Metazoa monophyletic in these trees?

<img Simple tree>

- Can you select which branches could be paraphyletic?

- Can you identify a polyphyly group?

- Which topology is correct?




 
 Reading and comparing trees (same species, different topologies). Identifying monophyly, paraphyly, polyphyly, sister groups, incongruences, etc.

Check if Orthofinder is still running, if it has stopped move to the next section of the tutorial.


## Pre-alignment and quality filtering
Orthofinder is a very powerful tool, but can be easily mislead users to approach it as a black box.

Check the output by exploring all the different folders that has created. 
In our case, the one we care the most is `/BigData_physalia/Proteomes/OrthoFinder/Results_Jun03/Single_Copy_Orthologue_Sequences` 

If you have good quality genomes, normally you will end up with a set of genes that orthofinder has identified as single copy orthologs (i.e. true orthologs) and you could ideally start you analysis from these set of genes. Supposedely, these genes do not contain paralogs. 
In case you are not so lucky (especially when using a lot of transcriptomes), most likely you will have to fish out you genes from the `Orthogroups` of `Orthologues` folders.

Here we were lucky enough that we have retrieve single copies, we will need to align them and infer the trees.

## Multiple sequence alignment

Copy the fasta files in a new directory in the home directory

```sh
mkdir Mollusca_orthologs
cd Mollusca_orthologs

cp /OrthoFinder/Results_Jun03/Single_Copy_Orthologue_Sequences/*fas .
```

You should now have all the fastas in the repository. Check with `ls *fas | wc -l `
it should tell you how many you have.

The next step is to infer multiple sequence alignments from orthogroups. Multiple sequence alignments allow us to *propose* which amino acids/ nucleotides are homologous. A simple yet accurate tool is [MAFFT](https://mafft.cbrc.jp/alignment/server/).

We will align gene files separately using a for loop:

```
for f in *fas; do mafft $f > $f.mafft; done
```

To check if all your fastas have aligned count how many *.mafft files have been created.



## Alignment trimming


Some gene regions (e.g., fast-evolving) are difficult to align and thus positional homology can be uncertain. It is unclear (i.e., problem-specific) whether trimming suspicious regions [improves](https://academic.oup.com/sysbio/article/56/4/564/1682121) or [worsens](https://academic.oup.com/sysbio/article/64/5/778/1685763) tree inference. However, gently trimming very incomplete positions (e.g. with >90% gaps) will speed up computation in the next steps without a significant loss of phylogenetic information.

To trim alignment positions we can use [TrimAl](https://trimal.readthedocs.io/en/latest/) but several other software are also available.


```
for a in *.mafft;

	do trimal -in $a -out $a.trim  -fasta -gappyout;

done

```

While diving into phylogenomic pipelines, it is always advisable to check a few intermediate results to ensure we are doing what we should be doing. Multiple sequence alignments can be visualized in [SeaView](http://doua.prabi.fr/software/seaview) or [AliView](https://github.com/AliView/AliView). Also, one could have a quick look at alignments using command line tools (`less -S`).



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
  