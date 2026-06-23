# Lab 2 Inferring Orthology

## Objective and data

We are ready to start our first phylogenomic study.
We want to study the evolution of molluscs, based on this [paper](https://www.science.org/doi/10.1126/science.ads0215).

Due to time and computational constrains, we will work on a subset of this dataset. We will focus on 20 proteomes.

First we connect to the server as explained in the instruction. 

```
#Copy the proteomes folder in your own personal folder 
cp -r /home/user_number/Share/Proteomes .
```

These are the proteomes for the molluscs species we want to infer the phylogeny for.
To check if you have all of them, go into the directory using `cd` and count how many files you have. You can use `ls * | wc -l`.

We are ready to infer the orthology, let's go to the following section!

## Orthofinder

To identify homologs among all the proteins, we will use [OrthoFinder](https://github.com/davidemms/OrthoFinder) with default parameters, just providing the folder containing the proteome files.

First we need to activate the conda env with Orthofinder

```
conda activate Orthofinder3
```
This is a virtual environment, with different softwares installed that will help us running Orthofinder and all its dependencies.
If you want to know more about conda you can have a look at the [cheat sheet here](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf).

Now we are ready to infer our orthologous genes, it will take a while to run. Lunch it and please keep going with the tutorial, while it runs to completion. We will go back to it at the end of this practical session.
```
orthofinder -f Proteomes/
```
It will take 40 minutes to run, let it run, you will go back to the output later.


## Comparing trees

 The main output of a phylogenomic study is a tree, so reading trees correctly is fundamental. 
 
 Here you have different trees, have a look at them and try to answer some of the questions 

![alt text](Nosenko.png)

1) List at least 3 differences you can see
2) Can you identify the paraphyletic group?
3) Is Metazoa monophyletic in these trees?

![alt text](treeexamples.png)

4) Can you select which branches could be paraphyletic?



![alt text](fishtree.png)

5) Can you identify a polyphyly group?
6) Among the taxa represented, can you identify a group that is famously paraphyletic?

![alt text](Redmond_competingtrees.webp)

7) Focusing on section b and c, which topology is correct?
8) Pay attention to Placozoa. Can you describe how this taxon move in all different trees?

![alt text](Lucatree.webp)

10) Pay close attention to the nodes and the branches of this tree. 
	- What are the authors of this paper trying to do here? 
	- Do you see any unxepected pattern depicted in this topology?

Ok we have now familiarised with trees and tree thinking. It's time to make our own. 
Let's move on the following section!

## Selecting orthologs
Orthofinder is a very powerful tool, but can be easily mislead users to approach it as a black box.

If it has finished running, check the output by exploring all the different folders that has created. 
If not, we have a folder prepared that you can access and look at the output. Kill the previous run by typing `CTRL+C`.


Our goal is to identify the orthologous genes that we will use to create the super alignment to infer the mollusc phylogeny. 
In our case, the repository where we can find them is `/home/user_number/Share/Single_Copy_Orthologue_Sequences` 

If you have good quality genomes, normally you will end up with a set of genes that orthofinder has identified as single copy orthologs (i.e. true orthologs) and you could ideally start you analysis from these set of genes. Supposedely, these genes do not contain paralogs. 
In case you are not so lucky (especially when using a lot of transcriptomes), most likely you will have to fish out you genes from the `Orthogroups` of `Orthologues` folders.

Here we were lucky enough that we have retrieve single copies, we will now align, trim, and infer the trees.


## Multiple sequence alignment

Due to time constraints we have selected a subset of the single copy orthologs that were identified. 

```sh
#First we deactivate the Orthofinder environment
conda deactivate

#We now activate the environment that we are going to use for the phylogenomic inference.
conda activate BigDataPhylo
#if you want to know which softwares are installed in the environment you can type
conda list #it will provide a list of all packages and dependencies you will use in the next practicals
#Copy the entire directory in your home directory
cp -r /home/user_number/Share/Single_Copy_Orthologue_Sequences .
```
You should now have all the fastas in the repository. Enter into the directory using `cd` and check with `ls *fa | wc -l `
It should tell you how many you have.

The next step is to infer multiple sequence alignments from orthogroups. Multiple sequence alignments allow us to *propose* which amino acids/ nucleotides are homologous. A simple yet accurate tool is [MAFFT](https://mafft.cbrc.jp/alignment/server/).


We will align gene files separately using MAFFT in a for loop:

```sh
for f in *fa; do mafft $f > $f.mafft; done
```

To check if all your fastas have aligned count how many *.mafft files have been created. It should be the same number as the fasta files. 

## Alignment trimming


Some gene regions (e.g., fast-evolving) are difficult to align and thus positional homology can be uncertain. It is unclear (i.e., problem-specific) whether trimming suspicious regions [improves](https://academic.oup.com/sysbio/article/56/4/564/1682121) or [worsens](https://academic.oup.com/sysbio/article/64/5/778/1685763) tree inference. However, gently trimming very incomplete positions (e.g. with >90% gaps) will speed up computation in the next steps without a significant loss of phylogenetic information.

To trim alignment positions we can use [TrimAl](https://trimal.readthedocs.io/en/latest/) but several other software are also available.

From the same directory run:

```sh
for a in *.mafft;

	do trimal -in $a -out $a.trim  -fasta -gappyout;

done

```
While diving into phylogenomic pipelines, it is always advisable to check a few intermediate results to ensure we are doing what we should be doing. Multiple sequence alignments can be visualized in [SeaView](http://doua.prabi.fr/software/seaview) or [AliView](https://github.com/AliView/AliView). Also, one could have a quick look at alignments using command line tools (`less -S`).

Well done on your first day! You have a first dataset prepared, tomorrow we will lear how to infer trees, with Maximum likelihood, Bayesian inference, and how to check for paralogs.
