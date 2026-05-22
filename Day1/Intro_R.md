Introduction to R (ape, ggtree)

by Marta Álvarez Presas

During the course we will use R to explore and visualise phylogenetic trees. 
R is a statistical programming environment widely used in bioinformatics, especially for analyses and visualisation. 
In this tutorial, we will focus on how to read, manipulate and visualise phylogenetic trees using the packages ape and ggtree.

We start by opening R (or RStudio). You will see a console where you can type commands.
>

This is the prompt, where you can type commands and run them by pressing Enter.

# Running commands
In R, you type commands line by line. You can also assign results to objects using the arrow <-.
For example:
a <- 3a
This assigns the value 3 to the variable a and prints it.

# Context and setup
Before working with trees, we need to install and load the required packages.

# Installation (only needed once)
install.packages("ape")
install.packages("ggplot2")
if (!requireNamespace("BiocManager", quietly = TRUE))  
	install.packages("BiocManager")
BiocManager::install("ggtree", update = TRUE, ask = FALSE)
	
# Load libraries (needed every session)
library(ape)
library(ggtree)
Packages in R are collections of functions.

ape is used to manipulate phylogenetic trees
ggtree is used to visualise them


# Working with phylogenetic trees
Creating and visualising a tree
We can create a simple random tree using the function rtree():

tree <- rtree(10)
plot(tree)

This will generate a tree with 10 tips (taxa) and display it.
What do you see?


How many tips does the tree have?
Can you identify the branches and nodes?


Reading a tree from file
In bioinformatics, trees are often stored in Newick format.
We can create a simple example:
tree <- read.tree(text="((A:1,B:1):2,C:3);")
plot(tree)

In this tree:

A, B, and C are taxa
numbers represent branch lengths


# Exploring the tree structure
Let's inspect the object:
tree

You should see information about:

number of tips
number of nodes
branch lengths
labels


# Basic operations on trees
Distances between taxa
We can compute pairwise distances:
d <- cophenetic(tree)
d

Example output:
  A B C
A 0 2 6
B 2 0 6
C 6 6 0

These values represent evolutionary distances between taxa.

Removing taxa
Sometimes we want to prune a tree:
pruned <- drop.tip(tree, "A")
plot(pruned)

What changed compared to the original tree?

# Rooting a tree
We can root the tree using a specific taxon:
rooted <- root(tree, "C")
plot(rooted)

Rooting changes the interpretation of evolutionary relationships.

# Exercise
Remove two species from a tree and compare the resulting topology.

Visualising trees with ggtree
The plot() function is simple but limited.
The ggtree package allows more flexible and publication-quality figures.

# Basic visualisation
ggtree(tree)

This creates a tree using the grammar of graphics.

Adding tip labels
ggtree(tree) +  
geom_tiplab(size = 4)

This shows the names of the taxa.

# Adding nodes and aesthetics
ggtree(tree) +  
geom_tiplab() +  
geom_nodepoint(color = "red") +  
theme_tree()

What do the red points represent?

# Changing layouts
We can visualise the same tree with different layouts:
ggtree(tree, layout = "circular") +  
geom_tiplab()

Other layouts are available (rectangular, slanted, etc.).

# Exercise
Create a circular tree with:

tip labels
highlighted nodes


# Useful tips

In R, you can get help using:

?rtree
?ggtree


Errors are normal — read them carefully and try to understand what is missing


Always check your objects:

class(tree)

# Final notes
This tutorial provides only the basic tools to work with phylogenetic trees in R.
During the course, you will combine these skills with real genomic datasets and more complex analyses.
Good luck and have fun! 
