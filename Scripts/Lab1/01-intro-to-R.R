### ============================================================
### Practical 1 - Introduction to R for Phylogenetics
### ============================================================

### ----1. R basics ----

## 1.1 R as a calculator and assignment
2 + 2
10 / 3
2^5          # exponent

x <- 5
y <- 10
x + y

## 1.2 Data types
a <- 4.5             # numeric
b <- "Homo sapiens"  # character (string)
d <- TRUE            # logical
e <- c(1, 2, 3, 4)   # vector of numbers
f <- c("ape", "ggtree", "phytools")  # vector of strings

class(a)
class(f)
length(e)

## 1.3 Vectors and indexing
species <- c("Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Gallus_gallus")
species[1]
species[2:3]
length(species)

## 1.4 Data frames
traits <- data.frame(
  species = c("Homo_sapiens", "Pan_troglodytes", "Mus_musculus"),
  body_mass_kg = c(70, 50, 0.02),
  genome_size_Mb = c(3100, 3300, 2700)
)

traits
str(traits)
traits$body_mass_kg
nrow(traits)

## 1.5 Functions and help
mean(traits$body_mass_kg)
sqrt(16)

# ?mean        # uncomment to open help page
# help(mean)

## 1.6 Working directory
getwd()
# setwd("path/to/folder")   # uncomment and edit if needed


### ---- 2. Installing packages ----

## 2.1 CRAN packages (run once per computer)
install.packages("ape")
install.packages("phytools")
install.packages("ggplot2")
install.packages("BAMMtools")

## 2.2 Bioconductor package
install.packages("BiocManager")
BiocManager::install("ggtree")

## 2.3 Load packages (run every session)
library(ape)
library(ggtree)
library(phytools)
library(ggplot2)


### ---- 4. Your first phylogenetic tree ----

set.seed(123)
tree <- rtree(10)        # simulate a random tree with 10 tips

class(tree)
tree
str(tree)

# Base ape plotting
plot(tree)

# ggtree plotting
ggtree(tree) + geom_tiplab()

### ============================================================
### End of Practical 1
### ============================================================
