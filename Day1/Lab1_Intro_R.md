# Introduction to R for Phylogenetics & Comparative Genomics

by Marta Álvarez Presas

During the course we will use R to explore and visualise phylogenetic trees. 
R is a statistical programming environment widely used in bioinformatics, especially for analyses and visualisation. 
In this tutorial, we will focus on how to read, manipulate and visualise phylogenetic trees using the packages ape and ggtree.

## Objectives
* Get comfortable with `RStudio` and `R basics`, install the packages we'll use during the course, and load your very first phylogenetic tree.

> If you already know R well, skim this fast and jump straight to [Section 4](#part-4-your-first-phylogenetic-tree-5-min). Everyone should run Section 2 (package installation) before the next practical.

---

## 0. Before we start

We start by opening **RStudio**. You should see four panels:

| Panel | What it's for |
|---|---|
| **Source** (top-left) | Your script (`.R` file) — code you write and save |
| **Console** (bottom-left) | Where code actually runs, line by line |
| **Environment** (top-right) | Variables and objects you've created |
| **Files / Plots / Packages / Help** (bottom-right) | File browser, plot viewer, package manager, documentation |

Open the script `01-intro-to-R.R` that comes with this practical. We'll run it together — line by line, using `Ctrl+Enter` (Windows/Linux) or `Cmd+Enter` (Mac) to send a line from the Source panel to the Console.

---

## 1. R basics 

In R, you type commands line by line. You can also assign results to objects using the arrow <-.
For example:

```r
a <- 3a
```

This assigns the value 3 to the variable a and prints it.

### 1.1 R as a calculator and assignment

```r
2 + 2
10 / 3
2^5          # exponent

# Assign values to variables with <-  (preferred) or =
x <- 5
y <- 10
x + y
```

### 1.2 Data types you'll meet constantly

```r
a <- 4.5          # numeric
b <- "Homo sapiens"  # character (string)
d <- TRUE         # logical
e <- c(1, 2, 3, 4)   # vector of numbers
f <- c("ape", "ggtree", "phytools")  # vector of strings

class(a)
class(f)
length(e)
```

`c()` ("combine") is the function you'll use constantly to build vectors — e.g., a list of species names.

### 1.3 Vectors and indexing

```r
species <- c("Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Gallus_gallus")
species[1]        # first element
species[2:3]      # second and third
length(species)
```

### 1.4 Data frames (your main data structure)

A data frame is a table — rows and columns, like a spreadsheet. Most trait/comparative data you'll use is a data frame.

```r
traits <- data.frame(
  species = c("Homo_sapiens", "Pan_troglodytes", "Mus_musculus"),
  body_mass_kg = c(70, 50, 0.02),
  genome_size_Mb = c(3100, 3300, 2700)
)

traits            # view it
str(traits)       # structure: column types
traits$body_mass_kg   # access one column with $
nrow(traits)
```

### 1.5 Functions and getting help

```r
mean(traits$body_mass_kg)
sqrt(16)

?mean        # opens help page
help(mean)   # same thing
```

### 1.6 Working directory and projects

R needs to know where to look for files. Always check this first:

```r
getwd()                 # where am I?
# setwd("path/to/folder")  # change it if needed
```

**Better habit:** use an RStudio **Project** (`File > New Project`): it sets the working directory for you automatically and keeps everything for the course in one folder.

---

## 2. Installing packages 

R's base installation is small on purpose. Almost everything useful comes from **packages**, and for this course, that means mainly phylogenetics packages.

There are two main repositories:

- **CRAN**: the standard R package repository. Install with `install.packages()`.
- **Bioconductor**: repository for bioinformatics/genomics packages (like `ggtree`). Needs a small setup step first.

### 2.1 Installing from CRAN

```r
install.packages("ape")
install.packages("phytools")
install.packages("ggplot2")
install.packages("BAMMtools")
```

You only need to **install** a package *once per computer*. Run these lines now if you haven't.

### 2.2 Installing from Bioconductor

`ggtree` lives on Bioconductor, not CRAN, so it needs the `BiocManager` bridge:

```r
install.packages("BiocManager")   # one-time setup
BiocManager::install("ggtree")

# There is another way to do it:
if (!requireNamespace("BiocManager", quietly = TRUE))  
	install.packages("BiocManager")
BiocManager::install("ggtree", update = TRUE, ask = FALSE)
```

If you're asked `Update all/some/none? [a/s/n]`, typing `n` (none) is fine and faster for class purposes.

### 2.3 Loading packages

Installing puts the package on your computer. **Loading** with `library()` makes it available in your current session (you need to do this every time you start a new R session).

```r
library(ape)
library(ggtree)
library(phytools)
library(ggplot2)
```

Packages in R are collections of functions.

**Example**: `ape` is used to manipulate phylogenetic trees.

**Example**: `ggtree` is used to visualise them.

If any of these throws an error, re-run the matching `install.packages()` / `BiocManager::install()` line above and try again.

### 2.4 The packages we'll use this course

| Package | What it's for | Source |
|---|---|---|
| `ape` | Core package: reading/writing/plotting trees, the `phylo` object format everything else builds on | CRAN |
| `ggtree` | Publication-quality tree plotting using `ggplot2`'s "grammar of graphics" | Bioconductor |
| `phytools` | Comparative methods — ancestral state reconstruction, trait evolution models, tree manipulation | CRAN |
| `ggplot2` | General-purpose plotting, used under the hood by `ggtree` | CRAN |
| `BAMMTools` | Analysis and Visualization of Macroevolutionary Dynamics on Phylogenetic Trees | CRAN |

You'll probably meet a few more specialised packages (e.g. `phangorn`, `treeio`, `caper`) later in the course as needed (same install pattern applies).

---

## 3. Quick troubleshooting tips

- **"there is no package called X"** → you forgot to install it. Run `install.packages("X")`.
- **"could not find function X"** → installed but not loaded. Run `library(X)`.
- **Package install fails with a compilation error** → try `install.packages("X", type = "binary")`.
- Stuck? Select the broken line and ask: *"what does this error mean?"* — but also just try Googling the exact error message, it's a real R skill.

---


## 4. Working with phylogenetic trees 

A phylogenetic tree in R is just an object of class `phylo` (defined by `ape`). Let's simulate one and plot it — no biology theory needed yet, just to see what we'll be working with for the rest of the course.

```r
library(ape)

# Simulate a random tree with 10 tips ("species")
set.seed(123)
tree <- rtree(10)

# What is this object?
class(tree)
tree
str(tree)

# Plot it - base R / ape style
plot(tree)
```

**Notice**: same underlying `phylo` object, two different plotting philosophies. `ape::plot.phylo()` is quick and built-in; `ggtree()` gives you the full `ggplot2` layering system (colors, themes, combining with data) which we'll lean on heavily later in the course.

<img width="1341" height="762" alt="image" src="https://github.com/user-attachments/assets/267bef6c-839d-4b01-bdb5-f72fcb11530e" />

## 5. Questions
1. What do you see? Describe the phylogeny in as much detail as possible.
2. How many tips does the tree have?
3. Can you identify branches and nodes?
---

## 6. Reading a tree from a file
In bioinformatics, trees are often stored in Newick format. We can create a simple example:

```r
tree <- read.tree(text="((A:1,B:1):2,C:3);")
plot(tree)
```
In this tree:
<img width="1275" height="727" alt="image" src="https://github.com/user-attachments/assets/079f7f21-3f0b-4ae5-90d7-09e013a92c7c" />


- A, B, and C are taxa
- numbers represent branch lengths

### 6.1 Exploring the tree structure
Let's inspect the object:

```r
tree
```

You should see information about:
- number of tips
- number of nodes
- branch lengths
- labels
- root

### 6.2 Basic operations on trees
#### 6.2.1 Distances between taxa
We can compute pairwise distances:

```r
d <- cophenetic(tree)
d
```
Example output:
```r
  A B C
A 0 2 6
B 2 0 6
C 6 6 0
```

*These values represent evolutionary distances between taxa.*

#### 6.2.2 Removing taxa
Sometimes we want to prune a tree:

```r
pruned <- drop.tip(tree, "A")
plot(pruned)
```
<img width="1258" height="637" alt="image" src="https://github.com/user-attachments/assets/a1fdf793-d5e6-4376-a5c0-360bd14e8b0b" />

4. What changed compared to the original tree?

#### 6.2.3 Rooting a tree
We can root the tree using a specific taxon:

```r
rooted <- root(tree, "C")
plot(rooted)
```
<img width="1236" height="592" alt="image" src="https://github.com/user-attachments/assets/6b360081-6855-481d-9e0f-2a4a66209e61" />

*Rooting changes the interpretation of evolutionary relationships.*

#### 6.2.4 Exercise
Remove two species from a tree and compare the resulting topology.

### 7. Visualising trees with ggtree
* The `plot()` function is simple but limited.
* The `ggtree` package allows more flexible and publication-quality figures.

#### 7.1 Basic visualisation

```r
ggtree(tree)
```
<img width="1380" height="572" alt="image" src="https://github.com/user-attachments/assets/19b7ceaa-f69b-4d91-bc16-fbe89921eea6" />

*This creates a tree using the grammar of graphics.*

#### 7.2 Adding tip labels:

```r
ggtree(tree) +  
geom_tiplab(size = 4)
```
<img width="1307" height="570" alt="image" src="https://github.com/user-attachments/assets/3f967b4b-8409-4f36-b99e-c61de50b4bdd" />

*This displays the names of the taxa.*

#### 7.3 Adding nodes and aesthetics

```r
ggtree(tree) +  
geom_tiplab() +  
geom_nodepoint(color = "red") +  
theme_tree()
```
<img width="1292" height="510" alt="image" src="https://github.com/user-attachments/assets/4a596607-80c7-48b1-b33d-3e2c726634d4" />

5. What do the red points represent?

#### 7.4 Changing layouts
We can visualise the same tree with different layouts:

```r
ggtree(tree, layout = "circular") +  
geom_tiplab()
```
<img width="472" height="441" alt="image" src="https://github.com/user-attachments/assets/a58578cf-9de2-4284-8d0b-13e0ed951d4d" />

*Other layouts are available (rectangular, slanted, etc.).*

#### 7.5 Exercise
Create a circular tree with:
- tip labels
- highlighted nodes


## Useful tips
In R, you can get help using:

```r
?rtree
?ggtree
```

## Errors are normal, read them carefully and try to understand what is missing.
Always check your objects:

```r
class(tree)
```
---

## Recap
This tutorial provides only the basic tools to work with phylogenetic trees in R. During the course, you will combine these skills with real genomic datasets and more complex analyses.

By now you should have:

- [ ] Run basic R commands (arithmetic, variables, vectors, data frames)
- [ ] Installed `ape`, `phytools`, `ggplot2`, `BAMMTools` (CRAN) and `ggtree` (Bioconductor)
- [ ] Loaded all four with `library()` without errors
- [ ] Simulated and plotted your first `phylo` tree object

Good luck and have fun! 
