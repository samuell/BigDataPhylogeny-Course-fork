# Gene Family Expansion & Contraction with CAFE5

## 0. Introduction

In this practical we will explore our dataset at the phylum-level, searching for the gene families *expanded* or *contracted* along branches, exploring the rate at which gene content changes.

*CAFE5* models gene family size as a birth-death process along a phylogeny: at every point on every branch, a gene family can gain a copy ("birth") or lose one ("death") at some rate. The single parameter *λ (lambda)* is the rate of this process (how fast gene family sizes tend to change per unit branch length, genome-wide).

Given:
* a tree with branch lengths in comparable (ultrametric) time units, and
* a table of how many gene copies each species has, for each gene family,

*CAFE5* finds the value of λ that makes the observed gene counts most probable under the birth-death model (maximum likelihood), then uses that fitted model to:
* Reconstruct the most likely ancestral gene count at every internal node of the tree, for every family.
* Flag gene families whose size changes are statistically unlikely under the genome-wide λ (i.e., families evolving unusually fast on specific branches). These are your expansion/contraction candidates.

### Objectives

By the end of this practical you'll be able to:
* Understand what a birth-death model of gene family evolution is estimating (λ), and what assumptions it makes (e.g. why CAFE needs an ultrametric, time-calibrated tree, not just a topology).
* Convert real-world bioinformatics output (OrthoFinder) into CAFE5's required input format, a routine but essential reformatting skill.
* Run CAFE5 and interpret its core output: the global λ, and which gene families are significantly rapidly evolving.
* Identify, for at least one gene family, on which branch(es) of the tree expansion or contraction is inferred to have happened.
* Map CAFE5's results back onto the phylogeny and produce a figure of gene family size evolution.

### Datasets and software

For this practical you'll need to use the following software (already installed in the server):
* [CAFE v5.1.0.](https://anaconda.org/channels/bioconda/packages/cafe/overview)
* [biopython](https://anaconda.org/channels/conda-forge/packages/biopython/overview)
* R, with the packages `ape`, `ggtree`, `ggplot2`, `stringr`, and `tidytree` (for Section 2, visualisation)

We'll be using the following datasets:
* `cafe_input.tsv`: this corresponds to the 'Orthogroups.GeneCount.tsv' reformatted for CAFE5 (we already prepared it for you).
* `cafe_tree_ultrametric.nwk`: the topology from previous practicals, made ultrametric. We use `ape::chronos()` (penalized likelihood rate-smoothing) to produce a relative-time ultrametric tree directly from branch lengths.

Our tree covers 20 molluscan (+1 brachiopod outgroup) species:

```
((((((((S_philippinensis,V_penis)N8,S_velum)N7,(P_vernedei,S_dalli)N9)N6,((A_californica,C_concholepas)N11,((T_virginea,L_scabra)N13,H_cracherodii)N12)N10)N5,((S_atlantica,O_vulgaris)N15,N_pompilius)N14)N4,V_oligotropha)N3,(((A_discrepans,C_septemvalvis)N18,D_sirenkoi)N17,(E_babai,N_megatrapezata)N19)N16)N2,L_anatina)N1;
```

Note that the internal nodes are already labeled (`N1`–`N19`). **Keep these labels** because they are how you will cross-reference CAFE5's ancestral state reconstruction back onto specific branches later in Section 2.

---

### Why does CAFE5 need an *ultrametric* tree? (and not just any tree with branch lengths)

A regular phylogeny's branch lengths usually represent **amount of evolutionary change** (e.g., substitutions per site). Two sister species can sit at very different distances from their common ancestor if one lineage evolved faster than the other at the sequence level. That's normal and expected for a substitution-based tree.

CAFE5, however, is not modelling sequence substitutions. It's modelling **gene gain/loss as a function of time**, because birth and death are processes that happen *per unit of time*, not per unit of sequence divergence. For the model to compare rates across branches meaningfully, every tip needs to be the same temporal distance from the root (i.e., all species are sampled "now," and the branch lengths must all be in units of **time elapsed**, not substitutions accumulated.

A tree is **ultrametric** when the total branch length from the root to *every* tip is equal. Think of it like this:

> Take a real family tree of living relatives. Whatever happened in each branch of the family, all the living descendants are alive *today* (they are all the same "distance in time" from their common ancestor), even if their genealogies are more or less eventful in between. That's ultrametricity: the *tips* (today) are level, even though the topology and the events inside each branch can be totally different.

If you fed CAFE5 a tree with substitution-based branch lengths instead, you'd be implicitly telling it that a fast-evolving lineage (long branch, lots of substitutions) has had *more time* for genes to be gained/lost than a slow-evolving one, which conflates two unrelated things: rate of sequence evolution and elapsed time. This would systematically bias your estimate of λ and could make a perfectly normal lineage look like a gene-family-evolution outlier just because its DNA happens to substitute quickly.

This is exactly why we ran the topology through `ape::chronos()` beforehand: chronos takes a tree with substitution-based branch lengths and, using penalised likelihood rate-smoothing, finds a set of *relative time* branch lengths that (a) are internally consistent with an ultrametric tree, and (b) deviate as little as is statistically defensible from the original branch lengths, smoothing out local rate variation rather than blindly forcing equal branch lengths everywhere.

You'll find more information in this [Tutorial](https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md).

**Quick self-check before moving on:** if you `cat` the tree file and you see that all root-to-tip path lengths sum to (approximately) the same number, you have an ultrametric tree. If you're ever unsure whether a tree you've been handed is ultrametric, in R you can check directly:

```r
library(ape)
tree <- read.tree("data/cafe_tree_ultrametric.nwk")
is.ultrametric(tree)  # should return TRUE
```

## 1. Run CAFE5

### 1.1 Check your input files

```bash
ls -la data/cafe_input.tsv data/cafe_tree_ultrametric.nwk
head -3 data/cafe_input.tsv
cat data/cafe_tree_ultrametric.nwk
```

**Questions**:
- How many gene families are in `cafe_input.tsv`?
- How many species columns?
- Does the species count match the number of tips in the tree?

> *Tip: `wc -l data/cafe_input.tsv` gives you the family count (minus 1 for the header). `head -1 data/cafe_input.tsv | tr '\t' '\n' | wc -l` gives you the column count, subtract the non-species columns (typically "Orthogroup"/family ID, and sometimes a "Desc" column) to get the species count. Compare that to the 20 tips in our tree above.*

### 1.2 Running CAFE5: single global λ

This is the simplest CAFE5 model: one rate of gene family evolution for the entire tree, no exceptions. Don't forget to activate the environment.

```bash
mkdir -p results/cafe_base
cafe5 -i data/cafe_input.tsv -t data/cafe_tree_ultrametric.nwk -o results/cafe_base
```

*This should run in well under a minute for a dataset this size. Watch the terminal: CAFE5 prints its optimisation progress as it searches for the best-fitting λ.*

<img width="356" height="282" alt="image" src="https://github.com/user-attachments/assets/2085c048-75ec-4aa6-a128-4539aca586fb" />

**Checkpoint**: Open `results/cafe_base/Base_results.txt`:

```bash
cat results/cafe_base/Base_results.txt
```

You should see something like:

```
Model Base Result: <some log-likelihood score>
Lambda: <a small number, e.g. 0.001-0.01 range typical>
```

#### Interpreting λ

λ is a **rate per unit branch length**. In our relative-time tree, think of it as "expected number of gene gain/loss events per gene family, per unit of relative time, per lineage." A handful of points to anchor your intuition:

- λ on its own has no fixed "good" or "bad" value, it's only interpretable *relative to* the branch lengths in your specific tree (which are relative-time units from `chronos()`, not millions of years), and relative to λ estimated in other comparable studies using the same kind of tree.
- A **single global λ** is a strong simplifying assumption: it says every branch in the tree, from your fastest-evolving lineage to your slowest, shares the *same* underlying rate of gene family turnover. Real genomes rarely cooperate with this assumption uniformly. Some lineages undergo bursts of gene family expansion (e.g., after whole-genome duplication, or in venom/immune gene families under strong selection) that a single genome-wide rate will smooth over.
- The global λ tells you the *average* behaviour, but the per-family and per-branch tests are how you find the *exceptions* to that average, which are usually the biologically interesting bit.
- If you want to go further after this practical, CAFE5 also supports a **multi-λ ("k") model**, where you assign different parts of the tree (e.g., different clades) to different rate categories, and an **error model**, which accounts for genome assembly/annotation error inflating or deflating apparent gene counts. Both are beyond today's scope, but worth knowing they exist. The single global λ model here is a deliberately simple starting point, not the end of the analysis.

### 1.3 What other files did CAFE5 produce?

```bash
ls results/cafe_base/
```

You should see several files; the ones that matter most for this practical:

| File (name may vary slightly) | What it contains |
|---|---|
| `*_results.txt` | The headline: model score and estimated λ |
| `*_asr.tre` | **A**ncestral **s**tate **r**econstruction: a tree (Nexus format) where every node, including internal/ancestral ones, has an inferred gene copy number attached |
| `*_clade_results.txt` | Per-branch summary: how many gene families expanded vs. contracted on that branch |
| `*_family_results.txt` | Per-gene-family table with significance values (see below) |

### 1.4 Find the gene families CAFE5 thinks are interesting

Look at the per-family results table (`*_family_results.txt`). You're looking for two columns:

- **Family-wide P-value**: is this family's overall pattern of size change unlikely under the fitted genome-wide λ?
  (Low p-value, conventionally < 0.01, = yes, interesting.)


### 1.5 Questions

1. How many gene families come out as significant (family-wide p < 0.01)?
2. As a fraction of all families tested, does that number seem high, low, or about what you'd expect?
   (There's no universally "correct" answer here, it depends entirely on the data, but having an intuition for the proportion matters more than memorising a threshold.)

---

## 2. Visualising gene family evolution on the tree

So far everything has lived in text files. The real payoff of CAFE5 is seeing *where on the tree* expansions and contractions happened. This section walks you from a single significant family to a publication-style annotated phylogeny.

### 2.1 Pick a family to visualise

From your Section 1.4–1.5 exploration, pick one family with a low family-wide p-value (ideally one with a clearly significant p-value). Note its Orthogroup ID (e.g. `OG0001174`); you'll need it later.

```bash
grep "OG0001174" results/cafe_base/Base_family_results.txt
```

### 2.2 Extract the ancestral states for your family

The `*_asr.tre` file contains, for **every** family, a tree annotated with inferred copy number at every node (tip and internal). It's a multi-tree Nexus file — one tree per family, in the same order as your input table. You need to pull out the one tree that corresponds to your chosen family.

```bash
grep -n "OG0001174" results/cafe_base/Base_asr.tre
```

This should point you to the right line. Each tree in this file looks like a normal Newick/Nexus tree, except every node label carries an extra `_<count>` suffix (e.g. `N8_3` means "node N8 was reconstructed to have 3 gene copies"). 

### 2.3 Plot it on the phylogeny with `ggtree` (R)

`ggtree` is the field-standard tool for exactly this kind of "decorate a phylogeny with extra data" figure, and it plays nicely with CAFE5/CAFE-family output. Below is a minimal, adaptable script. Copy it into an R script or RStudio and adjust the family ID and counts to match your own data.

```r
library(ape)
library(ggtree)
library(ggplot2)
library(tidytree)
library(stringr)

# 1. Read in the base, ultrametric tree (topology + node labels we kept)
tree <- read.tree("cafe_tree_ultrametric.nwk")

# 2. Build a small data frame of ancestral + tip copy numbers for YOUR chosen family,
#    read off from Base_asr.tre (Section 2.1). Replace these example values with
#    the real counts for your family. To do so, you can extract the information from the asr_tree 
#    file:

tree_text <- readLines("results/cafe_base/Base_asr.tre")
og <- tree_text[grep("OG0001174", tree_text)]

# 3. Copy tips
tips <- str_extract_all(og, "[A-Za-z_]+<\\d+>_\\d+")[[1]]
tip_df <- data.frame(
  label = str_extract(tips, "^[A-Za-z_]+"),
  copies = as.numeric(str_extract(tips, "(?<=_)\\d+$"))
)

# 4. Copy values at nodes
nodes <- str_extract_all(og, "<\\d+>_\\d+")[[1]]

node_df <- data.frame(
  label = paste0("N", seq_along(nodes)),
  copies = as.numeric(str_extract(nodes, "(?<=_)\\d+$"))
)

copy_numbers <- rbind(tip_df, node_df)

# 5. Attach the data to the tree object
p <- ggtree(tree) %<+% copy_numbers

# 6. Plot: tip labels, tip+node points sized/colored by copy number, node labels for reference
p +
  geom_tiplab(size = 3, offset = 0.02) +
  geom_point(aes(color = copies, size = copies)) +
  scale_color_viridis_c(name = "Gene copies") +
  scale_size_continuous(name = "Gene copies", range = c(2, 8)) +
  theme_tree2() +
  ggtitle("OG0001174 — ancestral gene copy number reconstruction")
```
You'll see something like this:
<img width="1446" height="851" alt="image" src="https://github.com/user-attachments/assets/630098c7-1cf5-4bc0-84be-28b097273ada" />


**Some tips:**
- `%<+%` is `ggtree`'s operator for attaching a metadata data frame to a tree object by matching the `label` column.
- Coloring/sizing points by copy number turns "family-wide p-value is significant" into something visual: you can immediately see *where* the copy number jumps happen along the tree, which is the entire point of doing ancestral state reconstruction in the first place.

### 2.4 Highlighting the branch(es) CAFE5 flagged as significant

Once you've identified, from Section 1.4, *which branch* had a significant p-value for your chosen family, you can highlight that branch directly rather than relying on the reader to compare copy numbers visually:

```r
# Suppose Section 1.4 told you the significant branch leads to node N16
# (i.e., the contraction/expansion happened on the branch ABOVE N16)
p +
  geom_tiplab(size = 3, offset = 0.02) +
  geom_point(aes(color = copies, size = copies)) +
  geom_hilight(node = MRCA(tree, c("D_sirenkoi", "E_babai")), fill = "firebrick", alpha = 0.15) +
  scale_color_viridis_c(name = "Gene\ncopies") +
  scale_size_continuous(name = "Gene\ncopies", range = c(2, 8)) +
  theme_tree2() +
  ggtitle("OG0001174 — highlighted branch shows significant p-value")
```

`geom_hilight()` shades the clade descending from a given node. Here we find that node with `MRCA()` using two tips that define it, rather than hardcoding a node number, since internal `ggtree`-assigned node numbers don't always match your original `N1`–`N19` labels one-to-one (`ggtree` numbers nodes internally for plotting purposes, separately from any labels stored in the original tree file).

<img width="1431" height="792" alt="image" src="https://github.com/user-attachments/assets/630aa419-8a83-4991-98a8-5c679ddf0cf7" />


### 2.5 A tree-wide expansion/contraction summary

For a "big picture" figure across *all* families rather than just one, `*_clade_results.txt` (Section 1.3) gives you, per branch, a count of how many families expanded vs. contracted. This is naturally a **bar plot per branch**, or, more strikingly, a tree where each branch is colored by net expansion/contraction:

```r
library(ape)
library(ggtree)
library(ggplot2)

# 1. Read data
change <- read.table(
  "results/cafe_base/Base_change.tab",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = ""
)

tree <- read.tree("species_tree_topology.nwk")

# 2. Net change per node
node_cols <- grep("<", colnames(change), value = TRUE)
net_mat <- change[, node_cols]
names(net_mat) <- node_cols
net_change <- colSums(net_mat, na.rm = TRUE)

# 3. Select nodes from the tree
p <- ggtree(tree)

tree_df <- p$data %>%
  mutate(
    cafe_node = ifelse(isTip, label, paste0("<", node, ">"))
  )

tree_df$net_change <- net_change[tree_df$cafe_node]

# 4. Remove NAs
tree_df_clean <- tree_df %>%
  filter(!is.na(net_change))

# 5. Define extreme clades
thr <- quantile(abs(tree_df_clean$net_change), 0.95)

tree_df$highlight <- ifelse(
  !is.na(tree_df$net_change) & abs(tree_df$net_change) > thr,
  TRUE,
  FALSE
)

# 6. Final plot
ggtree(tree) %<+% tree_df +
  geom_tree(aes(color = net_change), linewidth = 1) +
  geom_tiplab(size = 3) +
  scale_color_gradient2(
    low = "firebrick",
    mid = "grey90",
    high = "steelblue",
    midpoint = 0,
    name = "Gene family\nexpansion / contraction"
  ) +
  geom_point(aes(size = abs(net_change), color = net_change)) +
  scale_size_continuous(
    name = "CAFE5 inferred\nfamily size shift"
  ) +
  theme_tree2() +
  ggtitle("CAFE5 gene family expansions/contractions")
```

<img width="1397" height="702" alt="image" src="https://github.com/user-attachments/assets/71afccce-b8c7-442a-8d2e-33e15de13d15" />


This is the figure that best answers the high-level question this whole practical is built around: **which lineages have been gaining gene families fastest, and which have been losing them?**
This figure shows how gene families have changed throughout the evolution of animals in your phylogenetic tree, estimated with CAFE5.
Specifically, we are looking at gene families **expansions** and **contractions** in different evolutionary lineages. The colors of the branches correspond to: blue (or light gray) → loss of gene families (contractions) and red → gain of gene families (expansions). The color indicates the estimated global change at that node/branch. The size of the dots indicates the magnitude of the evolutionary change. Small circle → few changes, large circle → many changes of gene families. **Important**: it is not the number of individual genes, but entire families.

We talk about gene family expansion when a gene family duplicates and more copies appear. It can indicate functional adaptation, new biological functions, specialisation, etc. Contraction of gene families, on the other hand, indicates that genes are lost and the family becomes smaller. It can indicate genomic simplification, loss of unnecessary functions, or adaptation to new environments. What we have done with CAFE5 here is to model the evolution of the number of genes per family and reconstruct what the ancestral genomes were like where genes have been gained or lost. It is observed that some lineages have experienced more genomic dynamism (expansions/contractions) than others. 
This particular tree shows how the number of gene families has changed throughout evolution. Red indicates gene loss, blue indicates gain, and the size of the dots indicates how many changes there have been. This allows us to see which lineages have had the most genomic innovation.

### 2.6 Questions

1. For the family you chose in Section 2.1, which branch shows the largest jump in inferred copy number? Does this match the branch CAFE5 flagged with a significant p-value?
2. Looking at your tree-wide summary: are expansions/contractions spread evenly across the tree, or concentrated in a few lineages? If concentrated, can you think of a biological hypothesis (ecology, life history, known whole-genome or tandem duplications in that group) that might explain it?

## 3. Wrap-up

You've now gone from raw OrthoFinder gene counts, through a maximum-likelihood birth-death model fitted with CAFE5, to ancestral state reconstructions and a tree-mapped visualisation of where gene family content has expanded or contracted across your phylogeny. The single global λ model you ran here is the simplest version of this analysis. In your own research, the natural next steps are exploring multi-rate models (different λ for different clades) and incorporating an error model to account for genome annotation noise, both built into CAFE5 but outside today's scope.
