# Species Delimitation & Networks

  ## 1. Introduction
  The Southern Ocean nudibranch genus Tritonia has a messy taxonomic history: several named "species" turned out to be colour morphs of the same biological species, and at least one specimen had been misidentified for   decades. Rossi et al. (2021) used three markers  (*COI* (mitochondrial, protein-coding), *16S rRNA* (mitochondrial, non-coding), and *H3* (nuclear, protein-coding)), together with two species delimitation methods (*ABGD* and *GMYC*) to show that Antarctic/Weddell Sea specimens and Bouvet Island specimens form two clearly separated species: *T. challengeriana* and *T. dantarti*, regardless of their orange or white colouration.

The molecular and morphological/colour stories disagree, which makes species delimitation results meaningful and discussable rather than a foregone conclusion. It pairs naturally with the mollusk datasets used in your other practicals. A follow-up study (Rossi et al. 2023, Cladistics; see 'further_reading' section) applied mPTP directly to this group and even resurrected an old genus name (*Myrella*) based on molecular species delimitation.

  ### Objectives
    By the end of this session, you should be able to:


  * Explain what a single-locus species delimitation method does, and what it does not do.
  * Run ASAP and mPTP on real data and interpret their output.
  * Build a quick ML tree as input for tree-based delimitation (mPTP).
  * Build a phylogenetic/haplotype network and explain how it differs from a tree.
  * Critically compare tree and network representations of the same data, and identify when a network is showing you something a tree cannot.

  ### Datasets and software
  The original paper analyses 41 *Tritonia* ingroup specimens + 17 outgroup taxa (58 sequences) across three markers. We use here a reduced dataset: the full *Tritonia* ingroup, plus 2–3 outgroups only (one close relative for rooting, e.g. another Tritoniidae/Dendronotoidea species, rather than all four Proctonotoidea rooting taxa). 

  For this session you'll need to use software pre-installed on your local computer:
  * PopArt (https://popart.maths.otago.ac.nz/download/)
  * SplitsTree (https://software-ab.cs.uni-tuebingen.de/download/splitstree6/welcome.html)

  Also software hosted in online servers:
  * mPTP ([https://link.springer.com/article/10.1007/s00300-021-02813-8](https://mptp.h-its.org/#/tree))
  * ASAP
  * GMYC (https://species.h-its.org/gmyc/)
  * 


## 2. BAMM analysis

There are many ways to dected speciation events or diversification in evolutionary rates. For this practical we are going to start with a BAMM analysis.


```R
##BAMM analysis
library(ape)
library(phytools)

#First we make the tree ultrametric. This is just for this practical, normally you will have a time calibrated topology, that is already ultrametric
# Load tree
tree <- read.tree("coi_sequences.fasta.mafft.tree")

# Basic checks
is.rooted(tree)
is.binary.tree(tree)
is.ultrametric(tree) #Bamm wants an ultrametric tree, if returned FALSE we need to convert it in ultrametric

# Plot initial tree
plot(tree, cex=0.4)
title("Original COI tree (non-ultrametric)")

# Convert to ultrametric using penalized likelihood
# (chronos = fast relaxed clock method)
chrono_tree <- chronos(tree)

# Check ultrametricity
is.ultrametric(chrono_tree)

# Compare before/after
par(mfrow=c(1,2))
plot(tree, cex=0.4)
title("Original tree")

plot(chrono_tree, cex=0.4)
title("Chronos ultrametric tree")

# Optional: rescale root age (e.g. 100 time units)
chrono_tree$edge.length <-
  chrono_tree$edge.length *
  (100 / max(node.depth.edgelength(chrono_tree)))

# Final validation
is.ultrametric(chrono_tree)

# Save tree # This is the tree we are going to use as input for the BAMM
write.tree(chrono_tree, file="nudibranch_ultrametric_chronos.tre")

# BAMM PIPELINE: Diversification rates in nudibranchs
# Requires ultrametric tree from chronos or BEAST

library(ape)
library(BAMMtools)
library(coda)

# Load ultrametric tree
tree <- read.tree("nudibranch_ultrametric_chronos.tre")
#root tree
rooted_tree <- root(
  tree,
  outgroup = c(
    "DQ026831.1_Dirona_picta_cytochrome_oxidase_subunit_1-like",
    "KF643932.1_Dirona_albolineata_voucher_10BCMOL-00344_cytochrome_oxidase_subunit_1"
  ),
  resolve.root = TRUE
)
# Check tree properties
is.ultrametric(rooted_tree)
is.binary.tree(rooted_tree)

plot(rooted_tree, cex=0.4)
title("Ultrametric tree for BAMM")

#Calculate the sampling fraction
#Based on species we have in our dataset, we need to calculate the sampling fractions (i.e. ratio species in tree/species described)
#Remember the BAMM analysis needs to be aware of the taxon sampling!!
sampling_fractions <- c(
  Dirona = 2/3,
  Tritonia = 8/20,
  Myrella = 3/4,
  Tritonicula = 2/5,
  Candiella = 2/8,
  Tritoniella = 1/2,
  Marionia = 5/20,
  Bornella = 5/10,
  Leminda = 1/1,
  Charcotia = 1/2
)

# Global sampling fraction (approx)
global_sampling <- mean(sampling_fractions)
writeLines(paste(global_sampling), con = "sampling_fractions.txt")


# Set the BAMM priors
setBAMMpriors(tree,outfile = "bamm_priors.txt")
##You wiil need to paste the outfile in a specific section of the control file
## Now we move on the terminal, to run the BAMM. We'll be back in R to analyse the output
```



```R
# Load BAMM results
library(coda)
#Check convergence
mcmc <- read.csv("Mollusca_bamm_mcmc_out.txt")

plot(mcmc$logLik, type="l",
     main="MCMC trace: log-likelihood")
#Remove burnin
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
#Check lolik before and after removing burnin
effectiveSize(mcmc$logLik)
effectiveSize(postburn$logLik)
effectiveSize(postburn$N_shifts)

#Plot diversification rates
# Phylogenetic rate map
plot.bammdata(edata,lwd = 2, spex = "s")
addBAMMshifts(edata)

#Identify rate shifts
edata <- getEventData(rooted_tree, eventdata = "event_data.txt", burnin=0.25, type = "diversification")
shift_probs <- summary(edata)
computeBayesFactors(mcmcout, expectedNumberOfShifts=1, burnin=0.25)

#Plot the shift identified by the BAMM
q <- plot.bammdata(edata, legend=TRUE, lwd=2, method="phylogram", pal="BrBG", breaksmethod='jenks')
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
summary(css)
plot.credibleshiftset(css, lwd=1.5, pal="BrBG")##compute the marginal shift probabilities

# Plot the diversification through time

plotRateThroughTime(edata)
plotRateThroughTime(edata, rate = "speciation")
plotRateThroughTime(edata, rate = "extinction")
plotRateThroughTime(edata, rate = "netdiv")


##plot a mean phylorate plot depicting net diversification rates
q2 <- plot.bammdata(edata, legend=TRUE, lwd=2, method="phylogram", pal="BrBG", breaksmethod='jenks', spex = "netdiv")
##plot histogram of rates
ratesHistogram(q, plotBrks = TRUE, xlab = 'speciation rates')
ratesHistogram(q2, plotBrks = TRUE, xlab = 'net diversification rates')




```