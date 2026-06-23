# Phylogeography, Species Delimitation, Networks, Hybridisation & Diversification Rates


  ## 0. Introduction
  This practical session covers the conceptual and computational tools needed to move beyond simple phylogeny reconstruction. You will work with two complementary datasets: (1) real sequence data from the Antarctic nudibranch genus *Tritonia* ([Rossi et al. 2021](https://link.springer.com/article/10.1007/s00300-021-02813-8), Polar Biology) for species delimitation and network analyses; and (2) a simulated mollusc-like dataset to explore hybridisation, coalescent methods, and diversification-rate analysis. The BAMM section uses the COI tree from the GitHub repository.

  The Southern Ocean nudibranch genus *Tritonia* has a messy taxonomic history: several named "species" turned out to be colour morphs of the same biological species, and at least one specimen had been misidentified for   decades. Rossi et al. (2021) used three markers  (*COI* (mitochondrial, protein-coding), *16S rRNA* (mitochondrial, non-coding), and *H3* (nuclear, protein-coding)), together with two species delimitation methods (*ABGD* and *GMYC*) to show that Antarctic/Weddell Sea specimens and Bouvet Island specimens form two clearly separated species: *T. challengeriana* and *T. dantarti*, regardless of their orange or white colouration.

The molecular and morphological/colour stories disagree, which makes species delimitation results meaningful and discussable rather than a foregone conclusion. It pairs naturally with the mollusc datasets used in your other practicals. 

  ### Objectives
  
By the end of this session, you should be able to:

  * Interpret geographic patterns in phylogenetic trees (phylogeography).
  * Explain what a single-locus species delimitation method does, and what it does not do.
  * Run ASAP and mPTP species delimitation on real data and critically evaluate their output.
  * Build and interpret phylogenetic/haplotype networks and explain how it differs from a tree.
  * Visualise reticulate evolution in SplitsTree and articulate when a network outperforms a tree.
  * Understand how hybridisation, introgression and recombination produce non-tree-like signals.
  * Run ASTRAL and ASTRAL-Pro to estimate species trees from gene trees (coalescent framework).
  * Interpret BAMM and MEDUSA diversification-rate analyses.

---

## 1. Simple Phylogeographic Example: Interpreting Geographic Patterns in Trees
  ### 1.1 Background
  The nudibranch genus *Tritonia* provides a textbook case of how molecular phylogenetics resolves taxonomic puzzles that morphology alone cannot. [Rossi et al. (2021)](https://link.springer.com/article/10.1007/s00300-021-02813-8) collected specimens from the Weddell Sea (high Antarctic, ~70–75°S) and Bouvet Island (Sub-Antarctic, ~54°S), covering depths from 130 to 789 m. Both sites harbour orange and white colour morphs, which were previously considered distinct species or at least diagnostically useful characters.

**Key finding**: COI + 16S + H3 data placed Weddell Sea specimens in one strongly supported clade (PP = 1, BS = 100) and Bouvet Island specimens in another (PP = 1, BS = 98), regardless of colour. Two species were supported by both ABGD and GMYC: *T. challengeriana* (circum-Antarctic) and *T. dantarti* (Bouvet Island endemic).

  ### 1.2 Dataset
  You have downloaded the data from Rossi et al. (2021). The dataset comprises 41 *Tritonia* ingroup specimens + 2–3 outgroup taxa, sequenced for three markers:
•	COI — mitochondrial, protein-coding, ~601 bp (3rd codon position included)
•	16S rRNA — mitochondrial, non-coding, ~486 bp
•	H3 — nuclear, protein-coding, ~328 bp (3rd codon position included)

  ### 1.3 Phylogeographic interpretation
  Open the published ML/BI tree (Fig. 2 of Rossi et al. 2021, or load your own RAxML run output). Answer the following:
1.	Specimens are colour-coded: red = Bouvet Island, blue = Weddell Sea. Do the two colours form reciprocally monophyletic groups, or do they mix across the tree?
2.	One specimen labelled *T. antarctica* (voucher CASIZ171177) clusters with *T. dantarti* rather than *T. challengeriana*. What does this suggest about historical taxonomy?
3.	T. challengeriana spans the Weddell and Ross Seas. What oceanographic feature could explain this circumpolar distribution? (Hint: think about the Antarctic Circumpolar Current.)
4.	The interspecific K80 distance between the two species is 12–14%, while intraspecific variation is only ~1.7–1.9%. What does this barcode gap tell you about the strength of species boundaries here?
5.	Would you expect a similar pattern in the nuclear H3 marker? Why or why not, given the different effective population sizes of mitochondrial vs. nuclear loci?

<img width="1024" height="703" alt="image" src="https://github.com/user-attachments/assets/bfba04ce-ef69-421c-994a-05f8323d964e" />


| Think About It|
| --- |
| A **barcode gap** exists when the minimum interspecific distance is larger than the maximum intraspecific distance. ABGD automates its detection. |
| Is there a gap here? What biological factors could erode it (e.g. recent speciation, gene flow, incomplete lineage sorting)? |

---

## 2. Species Delimitation
  ### 2.1 Background
  Single-locus species delimitation methods look for a transition between inter- and intraspecific branching patterns. They differ in what signal they use:
| Method	| What it models	| Input needed	| Key assumption |
| --- | --- | --- | --- |
| ABGD	| Barcode gap in pairwise distances	| Aligned sequences	| Gap exists between intra- and interspecific variation |
| GMYC	| Switch from Yule (speciation) to coalescent branching rate	| Ultrametric tree (BEAST)	| Single threshold in time |
| mPTP	| Per-species coalescent model on branch lengths |	Non-ultrametric ML tree	| Each species has its own coalescent rate |
| ASAP	| Barcode gap, greedily assembled	| Pairwise distance matrix	| Simple, fast; ranks partitions by score |

| What these methods do NOT do |
| They work on a single locus, they cannot detect species whose gene trees disagree (ILS, hybridisation). |
| They cannot confirm that delimited entities are reproductively isolated. |
| Their output is a hypothesis, not ground truth. Always combine with morphology and ecology. |

  ### 2.2 Running mPTP
  `mPTP` (**multi-rate Poisson Tree Processes**) requires an unrooted ML tree with branch lengths. We will use the COI alignment.

  You have the data in the server for genes 16S and h3, but if you want to practice, you can use the script provided to fetch coi sequences from accession numbers, it is a good practice.
  
  ```bash

  conda activate lab7
  # Fetch COI sequences
  python scripts/fetch_sequences.py --email your_email@whatever.com --accessions data/accessions_coi.tsv

  # Align sequences first
  mafft --auto tritonia_COI.fasta > tritonia_coi_aligned.fasta

  # Run RAxML-NG
  raxml-ng --msa tritonia_coi_aligned.fasta --model GTR+G+I \
         --prefix coi_tritonia --seed 42 --threads 2 --all --bs-trees 100  # spolier: this gives a very common error. Try to solve it on your own with everything you learnt already
  ```

  ### 2.3 Submit the job to mPTP server
  1.	Go to: https://mptp.h-its.org/#/tree
  2.	Upload your ML tree file (e.g. coi_tritonia.raxml.bestTree).
  3.	Select the outgroup, based on the publication, probably the outgroup here corresponds to the **Proctonotidae** clade. 
  4.	Select 'Single' (PTP) or 'Multi' (mPTP) rate try both and compare.
  5.	Set your email address to receive the results and a notice when the run is finished.
  6.	Click Run and download the result file and image.

  We try in parallel using mPTP installed in the server:

  ```bash
  mptp --ml --multi --outgroup 'HM162745_Leminda-millecra_S.A.-Western-Cape-Province','GQ292060_Curnon-granulosa_Antarctica-Ross-sea','KF643932_Dirona-albolineata_unknown','DQ026831_Dirona-picta_unknown' --outgroup_crop
  ```



## 2. BAMM analysis

There are many ways to dected speciation events or diversification in evolutionary rates, for this practical we are going to start with a BAMM analysis.
BAMM (Bayesian Analysis of Macroevolutionary Mixtures) is a method used to study how diversification rates have changed across a phylogenetic tree through evolutionary time.
The central question BAMM addresses is:
-   Have all lineages in a phylogeny diversified at the same rate, or have some groups experienced periods of accelerated or slowed diversification?

In evolutionary biology, diversification refers to the balance between speciation (the formation of new species) and extinction. Many traditional analyses assume that diversification rates are constant across a phylogeny. However, this assumption is often unrealistic because different lineages may evolve under different ecological, geographical, or evolutionary conditions.

BAMM was developed to detect these rate differences directly from a time-calibrated phylogeny. Rather than assuming a single diversification rate for the entire tree, BAMM allows diversification rates to vary among lineages and through time.

For each scenario, BAMM estimates:

- The number of diversification-rate shifts on the tree.
- The location of those shifts.
- Speciation rates through time.
- Extinction rates through time.
- Net diversification rates (speciation minus extinction).

The most important concept in the BAMM is the *diversification rate shift*, which occurs when a lineage begins diversifying at a rate that differs from its ancestors or sister groups.
This is what we are going to try identify in this excercise.


### Input files

You will need:
- A rooted phylogenetic tree
- Branch lengths proportional to time (an ultrametric tree).
- Information about taxon sampling completeness.

Important: The most complete the taxon sampling is more reliable are the estimations on the shifts!
Go on the Data directory in GitHub, in the BAMM folder you will find:

- Phylogenetic tree `coi_sequences.fasta.mafft.tree`
- Control file to use in the BAMM `control_file.txt`

First we use R to prepare the tree and the control file for the BAMM. 
Download those 2 files from the the github and put them in a directory, open RStudio, and set the directory containing those files as wroking directory.

First we will transform the tree ready for BAMM, then we will calculate the Prior and the Sampling fraction. 


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

#root tree
rooted_tree <- root(
  tree,
  outgroup = c(
    "DQ026831.1_Dirona_picta_cytochrome_oxidase_subunit_1-like",
    "KF643932.1_Dirona_albolineata_voucher_10BCMOL-00344_cytochrome_oxidase_subunit_1"
  ),
  resolve.root = TRUE
)
# Plot initial tree
plot(rooted_tree, cex=0.4)
title("Original COI tree (non-ultrametric)")

# Convert to ultrametric using penalized likelihood
# (chronos = fast relaxed clock method)
chrono_tree <- chronos(rooted_tree)

# Check ultrametricity
is.ultrametric(chrono_tree)

# Compare before/after
par(mfrow=c(1,2))
plot(rooted_tree, cex=0.4)
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
```
Ok now we have the tree we need to calculate the priors and the sampling fraction. 

```R
# BAMM PIPELINE: Diversification rates in nudibranchs
# Requires ultrametric tree from chronos or BEAST

library(ape)
library(BAMMtools)
library(coda)

# Load ultrametric tree
rooted_tree <- read.tree("nudibranch_ultrametric_chronos.tre")

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
##Copy and paste the sampling_fractions.txt file in the control file in:
#globalSamplingFraction = 


# Set the BAMM priors
setBAMMpriors(tree,outfile = "bamm_priors.txt")
##You wiil need to paste the values in the outfile in a specific section of # PRIORS the control file.
```
Have a look at the control file, you have many parameters that you can modify in it, plus the name of the different output files. The ones we care the most are the `event_data.txt` and tne `mcmc_out.txt`, which will load in R later to plot our results.

Now we move on the terminal, to run the BAMM. We'll be back in R to analyse the output.
Login in the terminal, and create a directory called BAMM. Load the `control_file.txt` and the `nudibranch_ultrametric_chronos.tre` in the `BAMM` directory. Run the following command:
```sh
./bamm -c control_file.txt
```
Let it run until the end, once it's done check which outputfiles have been created. Download them and move them into the previous directory we used to generate the tree and the priors. We now go back to Rstudio.

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

You should have generated many different plot, with the tree on credible shifts.

- Which lineage started to diversify the most?
- Looking at `plotRateThroughTime(edata, rate = "speciation")` plot, at what time the shift happened?

