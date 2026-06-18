# Lab 4 - Time calibration and ancestral state reconstruction

For calibrating our phylogeny we will use IQTree3 and MCCMTREE from the Paml package. You can find more informations [here](https://iqtree.github.io/doc/Dating).
We will calibrate the tree using node dating.

Here are the calibrations we will used on the nodes:
- B: Uniform distribution
- L: Cauchy distribution

| Node number |Node | Fossil Age| Calibration |
| --- | --- | --- | --- | 
| tn_23| Root | 561 | "B(5.6100,6.0900,1e-300,0.0250)" | 
|tn_24 |Mollusca | 542 | "L(5.4200,0.1000,0.5000,1e-300)" |
|tn_25|Aculifera |497 |"L(4.9700, 0.1000, 0.5000, 1e-300)"|
|tn_26| Polyplacophora |220 |"L(2.2000, 0.1000,0.5000,1e-300)"|
|tn_41| Cephalopoda |495 |"L(4.9500,0.1000,0.5000|,1e-300)"|
|tn_37| Scaphopoda |345 |"L(3.4500,0.1000,0.5000,1e-300)"|
|tn_31| Bivalvia |525 |"L(5.2500,0.1000,0.5000,1e-300)"|
|tn_34|Origin of Caenogastropoda |405 |"L(4.0500,0.1000,0.5000,1e-300)n"|


Log into the server and go to the TimeCalibration folder.
Part of the Paml packge is now installed in IQTree3, MCMCTREE needs an Hessian matrix to calculate the approximate likelihood and infer the divergence times. 

Here is the command:

```sh

iqtree3 -s FcC_supermatrix.fas -m LG+C60+G -te Mollusca_LGC60_rooted.tree --dating mcmctree --prefix MCMCtree_mollusca

```
After running you'll have the Hessian Matrix (.hessian), and the control file to run the prior and the posterior calibration with MCMCTree (.ctl). It takes 11 hr to run, for this reason we have already generated the matrix for you. Here you have a simplified version of the control file (it's better to use this one that the one IQTree generates, as it still under development).

You will have 2 folders `Prior` and `Posterior`, each with the files needed for the estimates.

### Estimate the prior

Look at the files present in the Prior directory. 

The tree with node calibrations is  `MCMCtree_mollusca.rooted.calibration.nwk`, the alignment `MCMCtree_mollusca.dummy.phy`, and the `in.BV` (this is our hessian matrix).
Please open the `MCMCtree_mollusca.rooted.calibration.nwk` in FigTree, by checking the node label section, you can see which nodes have been calibrated.

Have a look at the mcmctree.ctl file below. This is where you specify all the details needed for MCMCTree.
Pay attention to the `usedata` flag.

- `usedata`: this variable can take three different options:

    `usedata = 0`: this option means that the sequence data present in the input sequence file (i.e., alignment block/s) will not be used as data in the MCMC. Consequently, the likelihood is set to 1, and so the target distribution to be approximated during the MCMC is the prior, not the posterior.
    `usedata = 1`: this option means that the input sequence file will be used as data during the MCMC. The likelihood will be calculated using the pruning (or "peeling") algorithm of Felsenstein (Felsenstein 1981), which is exact but very slow for large genomic datasets, although feasible to use with small datasets. This option is available for nucleotide sequence data only, and the most complex model available is HKY85+G.

    `usedata = 2` <path_inBV> and usedata = 3: these two options enable the approximate likelihood calculation (dos Reis and Yang 2011). The main workflow consists of (i) running MCMCtree using option = 3 so that it calls BASEML (if nucleotide data) or CODEML (if amino acid data) to generate the in.BV file and then (ii) run again MCMCtree with option = 2 <path_inBV> for Bayesian timetree inference under the approximate likelihood calculation

Look briefly at the full explanation of the control file [here](https://github.com/abacus-gene/paml/wiki/MCMCtree#calibrations-how-to-set-up-node-age-constraints). 

Look at the clock options and at the seq type.

Copy the section below in a new file  with nano `nano mcmctree.ctl`. Save it and now you are ready to infer the prior.

```sh
seed = -1
seqfile =  MCMCtree_mollusca.dummy.phy * alignment
treefile =  MCMCtree_mollusca.rooted.calibration.nwk * tree with fossil calibration
mcmcfile = mcmc.txt * mcmc chain output
outfile = out.txt * out file
ndata = 1 
seqtype = 2 * 0 : nucleotides; 1: codons; 2: AAs 
usedata = 0 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
clock = 1 * 1: global clock; 2: independent; and 3: correlated rates
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
BDparas = 1 1 0.1 m/c  * birth, death, sampling
rgene_gamma = 2 20 1 * gammaDir prior for rate for genes
sigma2_gamma = 1 10 1 * gammaDir prior for sigma^2 (for clock=2 or 3)
finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, mixing...
print = 1 * 0: no mcmc sample; 1: everything except branch 2: ev...
burnin = 100000
sampfreq = 100
nsample = 20000

```

Type:
```sh
mcmctree mcmctree.ctl
```
It should take a couple of minutes, once it's done download the Figtree.tree file you have generted. Load the tree in Figtree, check the node labels and node bar option to see the prior estimates.

### Posterior distribution

Go to the Posterior folder and look at the file inside. Now MCMCTree will have to take into consideration the Hessian matrix, we will specify `usedata = 2` for the approximate likelihood (this is the step that save us a lot of time with big genomics datasets.)
This time we will also specify the clock.
It is good practice to run multiple chains (e.g. 5) with both clocks, meaning that you should run 5 chains specifying clock =2 for independent rates, and 5 chains specifying clock = 3 for correlated rates. 
Due to time constrains we will only run one, but in the future, you would just need to change the clock flag and run it multiple times in different folders.
The global option (clock =1) is biologically unlikely to happen to species that are distantly related, this is because it is unlikely that all genes evolve globally at the same rate.


Ok now, once again copy the section below and save it in a new mcmctree.ctl file. We will use the independent rate clock (e.g clock = 2) and we will specify usedata = 2

```sh
seed = -1
seqfile =  FcC_supermatrix.fas * alignment
treefile =  MCMCtree_mollusca.rooted.calibration.nwk * tree with fossil calibration
mcmcfile = mcmc.txt * mcmc chain output
outfile = out.txt * out file
ndata = 1 
seqtype = 2 * 0 : nucleotides; 1: codons; 2: AAs 
usedata = 0 * 0: no data; 1:seq; 2:approximation; 3:out.BV (in.BV)
clock = 2 * 1: global clock; 2: independent; and 3: correlated rates
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
BDparas = 1 1 0.1 m/c  * birth, death, sampling
rgene_gamma = 2 20 1 * gammaDir prior for rate for genes
sigma2_gamma = 1 10 1 * gammaDir prior for sigma^2 (for clock=2 or 3)
finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, mixing...
print = 1 * 0: no mcmc sample; 1: everything except branch 2: ev...
burnin = 100000
sampfreq = 100
nsample = 20000

```


Type:
```sh
mcmctree mcmctree.ctl
```

Once this is also done, please check the tree as you did with the Prior.
Do you notice any difference in the time estimates?


### Convergence and density plots

Download on your computer the `mcmc.txt` file and open it with Tracer, by clicking on the `+` sign. You can check the different factors such as ESS, and the trace plot of your run.

We can also check the calibrations from the Prior that we used on the nodes, these plots help visualize how we are interpreting the fossil information and how we are describing the minimum and maximum bound.


```R
# Open R on RStudio and, if you have already
# installed the `mcmc3r` R package you
# can run the following commands
#
####Density plot

devtools::install_github("dosreislab/mcmc3r")
library(mcmc3r)
library(ape)

#First we can visualise the nodes on the tree
# Read tree
tree <- read.tree("MCMCtree_mollusca.rooted.calibration.nwk")

# Remove branch lengths
tree$edge.length <- NULL

# Write tree without branch lengths
write.tree(tree, file = "Mollusca_NOBL.tre")

#Root tree
tr <- read.tree("Mollusca_NOBL.tre")
tr <- root(tr, outgroup = c("Urechis_unicinctus", "Lingula_anatina"), resolve.root = TRUE)

tr$edge.length <- NULL

write.tree(tr, file = "Mollusca_NOBL_rooted.tre")
#plot
plot(tr)
nodelabels()

#Now we can see the density plots
# Create 1 row with 3 columns
# NOTE: The x axis will be reversed (i.e., from
# older to younger) so that it is
# easier to interpret the calibration densities
par( mfrow = c( 1, 3 ) )
##> Plot the uniform calibration for the root
curve( mcmc3r::dB( x, tL = 560, tU = 609, pL = 1e-300, pU = 0.025 ), 
       from = 400, to = 700, n = 1e4, xlim = rev( c( 500, 650 ) ),
       xlab = "Time (Mya)", ylab = "Density",
       main = "Root | 'B (  5.6000,  6.0900,  0.0000,  0.0250 )'" )
abline( v = c( 0.06, 0.08 ), col = "#56B4E9" )

##> Plot the soft-bound calibration for node 27 Mollusca
curve( mcmc3r::dL( x, tL = 542, p = 0.1, c= 1, pL = 0.025),
       from = 540, to = 600, n = 1e4, xlim = rev( c( 530, 550 ) ),
       xlab = "Time (Mya)", ylab = "Density",
       main = "Mollusca Node 27 | 'L(5.4200,0.1,0.05,1e-300)'" )
abline( v = c( 0.5, 2 ), col = "#56B4E9" )
```
### Full Paml tutorial

The Paml package and MCMCTree are very powerful tools for node dating, however they can be difficult to tackle. If you want a more in depth knowledge plese refer to [this tutorial](https://github.com/abacus-gene/paml-tutorial/tree/main/mcmctree-approxlnL-aa).



# Ancestral state estimation

Ok we have our tree, we have estimated when molluscs have diverged, but we want to dig further. We want to know if the last common ancestor of molluscs had a shell. To do so we are going to do an ancestral state estimation using the MK model. 
We have compiled a little character matrix that you can open in excel, called `shell_matrix.csv`. Have a look at it. For this excercise we have 2 characters:
- Shell presence or abscence
- Type of shell 

Look at the character states, when coding for the states remember to always start from 0 (e.g. if I have 2 states only, state A will be 0, and state B will be 1). Missing data are normally coded as '?'. 
**IMPORTANT** State 0 does not mean no data/missing data. The absence of something (e.g. the shell) is still a character. Not knowing if a species present a charcater has to be coded as missing data and not as an absance.

Ok now we are ready to infer the ancestral state of the LCA of molluscs.

Open Rstudio and follow the script below.

```R

### ASE practical

#install package
if (!require(remotes)){
  install.packages("remotes")
}
remotes::install_github("evo-palaeo/treesurgeon")
library(treesurgeon)

#set the working directory appropriately where you have the matrix and the tree
setwd(path-to-your-folder)
###load matrix
matrix <- read.csv("shell_matrix.csv", row.names = 1)
tp <- get_tip_priors(matrix)
names(tp) <- colnames(matrix)
tree <- read.tree("Mollusca_LGC60_BL_rooted.tre")

##We are going to use the MK model, look at the manual to have a better look on how it works
?fitMk

#This is just a way to visualize rates matrices and see how they change when we increase the states.
#We have two characters. Shell presence/absence (2 states) and type of shell (5 states)
?get_index_matrix
get_index_matrix( states = 2, model = "ER", ordered = T) #with ER model you can see that we can either have state 1 to 2 or state 2 to 1

#matrices for character 2 with 5 states
get_index_matrix( states = 5, model = "ASYM", ordered = T)
get_index_matrix( states = 5, model = "ARD", ordered = F)
get_index_matrix( states = 5, model = "SYM", ordered = F)
get_index_matrix( states = 5, model = "ER", ordered = F)

#Ancestral state of presence/absence shell. Here we include fitzjhon as root prior (equally prior prob for any of the states, pi = fitzjhon)
fitERch1 <- fitMk(tree, x = tp$ch1_shell_abs0_pres1, model = "ER", pi="fitzjohn" )
print(fitERch1)
ancER <-ancr(fitERch1)
plot(ancER)
#What does the plot tell you?

#Now let's do it for the type of shell
fitARD <- fitMk(tree, x = tp$ch2_shelltype_0noshell_1external_2valves_3internal_4plates, model = "ARD", pi="fitzjohn" )
print(fitARD)
ancARD <-ancr(fitARD)
plot(ancARD)

fitERch2 <- fitMk(tree, x = tp$ch2_shelltype_0noshell_1external_2valves_3internal_4plates, model = "ER", pi="fitzjohn" )
print(fitERch2)
ancERch2 <-ancr(fitERch2)
plot(ancERch2)
```

