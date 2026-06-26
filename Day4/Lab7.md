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

* What these methods do NOT do 
  * They work on a single locus, they cannot detect species whose gene trees disagree (ILS, hybridisation).
  * They cannot confirm that delimited entities are reproductively isolated.
  * Their output is a hypothesis, not ground truth. Always combine with morphology and ecology. 

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

  We try the same analysis in parallel using mPTP installed in the server:

  ```bash
  mptp --ml --multi --outgroup 'HM162745_Leminda-millecra_S.A.-Western-Cape-Province','GQ292060_Curnon-granulosa_Antarctica-Ross-sea','KF643932_Dirona-albolineata_unknown','DQ026831_Dirona-picta_unknown' --outgroup_crop
  ```

  ### 2.4 Interpret results
  -	How many species does mPTP delimit?
  -	Do these correspond to the two clades (*T. challengeriana* / *T. dantarti*) recovered by ML/BI?
  -	Are any specimens placed in unexpected species entities? Check voucher numbers.


  ### 2.5 Running ASAP
  1. Go to the [ASAP website](https://bioinfo.mnhn.fr/abi/public/asap/). It seems that the website is down, but you can visit this [site](https://itaxotools.org/download.html#hyperlinkDelimit) and download the program. It is pretty simple to install.
  2. Upload the COI alignment in FASTA format.
  3. Choose distance model: Kimura 2-parameter (K80).
  4. In the app version, mark the `MEGA CSV` option instead of `Generate all files`, otherwise it will be stuck.<img width="94" height="32" alt="image" src="https://github.com/user-attachments/assets/5656262e-6c95-4c3b-892f-1ac19b904456" />
  5. Click `Analyse` (in the website) or `RUN` (in the app). ASAP presents ranked partitions with an ASAP-score. # If it takes too long, you can leave it running and move to the next step.
  6. The top-ranked partition is the most likely number of species. Note the score and the partition.

  ### Looking at the output
  * Histogram
<img width="737" height="543" alt="image" src="https://github.com/user-attachments/assets/883b49fe-3c55-4cb1-977a-b2a6383691a5" />

  The graph shows the distribution of genetic distances between all pairs of analysed sequences. Two groups of bars separated by a *gap* (the **barcoding gap**) are clearly visible:
  · The peak on the left (Low distances, 0.00–0.03): Represents intraspecific variation (differences between individuals of the same species).
  · The peak on the right (High distances, >0.12): Represents interspecific variation (differences between different species).
  
  The empty space corresponding to the area between 0.03 and 0.11 where there are almost no bars is what is called the `barcode gap`. This clear separation is the statistical evidence that allows the authors to affirm that the specimens belong to different species and are not simply variations of the same population.


  ### 2.6 Comparison Exercise
  Fill in the table below with your results:
| Method | N species delimited | Corresponds to morpho-species? |
| --- | --- | --- |
| ABGD (from paper) | 2 | Yes |
| GMYC (from paper) | 2 | Yes |
| mPTP (your run)   | ? | ? |
| ASAP (your run)   | ? | ? |

Do all methods agree? If there is disagreement, which specimens are problematic?

If you like testing more possible options, you can try [GMYC online](https://species.h-its.org/gmyc/). It is the same server as bPTP, so you can still taste another method.

  ### 2.7 Discussion questions
  - Why might mPTP (multi-rate) and GMYC (ultrametric, single-rate) give different answers?
  - Could intraspecific population structure within *T. challengeriana* (Weddell Sea vs. Ross Sea) fool a delimitation method into over-splitting?
  - What additional data (more loci, more specimens, morphology) would increase confidence?

  ### 2.8 Further analyses
  In the real life, if you have information for more molecular markers available, you should use it! We provided sequences for two more molecular markers (16S and h3). You can try using these markers for the delimitation analyses and then compare if the results are consistent with the ones you got for coi.

  ### 2.9 (Optional) Species Delimitation: Discovery vs. Validation

  #### Two philosophically different operations

  Everything you have done so far in this section (ABGD, GMYC, mPTP, ASAP) belongs to the category of **species discovery** (also called *species finding* or *unsupervised delimitation*). These methods take your sequences and ask: *"How many independently evolving units are here, with no prior hypotheses?"* They are fast, useful for poorly-known groups, and essential when you have no taxonomy to start from.

  But there is a second operation: **species validation**. Here you start with one or more explicit species hypotheses (e.g., "I think there are 2 species" vs. "I think there is 1 species") and ask which hypothesis the data support better, using a formal statistical test. Validation is more rigorous but requires prior biological knowledge.

| | Discovery | Validation |
|---|---|---|
| Prior hypothesis needed? | No | Yes |
| Output | Number of candidate species | Statistical support for a specific model |
| Examples | ABGD, GMYC, mPTP, ASAP, SODA | BPP, BFD\*, PTP validation mode |
| Main risk | Over-splitting (every structured population = species) | Depends on quality of prior hypotheses |
| Best for | Poorly known taxa, barcoding | Hypothesis testing in well-studied groups |

  The over-splitting problem you saw in the discussion questions where mPTP or ASAP may fragment *T. challengeriana* into Weddell Sea and Ross Sea entities is a discovery-method artefact. A validation step would ask: *"Is the 2-species model significantly better than the 1-species model?"*

  #### BPP: a coalescent-based validation framework

  **BPP** (Bayesian Phylogenetics and Phylogeography, [Rannala & Yang 2013](https://academic.oup.com/genetics/article/194/1/245/6065416), [2017](https://academic.oup.com/sysbio/article/66/5/823/2805857)) models the multispecies coalescent explicitly. You provide:
1. A guide tree (the species topology you want to test)
2. Sequence alignments for multiple loci
3. Prior distributions on population sizes (θ) and divergence times (τ)

  `BPP` then estimates the posterior probability that each node in the guide tree represents a real speciation event. A node with posterior probability > 0.95 is considered a supported species boundary.

For *Tritonia*, you could test:
- **H1**: 2 species (*T. challengeriana* + *T. dantarti*) guided by the COI tree
- **H2**: 3 species splitting *T. challengeriana* into Weddell and Ross Sea lineages
- **H3**: 1 species everything is one population

  #### Running BPP (conceptual walkthrough)

  BPP requires multi-locus data. Since you have three markers (COI, 16S, H3), this is feasible.

  > Install [BPP](https://github.com/bpp/bpp) . We didn't install `bpp` in the server, and it can take a while, so we recommend doing it when you have the time. Now, just have a look at the commands.


  > You need:
  * 1. A control file (Imap file mapping individuals to species)
  * 2. Sequence files for each locus
  * 3. The guide tree in Newick format

  A minimal control file looks like this:
  ```bash
  seed = -1

  seqfile = tritonia_3loci.txt     * sequence data (concatenated loci in BPP format)
  Imapfile = tritonia_imap.txt     * maps each specimen to a species
  outfile = bpp_out.txt
  mcmcfile = bpp_mcmc.txt

  speciesdelimitation = 1 1 2     * rjMCMC species delimitation
  speciestree = 0                  * fixed guide tree
  species&tree = 2 T_challengeriana T_dantarti
  30 12             * number of sequences per species
  (T_challengeriana, T_dantarti);

  thetaprior = 3 0.004 e          * inverse-gamma prior on theta
  tauprior = 3 0.002              * inverse-gamma prior on tau
  ```

  Then you can run bpp like this:
  
  ```bash
   bpp --cfile bpp_control.txt
  ```

  The output gives you the posterior probability for each delimitation model. If P(2 species) > 0.95, the 2-species hypothesis is validated.

  Here there is a [tutorial](https://inria.hal.science/PGE/hal-02536475) where you can go deep in this isue, if interested.

  #### Validation matters
  The discovery methods already told you there are 2 species. But suppose you found a third cluster: some Weddell Sea specimens that appear slightly differentiated in one marker. Should you name them a third species? A discovery method might say yes. A BPP validation run would ask whether that third lineage has a significantly distinct coalescent history across *all three markers*. If not, if the signal comes from one locus only, the 2-species model will have higher posterior probability and you would not split.

  This is the key message: **discovery generates hypotheses; validation tests them**. In a robust taxonomic study, you need both.

  #### Further reading
  - Rannala & Yang (2013) *Systematic Biology* 62:523 — BPP method paper
  - Leaché & Fujita (2010) *Proceedings Royal Society B* — classic BFD* application
  - Hausdorf (2025) *Molecular Ecology* — empirical comparison showing over-splitting in discovery methods

  ### 2.10 (Optional - informative) Species Delimitation at Scale: From Single Loci to UCEs

  #### Why single-locus methods don't scale

  Everything in Sections 2.2–2.5 (mPTP, ASAP, GMYC, ABGD) was designed for **one locus**. They look for a `barcode gap` or a `coalescent-to-speciation rate transition` in a single gene tree. When you have hundreds or thousands of loci — as in a UCE, RADseq, or exon-capture dataset you face a different situation:

  - Running mPTP on each locus separately gives you hundreds of conflicting answers.
  - Concatenating all loci and running mPTP on the concatenated tree ignores gene-tree discordance and treats the dataset as if it were a single marker.
  - Neither is appropriate.

  What you actually want is a method that **models gene-tree variation across loci** under the multispecies coalescent, the same framework used by ASTRAL for species tree estimation, which you will use further in this tutorial.

  #### SODA: species delimitation from gene trees

  **SODA** (Species delimitatiOn using Discordance Analysis; [Rabiee & Mirarab 2020, *Bioinformatics* 36:5623](https://academic.oup.com/bioinformatics/article/36/24/5623/6130817)) takes a set of gene trees as input and asks: at which nodes in the species tree does the branching pattern switch from within-species coalescence to between-species divergence? It does this by testing whether internal branches in the species tree have zero length under the MSC (a zero-length branch means the lineages diverged instantaneously and are not independently evolving).

  Input:  many gene trees (one per UCE locus)

  ↓

  ASTRAL-style quartet analysis

  ↓

  Output: which nodes are "real" species boundaries vs. population structure

  #### The UCE pipeline to SODA

  Assuming you have already processed UCE data through phyluce (assembly → locus extraction → alignment), the steps are:

  **Step 1 — Build a gene tree per UCE locus**

  ```bash
  # After phyluce alignment step, you have one FASTA per locus in:
  # mafft-nexus-edge-trimmed/

  # Build a gene tree for each locus with IQ-TREE
  mkdir -p gene_trees/

  for f in mafft-nexus-edge-trimmed/*.nexus; do
      base=$(basename $f .nexus)
      iqtree -s $f -m GTR+G \
              --prefix gene_trees/$base \
              -T 2 --redo -q   # -q = quiet mode
  done

  # Collect all best trees into one file
  cat gene_trees/*.treefile > all_UCE_gene_trees.nwk
  echo "Total gene trees: $(wc -l < all_UCE_gene_trees.nwk)"
  ```

  **Step 2 — (Recommended) Phase UCE alleles first**

  UCE loci are diploid (each individual has two alleles). Phasing separates them before tree-building, which significantly improves delimitation accuracy ([Andermann et al. 2019, *Systematic Biology*](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.faircloth-lab.org/assets/pdf/andermann-et-al-2019-systematic-biology.pdf)). `phyluce` implements this in its [Tutorial II](https://phyluce.readthedocs.io/en/latest/tutorials/) pipeline:

  ```bash
  # After edge-trimming, explode alignments by taxon
  phyluce_align_explode_alignments \
      --alignments mafft-nexus-edge-trimmed \
      --input-format nexus \
      --output mafft-nexus-edge-trimmed-exploded \
      --by-taxon

  # Then align raw reads back to individual-specific contigs
  # and call + phase SNPs (see phyluce Tutorial II for full steps)
  # Output: phased alignments where each individual has two sequences
  #         (taxon_0 and taxon_1 for each allele)
  ```

  If phasing is not done, SODA still works, but may over-split because heterozygosity within individuals is mistaken for population structure.

  **Step 3 — Run SODA**

  ```bash
  # Install SODA
  [SODA](https://github.com/maryamrabiee/SODA)
  
  # Run SODA with all gene trees
  # No guide tree required (unsupervised / discovery mode)
  soda -t all_UCE_gene_trees.nwk \
       -o soda_output/ \
       --threshold 0.05     # p-value threshold for branch length test
  ```

  SODA outputs a table listing each node and whether it is supported as a species boundary (p < threshold). The final species partition is written as a mapping of individuals to candidate species.

  **Step 4 — Compare to single-locus results**

  | Method | Data | N loci | Assumptions |
  |---|---|---|---|
  | ASAP | COI (1 locus) | 1 | Barcode gap exists |
  | mPTP | COI ML tree | 1 | Per-species coalescent rate |
  | SODA | UCE gene trees | 100s–1000s | MSC; no guide tree needed |
  | BPP | UCE alignments | 10s–100s | MSC; guide tree required |

 
  SODA and other MSC methods have a well-documented tendency to over-split. [Hausdorf (2025, *Molecular Ecology*)](https://pubmed.ncbi.nlm.nih.gov/40026292/) showed that SODA delimited 233% of the expected number of species across four studied species complexes, essentially treating every genetically structured population as a distinct species. This is not a bug, it is a consequence of the MSC assumption that any detectable divergence = independent evolution. The lesson: **use SODA output as a primary species hypothesis, then validate with BPP or geographic/morphological data**.

  #### A real example: UCE delimitation in velvet worms

  [Lord et al. (2026, *Invertebrate Systematics*)](https://connectsci.au/is/article-abstract/40/1/IS25038/266190/Genomic-species-delimitation-reveals-sympatry-in?redirectedFrom=fulltext) used a UCE probe set designed for Onychophora with exactly this pipeline: `phyluce → IQ-TREE gene trees → SODA`, and discovered that two morphologically similar species of the egg-laying velvet worm *Ooperipatellus* were in fact sympatric (living in the same place), something impossible to detect from morphology alone. This is a textbook case of how big-data delimitation reveals cryptic diversity invisible to traditional taxonomy.

  
---

## 3. Haplotype Networks
A rooted phylogenetic tree forces every bifurcation to be strictly hierarchical and every ancestral haplotype to be inferred (and invisible). A haplotype network:
•	Retains observed ancestral haplotypes as nodes.
•	Shows multiple connections to represent ambiguity or recombination.
•	Visualises within-species (population-level) relationships where ILS is common and a tree-like structure may not exist.
•	Maps geography, colour morph, or other metadata directly onto nodes.

A tree is appropriate when lineages have fully sorted (diverged long ago, no gene flow). A network is appropriate within or near the species boundary, where haplotypes may be shared across populations, ancestral haplotypes may still be sampled, and incomplete lineage sorting prevents a fully resolved tree.

  ### 3.1 Building a TCS/Median-Joining Network in PopArt
  **PopArt** (Population Analysis with Reticulate Trees) is pre-installed on your computer. If not, download it from [PopArt](https://popart.maths.otago.ac.nz/download/) if needed.

  #### Preparing the input file
  PopArt reads **NEXUS** format with a TRAITS block for metadata. Convert your COI FASTA to NEXUS in R:

  ```r
  library(ape)
  seqs <- read.FASTA('data/tritonia_coi_aligned.fasta')
  write.nexus.data(seqs, file='results/COI_tritonia.nex', format='DNA')
  ```
  Check that the format is correct. It should look like this:
  <img width="607" height="178" alt="image" src="https://github.com/user-attachments/assets/b9c07f5b-a636-4fb0-af2f-dc4dd1a1c39e" />

  Then manually you could add a TRAITS block with information from geolocalisation (GPS) or another trait you want to map to the network (see PopArt [documentation](https://popart.maths.otago.ac.nz/documentation/examplenex/) for syntax).

  #### Running the network
  1.	Open PopArt > File > Open > select your .nex file.

  A warning message with appear probably with information on your sequences:
  <img width="442" height="122" alt="image" src="https://github.com/user-attachments/assets/930c14d1-a55b-4cf1-ab8c-ae38ddd4a6b9" />
  Don't worry, you can click `accept`.
  A second message will pop-up:
  <img width="469" height="132" alt="image" src="https://github.com/user-attachments/assets/30ddee87-ab11-4364-a4b3-85328c126ab1" />
  Here you decide what to do, depending on your dataset.
  
  3.	Network > TCS Network  (or try Median-Joining Network for comparison).
  4.	Use Graph > Colour nodes by trait to map location (Weddell Sea / Bouvet Island) and colour morph (orange / white).
  5.	Export figure as PDF or SVG.

  #### Interpreting the network
  - Are the two species (*T. challengeriana* and *T. dantarti*) separated by many mutational steps or few?
  - Are orange and white colour morphs interleaved within each species' clade, or do they form separate sub-clusters?
  -	Are there shared haplotypes between the two species? If yes, what could explain this?
  -	Are Weddell Sea and Ross Sea specimens of *T. challengeriana* intermixed, or is there geographic structure?

  ### 3.2 Networks in SplitsTree
  #### Phylogenetic Networks and Reticulation
  `SplitsTree` implements neighbour-net and splits-graph algorithms that display incompatible phylogenetic signal as parallelograms (boxes) rather than forcing a single bifurcating topology. Large boxes indicate conflicting signal from:
•	Recombination (reshuffling of genetic material)
•	Hybridisation / introgression (gene flow between lineages)
•	Incomplete lineage sorting (ILS) — ancestral polymorphism retained across speciation events
•	Saturation / homoplasy in the sequences

  If you didn't do it yet, download [SplitsTree6](https://software-ab.cs.uni-tuebingen.de/download/splitstree6/welcome.html) and install.

1.	Open SplitsTree > File > Open > select your aligned COI FASTA or NEXUS.
2.	The default view uses Neighbour-Net. Examine the network.
3.	Go to Analysis > Bootstrap (100 replicates) to assess support for splits.
4.	Try the same with the 16S and H3 markers separately. Do they show the same topology?

  ### 3.3 Trees vs network comparison
  | Feature	| Tree	| Network (SplitsTree) |
  | --- | --- | --- |
  | Representation	| Strictly bifurcating	| Reticulate; boxes for conflicting signal |
  | Ancestral nodes	| Inferred, not observed	| Can be observed haplotypes |
  | Conflicting signal	| Hidden / collapsed	| Explicitly shown as parallelograms |
  | Best used for	| Well-diverged lineages, tree-like evolution	| Recent divergence, ILS, introgression, recombination |
  | Within-species use	| Poor — forces bifurcation	| Excellent — shows population structure |

  
- Does the SplitsTree network for COI show large boxes (conflicting signal) or is it mostly tree-like?
- Compare the COI network to the 16S network. Are the splits consistent? Discordance between markers can indicate ILS or past gene flow.
- If you see boxes between T. challengeriana and T. dantarti, does this mean they hybridise? What alternative explanations exist?
- In which scenarios would you publish a network rather than (or alongside) a tree?

---

## 4. Hybridisation, Introgression & Recombination (Optional)
### 4.1 Biological Background
In the Science 2025 mollusc dataset, phylogenomic analyses revealed extensive gene-tree discordance across the 8 major mollusc classes. Much of this discordance is attributed to deep coalescence (ILS) rather than hybridisation, but the basal topology, particularly the placement of Monoplacophora, remains contentious. 

The mollusc-like topology you will use as your reference species tree is:
```bash
((((((((Scintilla_philippinensis,Verpa_penis)N8,Solemya_velum)N7,(Pictodentalium_vernedei,Siphonodentalium_dalli)N9)N6,((Aplysia_californica,Concholepas_concholepas)N11,((Tectura_virginea,Lottia_scabra)N13,Haliotis_cracherodii)N12)N10)N5,((Sepiola_atlantica,Octopus_vulgaris)N15,Nautilus_pompilius)N14)N4,Veleropilina_oligotropha)N3,(((Acanthochitona_discrepans,Callochiton_septemvalvis)N18,Deshayesiella_sirenkoi)N17,(Epimenia_babai,Neomenia_megatrapezata)N19)N16)N2,Lingula_anatina)N1;
```

Notice the major groups: Bivalvia (*Scintilla, Solemya, Verpa*), Scaphopoda (*Pictodentalium, Siphonodentalium*), Gastropoda (*Aplysia, Concholepas, Tectura, Lottia, Haliotis*), Cephalopoda (*Sepiola, Octopus, Nautilus*), Monoplacophora (*Veleropilina*), Polyplacophora (*Acanthochitona, Callochiton, Deshayesiella*), Aplacophora (*Epimenia, Neomenia*), and the outgroup *Lingula* (Brachiopoda).

### 4.3 Detecting Gene-Tree Discordance
We will use a multi-gene approach. You can use the single-copy gene alignments from previous sessions.

```bash
mkdir -p lab7/data
cd lab7/data/
cp -r /home/ubuntu/Share/Single_Copy_Orthologue_Sequences .  #if you already aligned and inferred the trees, copy/paste them here as well
```

**Step 1 — Build individual gene trees (if already done, go to the next step)**
* Build ML trees for each gene if you don't have them yet
```bash
mkdir -p lab7/data/alignments
mkdir -p lab7/data/trees

cd lab7/data/
for f in *.fasta; do
  mafft --auto $f > alignments/"${f%.fasta}_aln.fasta"
done

for f in aligned/*.fasta; do
  iqtree -s $f -m TESTMERGE -T AUTO
done

# then move all the tree files to the `lab7/data/trees` folder
```

* Collect best trees into one file
```bash
cat trees/*.treefile > all_gene_trees.nwk
```

**Step 2 — Visualise discordance with DensiTree**
* Open all_gene_trees.nwk in `DensiTree` # (included in BEAST package, or standalone download from [here](https://www.cs.auckland.ac.nz/~remco/DensiTree/download.html))
  If you are using the last DensiTree version, to run it, from the command line use `java -jar DensiTree.jar` from the directory where you saved the DensiTree jar file.

* Don't be scared, this is your result:
  <img width="1611" height="699" alt="image" src="https://github.com/user-attachments/assets/498aa2d1-3ac3-4896-9f96-8fb4ef64c33f" />

  The image shows a lot of conflict. This is normal because you have put together a lot of gene trees. It is especially visible in the middle part. This indicates that:
- many nodes are not shared between genes
- there is a lot of gene-tree discordance

  The red lines represent the consensus tree (or the "target tree"). It is the tree that DensiTree takes as a reference. The green ones are the other trees. When the green ones coincide with the red ones:

  → there is support.

  When the green ones diverge a lot:

  → there is conflict.

  You can play around with the display options a bit if you want.

---

## 5. BAMM analysis

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

