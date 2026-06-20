# Concatenation · Model Selection · Tree Inference · Topology Tests · Tree Support · Paralog removal · Visualisation with FigTree & iTOL

## 0. Introduction

  ### Objectives
  By the end of this practical you will be able to:

  •	Select an appropriate substitution model using model-selection criteria (AIC, BIC, BF).
  •	Infer a Maximum Likelihood (ML) phylogenetic tree with IQ-TREE.
  •	Use Bayesian Inference (BI) to infer a phylogenetic tree.
  •	Remove paralogs from the orthofinder output.
  •	Evaluate branch support using bootstrap resampling and assess its meaning.
  •	Perform topology tests (SH, AU) to compare alternative tree hypotheses.
  •	Build consensus trees from a set of trees.
  •	Visualise and annotate phylogenetic trees in FigTree and iTOL.
  •	Critically compare bootstrap support with Bayesian posterior probability.

  ### Software and Data
  - Software required: IQ-TREE (≥ 2.2), PhyloPyPruner, FigTree (≥ 1.4.4), a web browser for iTOL.
  
  - Data files provided: Orthogroup sequences with paralogs, Mollusca_FcC_supermatrix.fas, alt_topo.nwk (alternative topology).
    
    All commands assume you are working in the directory where your data files are located.

## 1. Matrix concatenation


## 2. Model Selection
  ### 2.1 Why does the model matter?
  Phylogenetic inference requires a probabilistic model of sequence evolution. The model describes how nucleotides (or amino acids) change over time. Choosing a model that is too simple can cause systematic 
  biases in topology and branch lengths; an overly complex model wastes degrees of freedom. Model selection balances fit and complexity.

  ### 2.2. Key information criteria
  | Criterion	| Description & preference |
  | --- | --- |
  | AIC (Akaike)	| Penalises each free parameter by 2. Tends to select slightly richer models. Use when prediction is the goal. |
  | AICc	| AIC corrected for small sample size. Preferred when n / k < 40 (n = sites, k = parameters). |
  | BIC (Bayesian)	| Stronger penalty (ln n per parameter). Tends to prefer simpler models. Recommended for most phylogenomic datasets. |

  ### 2.3 Running ModelFinder in IQ-Tree
  IQ-TREE integrates ModelFinder, which evaluates hundreds of substitution models efficiently. We will use the matrix generated in the previous practical that is already aligned and curated.

  ```bash
  # Run ModelFinder (standalone — no tree inference).
  iqtree -s Mollusca_FcC_supermatrix.fas -m MF --prefix model_test -nt AUTO
  # -m MF = ModelFinder only.
  ```

  #### Inspect the results file:
    
  ```bash
  cat model_test.iqtree | grep -A5 'Best-fit model' # The best model according to BIC is printed here.
  ```
  
  *TIP*
  
  The '.iqtree' log file contains the full model comparison table.
  Look for lines labelled 'AIC', 'AICc', and 'BIC' — they may differ! Always justify which criterion you use.
  Common best-fit models for chloroplast data: GTR+F+I+G4, TVM+F+G4, SYM+I+G4.

  #### Understanding the model notation
  | Symbol	| Meaning |
  | --- | --- |
  | GTR	| General Time Reversible — 6 substitution rate categories, most parameter-rich. |
  | HKY	| Hasegawa–Kishino–Yano — distinguishes transitions from transversions only. |
  | +F	| Empirical base frequencies estimated from the data. |
  | +I	| Proportion of invariable sites. |
  | +G4	| Gamma-distributed rate variation across sites (4 rate categories). |
  | +R4	| FreeRate model with 4 categories — more flexible than +G4. |

* In this website you'll find information about models: https://iqtree.github.io/doc/Substitution-Models
  
  ### 2.4 Questions
  1.	Which model was selected under BIC? Is it the same as AIC?
  2.	How many free parameters does the selected model have? (Hint: look at the 'df' column.)
  3.	What does the proportion of invariable sites (+I value) tell you about the alignment?

## 3. Maximum Likelihood Tree Inference
  ### 3.1  The ML principle
  Maximum Likelihood (ML) inference finds the tree topology and branch lengths that maximise the probability of observing the data given the model. IQ-TREE uses stochastic perturbations (random NNI moves) from multiple      starting trees to escape local optima.

  ### 3.2 Run ML inference with the best-fit model and ultrafast bootstrap.
  You already inferred a ML tree in the previous practical. Let's consolidate what you learnt:

  ```bash
  iqtree -s Mollusca_FcC_supermatrix.fas -m GTR+F+I+G4 -bb 1000 -nt AUTO --prefix ml_tree #replace GTR+F+I+G4 with the model you found in Section 1.
  # -nt AUTO lets IQ-TREE select the optimal number of CPU threads.
  ```
  #### Check the log-likelihood of the best tree:
  ```bash
  grep 'Log-likelihood of the tree' ml_tree.iqtree
  ```
  *Higher (less negative) = better.*

  #### Locate your output tree file:
  ```bash
  ls ml_tree.*
  ```
  *Key files: 'ml_tree.treefile' (best ML tree with support), 'ml_tree.iqtree' (full report).*

  ℹ  NOTE
  *IQ-TREE output files:*
  '.treefile'    — the best ML tree in Newick format (with support values)
  '.iqtree'      — full analysis report (model, tree, AIC/BIC)
  '.log'         — run log (iterations, likelihood values)
  '.contree'     — consensus tree (50% majority rule)

  * You still have another option, that is to run both analyses at the same time, find which model is the best for your data, and infer the ML tree. You can do this with the 'TESTMERGE' option of IQ-tree:
  
  ```bash
  iqtree -s Mollusca_FcC_supermatrix.fas -m TESTMERGE -bb 1000 -nt AUTO --prefix model_ml
  ```

  ### 3.3 Non-parametric vs ultrafast bootstrap

  | Feature	| Non-parametric (-b) vs Ultrafast (-B) |
  | --- | --- |
  | Speed	| Non-parametric: very slow (100 reps ≈ 100 full ML searches). Ultrafast: ~100× faster. |
  | Bias	| Ultrafast can slightly overestimate support. Use -bnni to correct. |
  | Threshold	| ≥ 70 (non-param) ≈ ≥ 95 (ultrafast) as a rough guide. Always report which you used. |
  | Recommendation	| Ultrafast (-B 1000) + -bnni for most studies. Non-parametric for small datasets. |

  ### 3.4 Questions
  4.	What is the log-likelihood score of your best ML tree?
  5.	How many parsimony-informative sites does the alignment contain? (Look in the .iqtree report.)
  6.	Would you expect the tree topology to change if you used a simpler model (e.g., JC)? Why?



## 4. Topology tests
  ### 4.1 Why test topologies?
  Sometimes we want to know whether a specific tree (e.g., one reflecting a traditional classification, or a tree from a different dataset) is statistically worse than our ML tree. Topology tests compare log-likelihoods     and assess whether differences are significant.

  ### 4.2 Common topology tests
  | Test	| Description |
  | --- | --- |
  | SH (Shimodaira–Hasegawa)	| Non-parametric. Tests whether the ML tree is significantly better than candidate trees. Conservative (fewer false positives). |
  | AU (Approximately Unbiased) |	More powerful than SH. Recommended for most comparisons. p > 0.05 means the tree cannot be rejected. |
  | KH (Kishino–Hasegawa)	| Parametric. Only valid when the alternative tree is NOT chosen after seeing the data. |

  ### 4.3 Running topology tests in IQ-Tree
  * Create a file with both trees to compare (your ML tree + the alternative).

  ```bash
  cat ml_tree.treefile alt_topo.nwk > trees_to_test.nwk
  ```
  Both trees must have the same taxa. The alternative topology is provided in 'alt_tree.nwk'.

  * Run IQ-TREE 2 with the -z flag to test topologies.

  ```bash
  iqtree -s Mollusca_FcC_supermatrix.fas -m GTR+F+I+G4 -z trees_to_test.nwk -n 0 -zb 10000 -zw -au --prefix topo_test # use the model that best fits (from section 1)
  # -n 0 skips ML search; -zb 10000 = 10,000 RELL replicates for bootstrap; -au runs the AU test.
  ```

  * Read the topology test table.

  ```bash
  grep -A20 'USER TREES' topo_test.iqtree
  ```

  *Look at the 'p-AU' column. p > 0.05 = cannot reject that tree.*
  
  #### IMPORTANT
  * A topology test asks: 'Can we reject this tree?' (NOT 'Is this the true tree?')
  * Failing to reject an alternative tree does NOT mean it is equally supported.
  * The KH test is invalid if you chose the alternative tree because it looked good (this is known as 'topology fishing').

  * You'll find more information regarding topology tests here: https://iqtree.github.io/doc/Advanced-Tutorial
    
  ### 4.4 Questions
  7.	What is the AU p-value for your ML tree? For the alternative tree?
  8.	Can the alternative topology be rejected at the 5% significance level?
  9.	What biological conclusion can you draw from the topology test result?


## 5. Optional: Bayesian Inference

  ### 5.1 Phylobayes

  ### 5.2 Checking convergence

  
## 6. Support methods: Bootstrap vs Posterior Probability
  ### 6.1 Conceptual comparison
  Both bootstrap support (BS) and Bayesian posterior probability (PP) measure confidence in a bipartition (clade), but they estimate fundamentally different quantities:
  
  | Feature	| Bootstrap (ML/MP) vs Posterior Probability (Bayesian)| 
  | --- | --- |
  | Definition	| BS: frequency of clade across pseudoreplicate datasets. PP: probability of clade given data and priors. |
  | Scale	| 0–100% (treat ≥ 70% as moderate support). PP: 0–1 (≥ 0.95 generally considered strong). |
  | Interpretation	| BS is a measure of repeatability. PP directly quantifies probability under the model. |
  | Common criticism	| BS can be conservative (underestimates). PP can be overconfident (inflated by model misspecification). |
  | Software	| IQ-TREE, RAxML (BS). MrBayes, BEAST, PhyloBayes (PP). |

  ### 6.3 Non-parametric bootstrap
  It is pretty much slower than the UF bootstrap:
  ```bash
  iqtree -s Mollusca_FcC_supermatrix.fas -m GTR+F+I+G4 -b 100 -nt AUTO --prefix ml_np_boot # change the model for the one selected in section 1
  ```

  ### 6.2 Building a consensus tree
  A *consensus tree* summarises a set of trees (e.g., bootstrap replicates, MCMC samples) into a single tree. IQ-TREE automatically generates a 50% majority-rule consensus tree ('.contree').

  * View the consensus tree generated by IQ-TREE.

  ```bash
  cat ml_tree.contree
  ```
  *This is the 50% majority-rule consensus from the bootstrap replicates, with support values on branches.*

  * Build a strict consensus from multiple tree files.

  ```bash
  iqtree -con -t ml_tree.treefile --prefix strict_con
  ```
  *A strict consensus retains only clades present in ALL trees.*

  * Extract bootstrap trees for inspection.

  ```bash
  iqtree -s chloroplast_aln.fasta -m GTR+F+I+G4 -B 1000 --wbt --prefix ml_tree_bt
  # --wbt writes the full set of bootstrap trees to ml_tree_bt.ufboot.
  ```

  *TIP*
  The '.contree' file already has bootstrap values embedded — use this file for visualisation.
  In a 50% majority-rule consensus: a clade appears if it was present in > 50% of bootstrap trees.
  A strict consensus is very conservative — one conflicting tree removes a clade entirely.

  ### 6.3 Threshold guidelines
  | Support value |	Interpretation guideline |
  | --- | --- |
  | Bootstrap ≥ 95% (ultrafast) / ≥ 70% (non-param)	| Strong support — clade is well-supported. |
  | Bootstrap 50–69%	| Moderate support — treat with caution. |
  | Bootstrap < 50%	| Weak / no support — clade may not be real. |
  | Posterior Probability ≥ 0.95	| High confidence under the model. |
  | Posterior Probability 0.75–0.94	| Moderate confidence. |
  | Posterior Probability < 0.75	| Low confidence — interpret carefully. |

  ### 6.4 Questions
  10.	List the three best-supported clades in your ML tree (highest bootstrap values). Do they make biological sense?
  11.	What are the risks of publishing a tree based solely on posterior probability without checking for model adequacy?
  12.	If a clade has 98% ultrafast bootstrap but only 0.72 posterior probability (from a Bayesian analysis), how would you interpret the conflict?

## 7. Paralog removal

**Goal:** Understand why paralogs are a problem for species-tree inference, and remove them from the orthogroups resulting from the OrthoFinder run using PhyloPyPruner.

An **orthogroup** from OrthoFinder is a set of sequences descended from a single gene in the last common ancestor of all the species in your dataset. But orthogroups often contain more than one sequence per species. These extra copies are **paralogs**, sequences related by a **gene duplication** event rather than a speciation event.

If you concatenate or build a tree from sequences that include undetected paralogs, you risk inferring a wrong topology: the tree will reflect *when genes duplicated*, not *when species diverged*.

> **Ortholog** → related by speciation → useful for species tree
> **Paralog** → related by duplication → can mislead species tree inference if not removed

[**PhyloPyPruner**](https://github.com/fethalen/phylopypruner) solves this by building a tree for each orthogroup, finding clusters of paralogs (sequences from the same species that are more similar to *each other* than to other species), and keeping only one representative per species.

## 7.1 The dataset

For this exercise, we'll use some pre-selected orthogroups from the OrthoFinder output: **`OG0001135`**.

Files provided:

- `OG0001135.fasta` protein sequences 

Look at the FASTA headers:

```
>Lottia_gigantea|A0ABT3E9M6
>Lottia_gigantea|A0ABT3EDP3
>Octopus_vulgaris|Q8DWU9
>Octopus_vulgaris|Q8DX40
>Octopus_vulgaris|Q8DX66
...
```

Notice the format: **`Species|AccessionID`**. This is required (PhyloPyPruner uses the part before the `|` as the **OTU** (operational taxonomic unit / species) and needs it to match exactly between the FASTA and the Newick tree file.)

**Question:** 
1. Which species have more than one sequence here? What might that mean biologically?

*(Answer: `Lottia_gigantea`, `Crassostrea_gigas`, `Octopus_vulgaris` (×3), `Acanthopleura_granulata`, and `Neopilina_galatheae` all have 2-3 sequences — these are our candidate paralogs.)*

## 7.2 Alignment and tree inference
In order to apply PhyloPyPruner, we need to align orthogroup sequences and then infer the gene trees. The program takes one or more multiple sequence 
alignments (MSAs) and corresponding trees as input. The alignments should be in `Fasta` format and trees in `Newick` format. The name should be the same for each corresponding alignment and tree file.

> For this exercise we already provide you with the correct headers, but you should keep in mind that **the format is very important** for everything to work smoothly with `PhyloPyPruner`.
> With the sequence files, you will need to generate the alignment, trim it, and infer the gene trees with ML in the same way you did before.

## 7.3. Look at the tree first

Open one of the trees in FigTree (or paste into [iTOL](https://itol.embl.de/)) and look at where the multi-copy species fall.

You should see something like this pattern:

```
            ┌── Crassostrea_gigas (copy 1)
        ┌───┤
        │   └── Crassostrea_gigas (copy 2)
   ─────┤
        │       ┌── Octopus_vulgaris (copy 1)
        │   ┌───┤
        │   │   └── Octopus_vulgaris (copy 2)
        └───┤
            └── Octopus_vulgaris (copy 3)
```

The paralogs from the same species cluster tightly together as **sister branches with short branch lengths** — much shorter than the branches separating different species. This is the species-overlap signal PhyloPyPruner looks for automatically.

---

## 4. Run PhyloPyPruner (5 min)

In your terminal:

```bash
cd path/to/practical2_data

phylopypruner --msa OG0001135_mollusca.fasta \
              --tree OG0001135_mollusca.tre \
              --output OG0001135_pruned \
              --mask pdist \
              --prune MI
```

**What these options do:**

| Option | Meaning |
|---|---|
| `--msa` | input alignment |
| `--tree` | input gene tree (matching OTU\|ID labels) |
| `--output` | name of the output folder |
| `--mask pdist` | when multiple paralogs from the same species form a clade, keep the one with the smallest pairwise distance to its neighbors (most "typical" copy) |
| `--prune MI` | pruning algorithm: **M**aximum **I**nclusion — keeps as many species as possible while removing paralogy |

Run it, then look inside the output folder:

```bash
ls OG0001135_pruned/
```

You should find a **pruned alignment** and a **pruned tree** containing exactly **one sequence per species** — your dataset is now ortholog-only and safe to use for species tree inference.

```bash
grep ">" OG0001135_pruned/*.fa     # or similar, check the output file naming
```

**Check:** count the sequences before and after. You started with 14 sequences across 8 species; you should end with one representative per species (8 sequences).

---

## 5. Discussion (2 min)

- What would happen if we'd skipped this step and just concatenated all 14 sequences as if they were 14 different "species"?
- PhyloPyPruner used `--mask pdist` to pick *which* paralog to keep. What's an alternative criterion, and when might it matter? *(e.g. `--mask longest` — keeping the longest sequence makes sense if you're worried about incomplete/fragmented assemblies rather than true paralogy.)*
- In a real OrthoFinder run across your whole proteome set, you won't run this on one orthogroup by hand — you'd point `--dir` at a folder containing **all** orthogroup alignments + trees, and PhyloPyPruner processes them in batch. Same logic, just scaled up.

---

## Recap

- [ ] Understood the difference between orthologs and paralogs and why it matters for species tree inference
- [ ] Inspected a gene tree and visually identified paralog clusters
- [ ] Ran PhyloPyPruner on one orthogroup using `--msa`, `--tree`, `--mask`, `--prune`
- [ ] Confirmed the pruned output contains one sequence per species

**Next:** we'll take pruned, single-copy orthogroups like this one — but now across hundreds of genes — and concatenate them into a supermatrix for tree inference.

## 5. Visualising trees
  ### 5.1 Loading and orienting the tree in FigTree
  * Open FigTree and load your best ML tree.
  *File → Open → ml_tree.treefile*
  When prompted 'Define name for values', type: Bootstrap

  * Root the tree on an outgroup taxon.
  Select the outgroup *taxon label → right-click → Root on branch*
  Click the branch leading to your outgroup, then use *Tree → Root on branch*.

  * Show support values on branches.
  Tick 'Branch Labels' in the left panel → *Display: Bootstrap*
  Adjust font size under 'Font Size' for readability.

  * Colour-code taxa by group.
  *Select taxa → Highlight → pick a colour → repeat for each group*
  Use *Tip Labels → Colour By metadata* if you have a tab-separated metadata file.

  * Export a publication-ready figure.
  *File → Export → PDF* (or SVG for editing in Inkscape)
  Avoid PNG for publications, vector formats (PDF, SVG) scale without loss.

  * You can also import an annotation file to change labels or even add colours to your tips.
  *File → Import Annotations → 'codes.txt'*
  *Tip Labels → Display: Mollusc names*

  ### 5.2 FigTree Panels reference
  | Panel / option	| What it controls |
  | --- | --- |
  | Trees	| Switch between multiple trees if more than one is loaded. |
  | Layout	| Rectangular, Radial, Polar — choose based on the number of taxa. |
  | Appearance → Line weight	| Branch line width — increase for clarity. |
  | Scale Bar	| Toggle a scale bar showing substitutions/site. |
  | Node Bars	| Show node age uncertainty (95% HPD) — useful for time-calibrated trees. |
  | Colouring → Gradient	| Colour branches by rate or any continuous trait. |
  
  <img width="202" height="729" alt="image" src="https://github.com/user-attachments/assets/992b1724-37e2-4783-a1ba-9af10bfc18ca" />
  <img width="196" height="194" alt="image" src="https://github.com/user-attachments/assets/3dc72921-7228-4560-87bb-75064845342e" />
  <img width="194" height="167" alt="image" src="https://github.com/user-attachments/assets/54775393-0db2-4cdf-a207-a56abbd75d8a" />
  <img width="200" height="165" alt="image" src="https://github.com/user-attachments/assets/f0864c7a-0d73-466c-8c2f-401ae2082b4b" />


  ### 5.3 Visualising trees with iTOL
  
  * Go to the iTOL website.
    https://itol.embl.de
    - No account needed for basic use; register for free to save trees.

  * Upload your tree file.
    Click *'Upload' → drag ml_tree.contree* into the upload box → click 'Upload'
    - Newick format ('.treefile', '.nwk', '.tree') is accepted.

  * Click on your tree to open the interactive viewer.
    Click the tree name in your workspace
    - The tree opens with a default rectangular layout.

  * Explore display options.
    Use the top panel: Basic → Unrooted / Circular / Normal | Advanced → Node IDs, Branch lengths

  * Add a colour dataset (metadata annotation).
    Datasets panel → Click colour strip icon → paste or upload metadata TSV
    - Format: TAXON_NAME\tCOLOUR (e.g. #FF0000)

  * Export the figure.
    Export panel → SVG (recommended) or PNG
    - SVG can be edited in Inkscape or Adobe Illustrator for final touch-ups.

<img width="1909" height="588" alt="image" src="https://github.com/user-attachments/assets/d2780754-77d6-44bb-9cf6-e76f0ebe3ba9" />
<img width="412" height="662" alt="image" src="https://github.com/user-attachments/assets/a80769b4-34d3-4d30-bcc0-666d1cfbcf11" />
<img width="396" height="232" alt="image" src="https://github.com/user-attachments/assets/742ba1e4-a3db-4cc8-bf99-9cd3d43e53a2" />



  ### 5.4 iTOL annotation file formats
  
  | File type	| Purpose and example |
  | --- | --- |
  | DATASET_COLORSTRIP	| Colour bar at tip labels. One colour per taxon. Use for taxonomic groups. |
  | DATASET_BINARY	| Presence/absence matrix. Draws filled/empty circles at tips. |
  | DATASET_SIMPLEBAR	| Bar chart at tips. Use for quantitative data (e.g. genome size). |
  | DATASET_TEXT	| Custom text label next to tips. Use for accession numbers, common names. |
  | DATASET_SYMBOL	| Symbols (circles, squares) at tips or internal nodes. |
  | TREE_COLORS	| Change branch or label colours, styles, and widths. |

  ### 5.5 Comparing two trees in iTOL
  * Upload both trees to your iTOL workspace.
  *Upload 'ml_tree.contree' and 'alt_topo.nwk' separately*

  * Open the first tree, then use 'Compare' mode.
  In the tree viewer: *Advanced → Compare → select the second tree*
  - iTOL will highlight nodes that differ between the two topologies.

  * Interpret the differences.
  - Note which clades are resolved differently and relate this back to your topology test results from Section 3.

  *TIP*
  * iTOL is ideal for large trees (> 100 taxa) where FigTree becomes slow.
  * Save your iTOL tree with a free account to get a permanent shareable URL — useful for supplementary materials in publications.
  * The 'collapse' feature in iTOL lets you fold entire clades to simplify figures.

## 6. General discussion and questions
Answer the following questions individually or in groups. Be prepared to discuss your answers.

13.	What are the consequences of choosing a substitution model that is too simple for your data? Can you think of a biological example where this matters most?

14.	Your colleague argues that a clade with 100% ultrafast bootstrap support 'must be true'. How would you respond? What factors could still make the clade an artefact?

15.	You have two trees: Tree A (ML, bootstrap = 88%) and Tree B (Bayesian, posterior = 0.99) disagree on the placement of one taxon. How do you decide which to trust? What additional analyses could help?

16.	Explain in your own words the difference between a consensus tree and a single best ML tree. In which situations would you prefer to report each?

17.	Why is it important to apply a topology test when comparing an ML tree to an alternative tree based on morphological data?

18.	A phylogenomic study uses 500 genes but recovers very short internal branches with low bootstrap values. What phylogenetic phenomenon might explain this, and how could you test for it?




    





  








