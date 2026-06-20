# Gene Family Expansion & Contraction with CAFE5
  ## 0. Introduction
  In this practical we will explore our dataset at the phylum-level, searching for the gene families *expanded* or *contracted* along branches, exploring the rate at which gene content changes. 
  *CAFE5* models gene family size as a birth-death process along a phylogeny: at every point on every branch, a gene family can gain a copy ("birth") or lose one ("death") at some rate. The single parameter *λ (lambda)*     is the rate of this process (how fast gene family sizes tend to change per unit branch length, genome-wide).
  
  Given:
  * a tree with branch lengths in comparable (ultrametric) time units, and
  * a table of how many gene copies each species has, for each gene family,
  
  *CAFE5* finds the value of λ that makes the observed gene counts most probable under the birth-death model (maximum likelihood), then uses that fitted model to:
    - Reconstruct the most likely ancestral gene count at every internal node of the tree, for every family.
    - Flag gene families whose size changes are statistically unlikely under the genome-wide λ — i.e., families evolving unusually fast on specific branches. These are your expansion/contraction candidates.

  
  ### Objectives
  By the end of this practical you'll be able to:
  * Understand what a birth-death model of gene family evolution is estimating (λ), and what assumptions it makes (e.g. why CAFE needs an ultrametric, time-calibrated tree, not just a topology).
  * Convert real-world bioinformatics output (OrthoFinder) into CAFE5's required input format, a routine but essential reformatting skill.
  * Run CAFE5 and interpret its core output: the global λ, and which gene families are significantly rapidly evolving.
  * Identify, for at least one gene family, on which branch(es) of the tree expansion or contraction is inferred to have happened.

  ### Datasets and software
  For this practical you'll need to use the following software (already installed in the server):
  * CAFE v5.1.0. (https://anaconda.org/channels/bioconda/packages/cafe/overview)
  * biopython (https://anaconda.org/channels/conda-forge/packages/biopython/overview)

  We'll be using the following datasets:
  * 'cafe_input.tsv' # this corresponds to the 'Orthogroups.GeneCount.tsv' reformatted for CAFE5 (we already prepared it for you).
  * 'cafe_tree_ultrametric.nwk' the topology from previous practicals, made ultrametric. We use ape::chronos() (penalized likelihood rate-smoothing) to produce a relative-time ultrametric tree directly from branch lengths.

  ## 1. Run CAFE5
  ### 1.1 First check your input files:
  ```bash
  ls -la data/cafe_input.tsv data/cafe_tree_ultrametric.nwk
  head -3 data/cafe_input.tsv
  cat data/cafe_tree_ultrametric.nwk
  ```
  **Checkpoint**: How many gene families are in cafe_input.tsv? How many species columns? Does the species count match the number of tips in the tree?

  ### 1.2 Running CAFE5 — single global λ
  This is the simplest CAFE5 model: one rate of gene family evolution for the entire tree, no exceptions. Don't forget to activate the environment.
  ```bash
  mkdir -p results/cafe_base
  cafe5 -i data/cafe_input.tsv -t data/cafe_tree_ultrametric.nwk -o results/cafe_base
  ```
  *This should run in well under a minute for a dataset this size. Watch the terminal CAFE5 prints its optimization progress as it searches for the best-fitting λ.*
  **Checkpoint**: Open 'results/cafe_base/Base_results.txt':
  ```bash
  cat results/cafe_base/Base_results.txt
  ```

  You should see something like:

  ```
  Model Base Result: <some log-likelihood score>
  Lambda: <a small number, e.g. 0.001-0.01 range typical>
  ```

  ## 1.3 What other files did CAFE5 produce?

  ```bash
  ls results/cafe_base/
  ```

  You should see several files; the ones that matter most for this practical:

  | File (name may vary slightly) | What it contains |
  |---|---|
  | `*_results.txt` | The headline: model score and estimated λ |
  | `*_asr.tre` | **A**ncestral **s**tate **r**econstruction — a tree (Nexus format) where every node, including internal/ancestral ones, has an inferred gene copy number attached |
  | `*_clade_results.txt` | Per-branch summary: how many gene families expanded vs. contracted on that branch |
  | `*_family_results.txt` (or similar — check `*_results.txt` and look for a wider table) | Per-gene-family table with significance values (see below) |

  ### 1.4 Find the gene families CAFE5 thinks are interesting

  Look at the per-family results table (it may be a wider section within `Base_results.txt`, or a separate file depending on your CAFE5 build — check what's in your `results/cafe_base/` directory and open the most 
  likely candidate). You're looking for two columns:

  - **Family-wide P-value**: is this family's overall pattern of size change unlikely under the fitted genome-wide λ? (Low p-value, conventionally < 0.01, = yes,  interesting.)
  - **Viterbi P-value** (only meaningful when the family-wide p-value is significant): tells you *which specific branch(es)* the unusual change happened on.

  ### 1.5 Questions
  1- How many gene families come out as significant (family-wide p < 0.01)?
  2- As a fraction of all families tested, does that number seem high, low, or about what you'd expect? (There's no universally "correct" answer here, it depends entirely on the data — but having an intuition for the         proportion matters more than memorising a threshold.)


    
