# Gene Family Expansion & Contraction with CAFE5
  ## 0. Introduction
  In this practical we will explore our dataset at the phylum-level, searching for the gene families *expanded* or *contracted* along branches, exploring the rate at which gene content changes. 

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
