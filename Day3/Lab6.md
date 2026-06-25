# Reconciling Gene and Species Trees: DTL & Ancestral Sequence Reconstruction (ASR)
## 0. Before you start: what is this dataset?

Following previous practicals, today's dataset is a 20-species clade of molluscs (plus one outgroup), built on the real topology published in Chen et al.'s genome-based mollusc phylogeny in [*Science* (2025)](https://www.science.org/doi/10.1126/science.ads0215): clams, snails, chitons, aplacophorans, and cephalopods, all related exactly as that study reconstructed them. The protein sequences for this practical are simulated, not real genomic data, but the tree shape and the branch pattern you'll work with are taken from a real, peer-reviewed phylogeny rather than an arbitrary random tree, so the relationships you'll be reasoning about today (who is sister to whom, how deep the cephalopod radiation is, where the root sits) are biologically meaningful.

Two of the six core genes were also used as the starting point for two **constructed teaching families** [(Section 1.3)](https://github.com/MEleonoraRossi/BigDataPhylogeny-Course/blob/main/Day3/Lab6.md#13-the-engineered-dtl-teaching-families): one with an engineered duplication and loss, one with an engineered horizontal transfer. You are not discovering unknown biology in those two families, you already know, from this handout, exactly what was inserted into them. The point of working with a "known answer" gene family first is one of the most standard sanity checks in this kind of analysis: before trusting a reconciliation method on real, unknown data, you want to see it correctly recover a case where the right answer is known in advance.

## 1. Data overview

### 1.1 The species tree
Create a folder for this practical and copy the data in the folder created:
```bash
mkdir -p Lab6/data
cd Lab6/data/
cp /home/ubuntu/Share/Lab6/data/species_tree_rn.nwk .
```

`species_tree_rn.nwk` contains a fixed, rooted, bifurcating tree of 20 species, with *Lingula anatina* (a brachiopod, not a mollusc) as the outgroup, the same outgroup choice used in the source study. Internal nodes are pre-labelled N1 to N19 so you can refer to specific ancestors by name (for example, "the ancestor of the cephalopods is N14").

You can view the tree directly in the terminal at any point with:

```bash
conda activate lab6

nw_display species_tree_rn.nwk
```
<img width="1101" height="506" alt="image" src="https://github.com/user-attachments/assets/8f9de8d2-cc3f-4253-bfc7-778d4a693a65" />


or open it in any Newick-compatible tree viewer (FigTree, iTOL) for a graphical view.

A quick guide to the 20 species and where they sit on the tree:

| Clade (node) | Species |
|---|---|
| Cephalopods (N14) | *Sepiola atlantica* (bobtail squid), *Octopus vulgaris*, *Nautilus pompilius* |
| Gastropods (N10) | *Aplysia californica*, *Concholepas concholepas*, *Tectura virginea*, *Lottia scabra*, *Haliotis cracherodii* |
| Bivalves & scaphopods (N6) | *Scintilla philippinensis*, *Verpa penis* (the "watering-pot shell"), *Solemya velum*, *Pictodentalium vernedei*, *Siphonodentalium dalli* |
| Monoplacophoran (N3) | *Veleropilina oligotropha* |
| Chitons (N17) | *Acanthochitona discrepans*, *Callochiton septemvalvis*, *Deshayesiella sirenkoi* |
| Aplacophorans (N19) | *Epimenia babai*, *Neomenia megatrapezata* |
| Outgroup | *Lingula anatina* (brachiopod, not a mollusc) |

N16 groups the chitons and aplacophorans together (the Aculifera clade); N2 is the ancestor of all 19 molluscs in the tree, with *Lingula anatina* attaching at the root, N1.

### 1.2 The core gene set (single-copy orthologues)
Now copy the sequences, they are in two different folders:
```bash
cp -r /home/ubuntu/Share/Lab6/data/core_gene_set/ .
cp -r /home/ubuntu/Share/Lab6/data/dtl_engineered/ .
```

`core_gene_set/` contains six single-copy protein-coding gene families, each present in all 20 species exactly once. They were chosen because they are short, well-conserved, easy-to-align genes, exactly the kind of gene family you would use to first sanity-check a pipeline in a real research project, before tackling messier, larger gene families.

| Gene | File | Length (aa) | Notes |
|---|---|---|---|
| RS9 | `RS9_OG0000545.fasta` | 130, no indels | 30S ribosomal protein S9, used for the **ASR exercise** |
| RS7 | `RS7_OG0000433.fasta` | 156, no indels | 30S ribosomal protein S7 |
| SSRP | `SSRP_OG0000396.fasta` | ~155 | SsrA-binding protein |
| TRMD | `TRMD_OG0000263.fasta` | ~240, has small indels | tRNA (guanine-N1)-methyltransferase |
| SCPA | `SCPA_OG0000226.fasta` | ~235, has small indels | base gene used to build the engineered duplication family |
| RSMA | `RSMA_OG0000272.fasta` | ~290 | base gene used to build the engineered HGT family |

Because every gene here is a genuine single-copy orthologue (one gene, one species, no duplicates, no losses), reconciling any of them against the species tree will normally show you **zero duplications and zero transfers**. This is not a failure of the method, it is what a "clean," fully vertically-inherited gene family is supposed to look like, and it's an important reference point before looking at messier data.

### 1.3 The engineered DTL teaching families

Because the core gene set has no real duplications, transfers, or losses to find, `dtl_engineered/` contains two gene families built by deliberately editing the core genes. In short:

**FAMILY_A_DUP_LOSS** (built from the SCPA gene): a duplication was inserted on the branch leading to the cephalopod clade (*Sepiola atlantica*, *Octopus vulgaris*, *Nautilus pompilius*: node N14 on the species tree), so that *Sepiola* and *Octopus* now each carry two paralogous copies of this gene, while *Nautilus* only has one (the second copy was "lost"). All other 17 species are untouched, single-copy sequences.

**FAMILY_B_HGT** (built from the RSMA gene): *Lingula anatina*'s gene was replaced with a lightly mutated copy of *Octopus vulgaris*'s gene, simulating a horizontal transfer between two very distant lineages, *Lingula* sits at the very base of the tree as the outgroup, while *Octopus* is deep inside the cephalopod clade, on the opposite side of the tree. All other 19 species keep their original, unmodified sequence.

You already know the right answer for both families. Today's exercise is to check whether GeneRax recovers it.

---

## 2. DTL reconciliation with GeneRax 
### 2.1 What problem are we solving?

A gene tree built from a sequence alignment does not always have the same shape as the species tree. There are several reasons this can happen: a gene was duplicated at some point, so today's genomes carry two or more copies that have since evolved independently (duplication, D); a gene copy was lost in some lineages (loss, L); or a gene moved horizontally between two lineages that are not directly related by descent (horizontal transfer, T). "Reconciliation" is the process of explaining the differences between a gene tree and a species tree as a minimal or most-likely combination of speciation, duplication, loss, and transfer events.

`GeneRax` does this under a maximum-likelihood framework: given a multiple sequence alignment, a fixed species tree, and (optionally) a starting gene tree, it searches for the gene tree topology and the history of D/T/L events that together best explain both the sequence data and the reconciliation, jointly.

### 2.2 Align the sequences

Before any tree can be built, the sequences need to be aligned column-by-column so that homologous positions line up.

```bash
for f in *.fasta; do mafft --auto $f > aligned/"${f%.fasta}_aln.fasta"; done # move to the correct folder where the fasta sequences are in both 'core_gene_set' and 'dtl_engineered' datasets
```

This runs MAFFT with its automatic mode (`--auto`) on every FASTA file in `Lab6/data/core_gene_set/` and `Lab6/data/dtl_engineered/`, writing the alignments to matching `aligned/` subfolders. Open one or two of the resulting `_aln.fa` files in a text editor or alignment viewer (Jalview, AliView, or even just `less` in the terminal) and look at how similar the sequences already are.

### 2.3 Run GeneRax
For doing so, you have a script available that needs the alignments and the species tree in their correct place:

```bash
cd /home/user_number/Lab6/
mkdir -p scripts
mkdir -p data/species_tree
mkdir -p results
cp /home/ubuntu/Share/Lab6/scripts/generax.sh scripts/
cp /home/ubuntu/Share/Lab6/data/species_tree/species_tree.nwk data/species_tree/

bash scripts/generax.sh
```

This script runs three separate GeneRax jobs, all against the same fixed species tree (`Lab6/data/species_tree/species_tree.nwk`):

1. The **baseline** core gene set (RS7, SSRP, TRMD, SCPA, RSMA), under the full duplication+transfer+loss model (`UndatedDTL`).

```bash
generax \
    --families data/core_gene_set/families_baseline.txt \
    --species-tree data/species_tree/species_tree.nwk \
    --rec-model UndatedDTL \
    --per-family-rates \
    --reconcile \
    --strategy SPR \
    --prefix results/baseline_UndatedDTL
```

2. The **engineered** families, also under `UndatedDTL`.

```bash
generax \
    --families data/dtl_engineered/families_engineered.txt \
    --species-tree data/species_tree/species_tree.nwk \
    --rec-model UndatedDTL \
    --per-family-rates \
    --reconcile \
    --strategy SPR \
    --prefix results/engineered_UndatedDTL
```

3. The **engineered** families again, but under a more restrictive model that allows duplication and loss only, with no transfers (`UndatedDL`).

```bash
generax \
    --families data/dtl_engineered/families_engineered.txt \
    --species-tree data/species_tree/species_tree.nwk \
    --rec-model UndatedDL \
    --per-family-rates \
    --reconcile \
    --strategy SPR \
    --prefix results/engineered_UndatedDL
```


Each run produces a `results/<run_name>/reconciliations/` folder containing, for every gene family:

- `<family>_reconciliated.xml` — the reconciled gene tree, in recPhyloXML format (the standard exchange format for this kind of result, readable by several visualisation tools).
- `<family>_speciesEventCounts.txt` — a table of how many speciations, duplications, losses, and transfers happened on each branch of the species tree.
- `<family>_transfers.txt` — for every inferred transfer, which species was the donor and which was the receiver.

### 2.4 Read the results

Raw GeneRax output is not very readable at first glance, so a small helper script highlights the branches where something other than plain speciation happened:

```bash
python3 scripts/summarize_dtl.py results/baseline_UndatedDTL/reconciliations/RS7_speciesEventCounts.txt
python3 scripts/summarize_dtl.py results/engineered_UndatedDTL/reconciliations/FAMILY_A_DUP_LOSS_speciesEventCounts.txt
python3 scripts/summarize_dtl.py results/engineered_UndatedDTL/reconciliations/FAMILY_B_HGT_speciesEventCounts.txt
```

**Questions to answer (write down your reasoning, not just the numbers):**
- For the baseline gene set (RS7, SSRP, TRMD, SCPA, RSMA), how many duplication, loss, and transfer events did GeneRax infer in total, across all five families? Is this what you expected, given what you read in Section 1.2? Why?
- For `FAMILY_A_DUP_LOSS`, which branch of the species tree shows a duplication event? Does it match where we said we inserted one (the branch leading to N14, the cephalopod ancestor)? Which branch shows the loss?
- For `FAMILY_B_HGT` under the `UndatedDTL` model, is there a transfer event, and if so, between which two species/branches? Does it match the *Octopus vulgaris* → *Lingula anatina* transfer we engineered?

### 2.5 Comparing evolutionary scenarios: DTL vs DL-only

`FAMILY_B_HGT` was reconciled twice: once allowing transfers (`results/engineered_UndatedDTL`) and once forbidding them (`results/engineered_UndatedDL`). Compare the two:

```bash
python3 scripts/summarize_dtl.py results/engineered_UndatedDTL/reconciliations/FAMILY_B_HGT_speciesEventCounts.txt
python3 scripts/summarize_dtl.py results/engineered_UndatedDL/reconciliations/FAMILY_B_HGT_speciesEventCounts.txt
```

**Questions:**

- When transfers are not allowed, how does GeneRax explain the anomalous position of *Lingula anatina* in this gene family instead? Look at the duplication and loss counts: did they go up?
- Which of the two scenarios (DTL vs DL-only) is the more biologically parsimonious or plausible explanation for what we know was engineered into this family? Which one do you think a researcher analysing *real*, unknown data should prefer, and on what grounds (model fit / likelihood, not just "fewer events")?
- This is the general lesson behind "comparing evolutionary scenarios": the same gene tree and the same species tree can be reconciled under different modelling assumptions, and those assumptions can change the entire biological story you would tell about a gene family. Can you think of a real biological situation (not necessarily molluscs) where assuming "no horizontal transfer is possible" would lead a researcher to a wrong conclusion about gene duplication history?

### 2.6 (Optional, if time allows) Look at the reconciled tree picture

GeneRax's `.xml` output can be turned into a picture showing the gene tree embedded inside the species tree, with little symbols marking speciations (circles), duplications (squares), losses (crosses), and transfers (arrows). You can install `Notung`, that is a software java-dependent where you can reconcile your gene trees against the species trees:

**[Notung](https://amberjack.compbio.cs.cmu.edu/Notung/download.html#)**

<img width="1908" height="1010" alt="image" src="https://github.com/user-attachments/assets/8153069a-25c4-4ec3-a1a9-1faecf494601" />



You can also try [icytree](https://icytree.org/) to visualise the xml trees.

[Thirdkind](https://github.com/simonpenel/thirdkind) is another option for visualisation. It converts the reconciled trees into svg images.


There are many different options where you can apply reconciliation. One of them is using [SpeciesRax](https://academic.oup.com/mbe/article/39/2/msab365/6503503) to infer phylogenies. Because of time constraints, we are not using it here.


---

## 3. Ancestral Sequence Reconstruction with RAxML-NG 

### 3.1 What problem are we solving?

Given an alignment of present-day sequences and a tree relating them, can we infer what the sequence looked like at an ancestral node, for example, at the common ancestor of all 20 species in our tree (node N1), or at the ancestor of just the cephalopods (node N14)? This is **Ancestral Sequence Reconstruction (ASR)**. RAxML-NG's `--ancestral` mode computes, for every internal node of a fixed tree, the marginal probability of each possible amino acid at every alignment column, then reports both the full probability table and a single best-guess sequence (the amino acid with the highest probability at each site).

### 3.2 Why RS9?

We picked the RS9 gene (30S ribosomal protein S9) for this exercise because all 20 species have exactly the same sequence length (130 amino acids) with no insertions or deletions anywhere in the alignment. That means every column of the alignment is unambiguously homologous across all 20 species, and every one of the 130 reconstructed ancestral states corresponds to one single, unambiguous position. This lets you focus entirely on what ASR means and how to read its output, without first having to reason about gaps and indel placement, which is its own, harder topic.

### 3.3 Run the reconstruction

```bash
bash scripts/asr.sh
```

This runs:

```bash
raxml-ng --ancestral \
    --msa data/asr/RS9_OG0000545.fasta \
    --tree data/asr/species_tree.nwk \
    --model LG+G \
    --prefix results/ASR_RS9
```

Three output files matter:

- `results/ASR_RS9.raxml.ancestralStates` — one reconstructed sequence per internal node, with the node names matching the labels in the species tree (N1, N2, ... N19).
- `results/ASR_RS9.raxml.ancestralProbs` — for every site and every internal node, the full posterior probability of each of the 20 amino acids (so you can see *how confident* the reconstruction is, not just the single best guess).
- `results/ASR_RS9.raxml.ancestralTree` — the tree with the node labels used in the two files above (should match `species_tree.nwk` exactly).

### 3.4 Reading the output

Open `results/ASR_RS9.raxml.ancestralStates` in a text editor. Each line has a node name followed by a 130-character amino-acid sequence: the single most likely ancestral sequence for that node.
Watch out because raxml-ng changed the node labels, and now the root node and all the others are renamed. If you want to find out what is the new name in order to make comparisons, you can run this short script in R:

```r
library(ape)

tree <- read.tree("ASR_RS9.raxml.ancestralTree")
plot(tree)
nodelabels()
tiplabels()

#Find root node label
root_node <- getMRCA(tree, tree$tip.label)
root_node

# Find cephalopod node
tree$tip.label  # first make sure that you know how are the tips written
ceph_node <- getMRCA(tree, c("Octopus_vulgaris", "Sepiola_atlantica", "Nautilus_pompilius"))  
ceph_node
```

**Questions:**

1. Compare the reconstructed sequence at N1 (the root, the common ancestor of all 20 species, including the outgroup *Lingula anatina*) with the sequence at N14 (the ancestor of the three cephalopods). How many positions differ? Is this what you would expect, given how long ago these two ancestors lived relative to each other on the tree?
2. Pick any one of the 20 present-day species and compare its real sequence (from the alignment) to the reconstructed sequence at its immediate ancestral node (its parent in the species tree). How many differences are there? What does a small number of differences tell you about how conserved this gene is along that particular branch?
3. Open `results/ASR_RS9.raxml.ancestralProbs` and find one site (column) where the probability is split fairly evenly between two amino acids (for example, no single amino acid has more than 60% posterior probability), versus a site where one amino acid has more than 99% probability. What does the difference between these two situations tell you about how much we should trust the single "best guess" sequence reported in `.ancestralStates`, site by site?

### 3.5 ASR on a gene tree instead of the species tree

Everything above used the fixed species tree as the backbone for ASR, because we know this gene has no duplications or losses (Section 1.2), so its gene tree and the species tree should have essentially the same shape anyway. In a real research project working with a larger, messier gene family, you would normally use the gene family's own ML gene tree (for example, the one GeneRax inferred for you in [section 2](https://github.com/MEleonoraRossi/BigDataPhylogeny-Course/blob/main/Day3/Lab6.md#23-run-generax)) for the ASR step instead, since that is the tree that actually describes the gene's history, not the species' history. If you have time, try rerunning Section 3.3 using one of the gene trees produced by GeneRax in section 2 (look inside `results/baseline_UndatedDTL/` for an inferred gene tree file) instead of `data/asr/species_tree.nwk`, and see whether your reconstructed ancestral sequences change.

---

## 4. Wrap-up discussion

Be ready to discuss as a class:

- What is the practical difference between "this gene family shows zero D/T/L events" and "this method correctly detected a known duplication and a known transfer"? Why did today's lab need both kinds of gene families to teach the topic properly?
- ASR gave you one most-likely sequence per ancestral node, but also a full probability distribution at every site. When would a researcher care about the second part more than the first?
- Connecting the two halves of the lab: if a gene family had real, undetected duplications, why would that be a problem for an ASR analysis that (incorrectly) assumed the gene tree was identical to the species tree?

