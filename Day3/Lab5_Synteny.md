
# Synteny analysis

 **Synteny** is the conservation of gene content and gene *order* along chromosomes between species, inherited from a shared ancestor.

Synteny matters because:
- It is independent evidence for orthology (same gene by descent), on top of sequence similarity.
- Long conserved blocks tell us which ancestral chromosome arrangements have survived; broken synteny tells us where genome reshuffling happened.
- In molluscs specifically, bivalves, gastropods, and cephalopods diverged hundreds of millions of years ago and show very different amounts of genome reshuffling — making this a great system to *see* synteny erosion happen, not just read about it.

### Dataset and tutorial
Here we want to visualize the synteny of 3 mollusc species.
One per major molluscan lineage, all chromosome-level RefSeq assemblies:

| Species | Common name | Lineage | Assembly accession |
|---|---|---|---|
| *Crassostrea gigas* | Pacific oyster | Bivalvia | GCF_963853765.1 |
| *Pomacea canaliculata* | Golden apple snail | Gastropoda | GCF_003073045.1 |
| *Octopus bimaculoides* | Two-spot octopus | Cephalopoda | GCF_001194135.2 |

**IMPORTANT**: you will need chromosome-level genomes with the annotations and the proteomes, otherwise yyou want be able to infer the orthology and create the synteny blocks.

In this tutorial we are going to run DIAMOND to generate the all_vs_all output needed to run MCScanX.


```bash
cp -r /home/ubuntu/Share/mollusc_synteny/ .
cd mollusc_synteny
conda activate synteny
```
Here is the structure of the directory
```
mollusc_synteny/
├── data/        # downloaded genomes, proteomes, GFF
├── blast/       # all-vs-all DIAMOND/BLAST output
├── synteny/     # MCScanX input/output
```
In data you will find the genomes we have selected for you: E.g.  `data/Cgigas.faa`, `data/Cgigas.gff` , and for the other two species  The `.faa` should have 20,000-40,000 protein sequences (`grep -c ">" data/Cgigas.faa`).


### All-vs-all protein search with DIAMOND

We are now find similar protein pairs within and between species, this will be the evidence synteny detection is built on.
From the `mollusc_synteny` just copy paste the entire next code block in the terminal. It will take around 8 min to run. If it gets stuck for memory reason call us and we can provide the diamond output.

```bash
cd blast
declare -A CODE=( [Cgigas]=cgi [Pcanaliculata]=pca [Obimaculoides]=obi )
for sp in "${!CODE[@]}"; do
  sed "s/^>/>${CODE[$sp]}_/" ../data/${sp}.faa > ${sp}_prefixed.faa
done
cat *_prefixed.faa > all.faa

diamond makedb --in all.faa -d all
diamond blastp -q all.faa -d all -o all_vs_all.blast -e 1e-10 --max-target-seqs 5 --outfmt 6 -p 4
cd ..
```
The output file we care the most is the `all_vs_all.blast`

### Build MCScanX input and run synteny detection 

We first parse the gff files with `gff2mcscanx.py` script and the infer the synteny.

```bash
cd synteny
#this command will just create this script that will parse the gff files for MCScanX

cat > gff2mcscanx.py << 'EOF'
import sys, re
code, gff, out = sys.argv[1], sys.argv[2], sys.argv[3]
with open(gff) as fh, open(out, "w") as o:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene": continue
        m = re.search(r'GeneID:(\d+)', f[8])
        if not m: continue
        o.write(f"{code}_{f[0]}\t{code}_gene-{m.group(1)}\t{f[3]}\t{f[4]}\n")
EOF

declare -A CODE=( [Cgigas]=cgi [Pcanaliculata]=pca [Obimaculoides]=obi )
for sp in "${!CODE[@]}"; do
  python gff2mcscanx.py ${CODE[$sp]} ../data/${sp}.gff ${sp}.gff_mcscanx
done
cat *.gff_mcscanx > mollusc.gff

```
Build the protein-accession ID to gene-ID table
```sh 
cat > build_protein2gene.py << 'EOF'
import sys, re

code, gff_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
seen = set()

with open(gff_path) as fh, open(out_path, "w") as out:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = f[8]
        m_gene = re.search(r'GeneID:(\d+)', attrs)
        m_prot = re.search(r'protein_id=([^;]+)', attrs)
        if m_gene and m_prot:
            key = f"{code}_{m_prot.group(1)}"
            if key not in seen:
                seen.add(key)
                out.write(f"{key}\t{code}_gene-{m_gene.group(1)}\n")
EOF

for sp in "${!CODE[@]}"; do
  python build_protein2gene.py ${CODE[$sp]} ../data/${sp}.gff ${sp}.protein2gene
done
cat *.protein2gene > all.protein2gene
wc -l all.protein2gene   # sanity check: should be tens of thousands of lines

```
From the same directory run this:

```sh
awk 'NR==FNR{map[$1]=$2; next} ($1 in map) && ($2 in map) {print map[$1]"\t"map[$2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' \
    all.protein2gene ../blast/all_vs_all.blast > mollusc.blast

# Sanity check before running MCScanX, it shoud give you a number, if it's 0 it means the parsing failed
comm -12 <(cut -f2 mollusc.gff | sort -u) <(cut -f1 mollusc.blast | sort -u) | wc -l
```
Ready to run MCScanX
```sh
MCScanX mollusc
```

**Output:** `synteny/mollusc.collinearity` here you have the syntenic blocks. Check block count: `grep -c "^## Alignment" synteny/mollusc.collinearity`.


### Visualization

**SynVisio** 
Now follow this link and upload just upload your `mollusc.gff` and `mollusc.collinearity` files there to generate an interactive plot.
Do not get scared of the very low matching number. Go on the right side, in the filter panel, in the `source chromosome` and `target chromosome` click on select all, and then click go. These are the 3 species with an all vs all plot. 
You can try and compare Cassostrea gigas (cgi) vs Snail (pca) or Cassostrea vs Octopus (opi).

Play around Synvisio, you can zoom in specific chromosomes and check which blocks have moved.
Remember that these species have diverged more than 500 My ago, it is normal that synteny is low.
1. **Oyster vs. snail:** do you see any diagonal lines / colored ribbons indicating conserved gene order, even if patchy?
2. **Either species vs. octopus:** is synteny mostly or entirely absent? This is expected — cephalopod genomes are known to be heavily reorganized relative to other molluscs.
