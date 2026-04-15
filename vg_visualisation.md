# `vg` visualisation

To interpret and visualize a complex pangenome tangle from your `odgi` graph, follow these consolidated steps. This workflow moves from initial statistical identification to a clean, readable "train track" visualization.

---

### 1. Identify the Tangle
Use `odgi stats` to find paths with high private sequence or complexity. A "tangle" is often characterized by high copy-number repeats or divergent sequences that don't align to the reference.


```bash
# show path labels and their mean links length and sum of node distances in 2D space
odgi stats -i test_results/qc/pangenome/graphs/family.retained.og -c work/42/a1b79d4c4f64628e873931db12b208/family.retained.og.lay -p -s

# how much of each path belongs to specific Pangenome Sequence Classes. These classes tell you how "shared" a sequence is across the samples in your graph.
# The columns in your output represent:
# Path Name
# Core: Nucleotides shared by all samples.
# Shell: Nucleotides shared by some (but not all) samples.
# Unique/Private: Nucleotides found only in this specific sample/path.
odgi stats -i test_results/qc/pangenome/graphs/family.retained.og -a "#,0" | sort -k2 -nr | head -10
```

### 2. Precise Subgraph Extraction
To avoid "InvalidSize" errors in rendering, extract a tiny window (100–200bp) where the variation occurs. Using a BED file is the most robust way to handle path names containing colons.

```bash
# Create a BED file manually for a 50bp micro-region
printf "HG002_PAO83395.000091:1024-5358\t1050\t1100\n" > data/micro.bed

# Extract the region with a context of 1 node
odgi extract -i family.retained.og -b micro.bed -c 1 -o micro_slice.og -O
```

### 3. Format Conversion
Convert the `odgi` binary format into a GFA bridge, then into a `vg` graph file. This ensures compatibility with the `vg` visualization engine.

```bash
# Convert OG -> GFA -> VG
odgi view -i micro_slice.og -g > micro.gfa
vg convert -g micro.gfa -p > micro.vg
```

### 4. System Font Configuration
To render the path-tracking pictographs (emoji icons) used by `vg`, ensure your environment has the necessary Unicode fonts installed.

```bash
# Install emoji and symbol fonts (Ubuntu/Debian)
apt-get update && apt-get install -y fonts-noto-color-emoji fonts-dejavu
fc-cache -fv
```

### 5. Final Visualization (SVG)
Generate the visualization using the Path (`-p`) and Node (`-n`) flags. Outputting to **SVG** is recommended to bypass pixel-dimension limits and allow web browsers to handle font rendering for pictograph icons.

```bash
# check it has content
vg paths -L -v micro.vg | wc -l
# Render to SVG with explicit font support
vg view -dpn micro.vg | dot -Tsvg -o micro_final.svg
```

---

### Key Visual Interpretations
* **Linear Segments:** Flanking sequences where all paths (samples) are in consensus.
* **Bubbles (Split/Rejoin):** Indicates a SNP or small Indel.
* **Large Loops/Cycles:** Indicates a Tandem Repeat or Copy Number Variation (CNV). If a path loops back from the end of a node to an earlier node, it represents multiple copies of that sequence.
* **Unique Colors/Icons:** Each horizontal track represents a specific haplotype. Divergence between these tracks highlights structural differences between individuals in the pangenome.
