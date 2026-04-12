# ONT Long Reads → CHM13 Pangenome Alignment (chr22) — E2E Pipeline

Aligning to a graph reduces "reference bias" by accounting for structural variations in the population.
**Chemistry:** ONT R10.4.1 (PromethION, Kit 14, HAC basecalling)
**Reference pangenome:** HPRC Freeze 1 / minigraph-cactus / CHM13 backbone
**Mapper:** vg giraffe ≥ 1.73.0 (long-read mode)
**Visualization:** odgi viz

---

## Part 1 — Prepare the Graph

```bash
# Download Chr22 graph (.vg)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.chroms/chr22.vg

# Convert .vg → GFA
vg view chr22.vg > chr22.gfa

# Fix path name: add Haplotype index (#0) required by downstream tools
sed 's/CHM13#chr22/CHM13#0#chr22/g' chr22.gfa > chr22.fixed.gfa

# Create GBZ (tags CHM13 as the reference coordinate system)
vg gbwt -G chr22.fixed.gfa -g chr22.gbz --set-reference CHM13

# Convert GBZ → mutable .vg format (required for vg augment later)
vg convert chr22.gbz > chr22_mutable.vg

# Verify CHM13 path is present with correct name
vg paths -x chr22_mutable.vg -L | grep CHM13
# Expected: CHM13#0#chr22
```

---

## Part 2 — Build Long-Read Indexes

> **Critical:** HPRC pre-built indexes predate long-read support in giraffe.
> The existing `.min` and `.zip` files must be regenerated with `--workflow lr-giraffe`.
> Do NOT reuse short-read indexes for ONT data.

```bash
vg autoindex \
    --workflow lr-giraffe \
    --prefix chr22_lr \
    --gbz chr22.gbz

# Produces:
#   chr22_lr.dist
#   chr22_lr.longread.withzip.min
#   chr22_lr.longread.zipcodes
```

---

## Part 3 — Prepare FASTQ Reads

Zero-pad read names so lexicographic sort matches numeric sort in visualizations.
Adjust the digit count to match your read count (e.g. 4 digits for ≤9999 reads, 6 for ≤999999).

```bash
# Example: rename reads from a tiny test subset (1000 reads → 4-digit padding)
zcat input.fastq.gz | \
    awk 'NR%4==1 {printf "@SAMPLENAME.%04d\n", int((NR+3)/4)} NR%4!=1 {print}' | \
    gzip > sample_renamed.fastq.gz

# Verify ordering
zcat sample_renamed.fastq.gz | grep "^@" | head -5
# Expected: @SAMPLENAME.0001, @SAMPLENAME.0002 ...

# For full-scale runs (100k+ reads), use 6-digit padding:
# awk 'NR%4==1 {printf "@HG002.%06d\n", int((NR+3)/4)} NR%4!=1 {print}'
```

---

## Part 4 — Map Reads with Giraffe (Long-Read Mode)

Use `-b r10` for ONT R10.4.1. Use `-b hifi` for PacBio HiFi.

```bash
vg giraffe \
    -b r10 \
    -Z chr22.gbz \
    -d chr22_lr.dist \
    -m chr22_lr.longread.withzip.min \
    -z chr22_lr.longread.zipcodes \
    -f sample_renamed.fastq.gz \
    -p \
    --threads 8 \
    > sample.gam

# Check mapping quality
vg stats -a sample.gam
# Healthy indicators:
#   Total aligned: >95% of reads
#   Median mapping quality: 60
#   Softclips: <25% for chr22-subsetted reads
```

---

## Part 5 — Filter Unmapped Reads

Removes reads with MAPQ=0 (no seeds found, or ambiguous) before augmentation.
Unmapped reads produce empty paths in the graph which cause odgi build to fail.

```bash
vg filter \
    -q 1 \
    sample.gam \
    > sample_mapped.gam

# Confirm count (should be >95% of original)
vg stats -a sample_mapped.gam | grep "Total aligned"
```

---

## Part 6 — Augment Graph with Read Paths

> **Note:** `vg augment` requires a mutable graph format.
> GBZ is read-only — always use the `chr22_mutable.vg` created in Part 1.

```bash
vg augment \
    -t 8 \
    chr22_mutable.vg \
    sample_mapped.gam \
    --label-paths \
    > sample_augmented.vg

# Verify sample read paths were added
vg paths -x sample_augmented.vg -L | grep "SAMPLENAME" | wc -l
```

---

## Part 7 — Build odgi Graph

```bash
# Convert to GFA forcing P-lines (not W-lines) for odgi compatibility
vg convert sample_augmented.vg -f -W > sample_augmented.gfa

# Safety check: confirm no empty paths (would cause odgi build to fail)
awk '$1=="P" && $3==""' sample_augmented.gfa | wc -l
# Must return 0

# Verify CHM13 reference path survived conversion
grep "^P" sample_augmented.gfa | grep "CHM13"

# Build odgi graph
odgi build -g sample_augmented.gfa -o sample_augmented.og -t 8 -P

# Confirm path counts (HPRC haplotypes + sample reads)
odgi paths -i sample_augmented.og -L | wc -l
odgi paths -i sample_augmented.og -L | grep "SAMPLENAME" | wc -l
```

---

## Part 8 — Prepare Display Path List

Used to restrict visualizations to CHM13 reference + your sample reads only.

```bash
odgi paths -i sample_augmented.og -L \
    | grep -E "^CHM13|^SAMPLENAME" \
    | tr -d '\r' \
    | awk '!seen[$0]++' \
    > display_paths.txt

wc -l display_paths.txt
```

---

## Part 9 — Visualize

### 9a. Compressed Coverage Heatmap (recommended first view)

Collapses all read paths into a single heatmap row alongside reference haplotypes.
Shows where reads covered the pangenome. Best for a quick sanity check.

> **Note:** `-O` (compressed mode) cannot be combined with `-p` (path filter) — use `-I` to
> exclude unwanted paths instead.

```bash
odgi viz \
    -i sample_augmented.og \
    -o sample_heatmap.png \
    -O \
    -I "GRCh38" \
    -P \
    -x 2000 -y 300
```

### 9b. Per-Read Depth Coloring (coverage intensity across pangenome)

Colors bins black→blue by mean depth. Shows coverage uniformity.
Uses the display path list to restrict to CHM13 + sample reads only.

```bash
odgi viz \
    -i sample_augmented.og \
    -o sample_depth.png \
    -p display_paths.txt \
    -m \
    -P \
    -x 2000 -y 500
```

### 9c. Clustered View (reads ordered by similarity)

Orders reads by pairwise Jaccard similarity so reads covering the same regions cluster together.
Most structured view — reveals coverage patterns across haplotypes.

> **Note:** `-k` (cluster) cannot be combined with `-p` (path filter).
> Use `-I` to exclude GRCh38 paths instead.

```bash
odgi viz \
    -i sample_augmented.og \
    -o sample_clustered.png \
    -k \
    -I "GRCh38" \
    -P \
    -x 2000 -y 1000
```

### 9d. Strand Coloring (forward=red, reverse=blue)

Applies alignment strand motifs to sample reads only via `-A`.

```bash
odgi viz \
    -i sample_augmented.og \
    -o sample_strand.png \
    -p display_paths.txt \
    -S \
    -A "SAMPLENAME" \
    -P \
    -x 2000 -y 1000
```

---

## Known Limitations & Gotchas

| Issue | Cause | Fix |
|---|---|---|
| `vg augment` fails on GBZ | GBZ is read-only | Use `chr22_mutable.vg` instead |
| `odgi build` fails with id parsing error | Empty path from unmapped read | Run `vg filter -q 1` before augmenting |
| `-O` and `-p` conflict in odgi viz | odgi bug/limitation | Use `-I` to exclude paths instead of `-p` |
| `-k` and `-p` conflict in odgi viz | Documented odgi limitation | Use `-I` instead of `-p` |
| `-s` and `-S` conflict in odgi viz | Mutually exclusive color modes | Run as separate commands |
| Path ordering wrong (e.g. .1, .10, .100) | No zero-padding on read names | Use `%04d` / `%06d` in awk rename step |
| `vg pack -D` gives edge coverage | `-D` = edge dump, not node | Use `-d` (lowercase) for per-base, parse with downstream tools |
| `odgi extract` core dump on augmented graphs | Known odgi bug | Use `vg chunk` + convert instead |
| Short-read indexes cause tail-alignment warnings | Wrong index type for ONT | Rebuild with `vg autoindex --workflow lr-giraffe` |

---

## Coverage Expectations

| Reads | Approx coverage (chr22 ~51Mb) | Useful visualization |
|---|---|---|
| 1,000 reads | ~0.3x | Heatmap only — confirms pipeline works |
| 10,000 reads | ~3x | Depth coloring starts to show patterns |
| 100,000 reads | ~30x | Full analysis — clustered, depth, strand |

For meaningful pangenome visualization with HG001–HG004, target ≥30x coverage per sample on chr22 (~100k reads per sample with ONT R10.4.1 read lengths).

---

## Multi-Sample Extension (HG001–HG004)

Run Parts 3–8 independently per sample, then merge augmented graphs:

```bash
# Per-sample: produces HG001_augmented.vg, HG002_augmented.vg, etc.
for SAMPLE in HG001 HG002 HG003 HG004; do
    # Rename reads
    zcat ${SAMPLE}_chr22.fastq.gz | \
        awk -v s=$SAMPLE 'NR%4==1 {printf "@%s.%06d\n", s, int((NR+3)/4)} NR%4!=1 {print}' | \
        gzip > ${SAMPLE}_renamed.fastq.gz

    # Map
    vg giraffe -b r10 -Z chr22.gbz \
        -d chr22_lr.dist \
        -m chr22_lr.longread.withzip.min \
        -z chr22_lr.longread.zipcodes \
        -f ${SAMPLE}_renamed.fastq.gz \
        -p --threads 8 > ${SAMPLE}.gam

    # Filter
    vg filter -q 1 ${SAMPLE}.gam > ${SAMPLE}_mapped.gam

    # Augment (each sample augments from the base mutable graph)
    vg augment -t 8 chr22_mutable.vg ${SAMPLE}_mapped.gam \
        --label-paths > ${SAMPLE}_augmented.vg
done

# Merge all augmented graphs
# (combine paths from all samples into one graph for joint visualization)
# Note: joint augmentation from a single vg augment call with multiple GAMs
# is preferable when running at scale — consult vg augment --help for options
```
