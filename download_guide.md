# Data downloads

Essentially, there are **three** ways to cut this cake. Started with the hardest one (Method 1), troubleshoot a "Range Request" (Method 2), and landed on the "S3 Stream" (Method 3).

---

## Summary Comparison

| Metric | Method 1 (Raw FASTQ) | Method 2 (Remote `-X`) | Method 3 (S3 Stream) |
| :--- | :--- | :--- | :--- |
| **Download Size** | 60 GB - 100 GB | ~2 GB | **~2 GB** |
| **RAM Required** | >12 GB | <1 GB | **<1 GB** |
| **Total Time** | 5-10 Hours | 1-2 Hours (unstable) | **10-15 Minutes** |
| **Complexity** | High (Align + Filter) | Medium (Index match) | **Low (One command)** |

**Documentation Recommendation:** "For rapid prototyping (POC), Method 3 (S3 Streaming) is the preferred workflow. Method 1 should only be used when raw FASTQ data is the only available source."

---

## 1. The "Dead End" Log (Initial Attempts)

### **Method 1: The Full-Genome FASTQ Subset**
**The Approach:** Download the entire raw FASTQ (60GB-100GB) via `wget` and use `minimap2` to align it to the whole human genome.
* **Status:** **FAILED / DEAD END**
* **Why:** * **Memory Wall:** Mapping to the full GRCh38 index requires >12GB of RAM. In a standard Docker container, the process is killed (OOM) before it even starts.
    * **Inefficiency:** Spend 10+ hours downloading 100GB of data just to throw away 98% of it (everything not on chr22).
    * **Storage:** Requires massive disk overhead (100GB FASTQ + 100GB BAM + 3GB Ref).

### **Method 2: The Manual Remote Indexing (`-X` flag)**
**The Approach:** Use `samtools view -X` with a remote URL and a local `.crai` index to fetch specific byte ranges.
* **Status:** **INCONSISTENT / SLOW**
* **Why:** * **Handshake Overhead:** Every time `samtools` asks for a new "chunk" via HTTPS, there is a network handshake. For thousands of reads, this results in extreme latency.
    * **CRC Failures:** Network jitters often lead to the `Container header CRC32 failure` saw, causing the stream to crash mid-way.

---

## 2. The "Lightweight" FASTQ Extract

**The Strategy:** Instead of indexing the whole genome, index **only chr22**. 
* **Pros:** Low memory (200MB RAM), fast alignment.
* **Cons:** still had to wait for the 60GB download.

```bash
# download (any) of the large fastq file direct
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_1.fastq.gz

# download grch38 reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# unzip file
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Extract chr22 from reference to make a tiny index
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna chr22 > chr22_ref.fa

# Map the huge file against ONLY that tiny index (~3 hours)
minimap2 -ax map-ont -t 3 chr22_ref.fa GM24385_1.fastq.gz | \
samtools view -@ 4 -b -F 4 - > HG002_chr22_only.bam

# -F 4 removes unmapped reads
samtools view -@ 4 -b -F 4 HG002_chr22_only.bam > HG002_chr22_final.bam

# Now convert this tiny, clean BAM back to FASTQ
samtools fastq HG002_chr22_final.bam > HG002_chr22_final.fastq

# and compress it
gzip HG002_chr22_final.fastq
```

---

## 3. The "Gold Standard" Method (S3 Streaming)

This is the method used for HG001, HG003, and HG004. It is the most professional way to handle public datasets.

**The Strategy:** Access the data via the S3 protocol. This allows `samtools` to "peek" into the file and only pull the exact data packets needed for Chromosome 22.

* **Pros:** Extremely fast (minutes vs hours), no massive downloads, very low storage footprint.
* **Cons:** Requires a specific HTSlib build with S3/libcurl support.

---

# T2T-CHM13 & Pangenome Workflow: Chr22 Subset
**Goal:** Extract Chromosome 22 reads from the HG002 Ashkenazim Son dataset (ONT) and align them to both linear (T2T-CHM13) and pangenome references.

---

## 1. Reference Preparation (T2T-CHM13)
Use the **T2T-CHM13v2.0** assembly. To optimize for POC testing, subset the full reference to Chromosome 22.

```bash
# Download T2T-CHM13 reference
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

# Extract Chr22 and index
samtools faidx chm13v2.0.fa.gz chr22 > chm13_chr22.fa
samtools faidx chm13_chr22.fa
```

---

## 2. Read Extraction (The "Streaming" Method)
This method uses "Range Requests" to download only the byte-ranges corresponding to `chr22` from the S3 bucket.

### Prerequisites: Decoding Reference
CRAM files are reference-compressed. Even though are moving toward CHM13, must use the original **GRCh38** reference to decode the HG002 CRAMs.

```bash
# Download index for remote/local access
curl -C - -O https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/hg002/sup/PAO83395.pass.cram

# while `samtools view -X` allows remote streaming via Range Requests, the overhead of the HTTPS handshake often makes downloading the full file (or a large chunk) faster for repeated trials.
# By using -X, samtools is reading the local `.crai` index to find the exact byte addresses for chr22 on the Amazon servers. It is then making a "Range Request" to only download those specific chunks.
# Decoding: The `required_fields=0x1ff` flag is telling it to use the sequence data embedded directly in the CRAM, bypassing the need for a local GRCh38 fasta file.

samtools view -h -b \
  --input-fmt-option decode_md=0 \
  --input-fmt-option required_fields=0x1ff \
  -X https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/hg002/sup/PAO83395.pass.cram PAO83395.pass.cram.crai \
  chr22 > hg002_chr22.bam

# Extract chr22 reads to BAM using the GRCh38 reference (-T)
samtools view -h -b -T GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -o hg002_chr22_only.bam \
  PAO83395.pass.cram chr22

# then confirm BAM file (confirm read maps only to chr22)
samtools flagstat hg002_chr22_only.bam

# Convert BAM to FASTQ
samtools fastq hg002_chr22_only.bam > hg002_chr22_only.fastq
```

### Execution: S3 Stream to FASTQ

```bash
# GRCh38 reference as the "decoding key" (-T) to pull the chr22 slice from the S3 CRAM
export AWS_NO_SIGN_REQUEST=YES
export CRAM_S3="s3://ont-open-data/giab_2023.05/analysis/hg002/hac/PAO83395.pass.cram"

# Stream chr22 and convert directly to FASTQ
samtools view -@ 8 -u -T GCA_000001405.15_GRCh38_no_alt_analysis_set.fna "$CRAM_S3" chr22 | \
samtools fastq - > hg002_chr22.fastq

# To Download all flow cell samples and combine into one fastq.gz file
# but for test; only downloading one fastq file for each HG sample
export HG_SAMPLE=hg002
export REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# List CRAMs, then for each one: extract chr22 → FASTQ, and concatenate everything
for cram in $(aws s3 ls --no-sign-request --recursive "s3://ont-open-data/giab_2023.05/analysis/${HG_SAMPLE}/hac/" \
    | grep "\.pass\.cram$" \
    | awk '{print "s3://ont-open-data/" $4}'); do
    samtools view -@ 8 -u -T "$REF" "$cram" chr22 | samtools fastq -
done | gzip > hg002_chr22_merged.fastq.gz
```
* **Why `-T`?** The CRAM header contains a hardcoded path to the original author's local directory. `-T` forces `samtools` to use local copy for decompression.

---

## 5. Pangenome Alignment (`vg giraffe`)
Aligning to a graph reduces "reference bias" by accounting for structural variations in the population.

### Step A: Prepare the Graph (HPRC Freeze 1)
```bash
# Download Chr22 graph (.vg)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-chm13/hprc-v1.1-mc-chm13.chroms/chr22.vg

# 1. Create a chopped GFA
# vg mod -X 1024: Breaks long nodes into segments of max 1024bp.
# This is required because the Distance Index (Step 4) cannot handle 
# nodes longer than ~1024bp without crashing or becoming massive
# Use a "Hard GFA Rewrite" to force the reference path into the correct format (`Sample#Haplotype#Contig`).
vg view chr22.vg > chr22.gfa

# Physically rename the path in the GFA to satisfy the HPRC/Ducky metadata parser
# This changes "CHM13#chr22" to "CHM13#0#chr22" (adding the #0# haplotype)
sed 's/CHM13#chr22/CHM13#0#chr22/g' chr22.gfa > chr22.fixed.gfa

# Convert the fixed GFA back to a VG graph using the modern 'convert' tool
# -p specifies the Protobuf (.vg) output format
vg convert -g chr22.fixed.gfa -p > chr22.renamed.vg

# Verify the rename worked (Crucial: Must show CHM13#0#chr22)
vg paths -x chr22.renamed.vg -L | grep "CHM13"

# Chop the graph into 1024bp segments (Recommended for Giraffe/Surject)
vg mod -X 1024 chr22.renamed.vg > chr22.chopped.vg

# 2. Structural Indexing (The "Named" Chain)
# Build the XG first because it acts as the "Translation Table" source for the final GBZ. This preserves the segment-to-node name mapping
vg index -x chr22.chopped.xg chr22.chopped.vg

# Build the GBWT from the XG paths
# -E indexes the embedded reference paths
vg gbwt -x chr22.chopped.xg -E -o chr22.gbwt

# Build the GBZ (The final Giraffe graph)
# --set-reference CHM13 tags the coordinate system correctly
vg gbwt -x chr22.chopped.xg -g chr22.gbz --set-reference CHM13 chr22.gbwt

# Build the Distance Index
vg index -t 16 -j chr22.dist chr22.gbz

# Build the Minimizer Index with Zipcodes
# -z ensures "zipcodes" are stored in a file so Giraffe doesn't rebuild them later
vg minimizer -t 16 -d chr22.dist -o chr22.min -z chr22.zip chr22.gbz

# Map reads using Giraffe (Mapping to Node Space)
# Note: --named-coordinates is REMOVED to prevent the translation error
vg giraffe -t 10 \
    -Z /data/scratch/temp/chr22.gbz \
    -d /data/scratch/temp/chr22.dist \
    -m /data/scratch/temp/chr22.min \
    --zipcode-name chr22.zip \
    -f /data/data/HG001_chr22_extracted.fastq.gz \
    --output-format gam > test_map.gam

# Surject GAM to BAM (Mapping to Reference Space)
# --into-ref CHM13 finds the CHM13#0#chr22 path and provides linear coordinates
vg surject -t 4 \
    -x chr22.gbz \
    -b \
    --into-ref CHM13 \
    test_map.gam > out.surjected.bam
```

---

### 4. Final Validation
Verify your BAM is ready for downstream tools like GATK or Samtools.

```bash
# Check that the header has the correct SN (Sequence Name)
samtools view -H out.surjected.bam | grep "SN"

# Check that reads have valid chromosome names and positions
samtools view out.surjected.bam | head -n 10
```
```


### Step B: Align Reads to Graph
```bash
# Align using GraphAligner
GraphAligner -g chr22.d9.gfa -f hg002_chr22.fastq -a hg002_chr22_pangenome.gaf -x ont
```

---

## Directory Summary

| Category | File | Description |
| :--- | :--- | :--- |
| **Reference** | `chm13_chr22.fa` | Target CHM13 Subset |
| **Decoding** | `GCA_000...analysis_set.fna` | Required to read HG002 CRAMs |
| **Reads** | `hg002_chr22.fastq` | Extracted POC reads (Chr22 only) |
| **Graph** | `chr22.gbz` | Pangenome reference |
| **Alignment** | `hg002_chr22_to_chm13.sorted.bam` | Final Linear Result |
| **Alignment** | `hg002_chr22_pangenome.gaf` | Final Graph Result |

---

**Note on Data Scale:** For a POC, a single flowcell (~4GB FASTQ for chr22) is sufficient. For production-grade variant calling, aggregate all available flowcells listed in the S3 directory to reach $>50\times$ coverage.

#### Fetch T2T-CHM13 v2.0 (hs1) Annotations
- The T2T project uses the "hs1" (Homo Sapiens 1) naming convention for its finalized release. We must subset these global annotations to match chr22 POC.

```bash
# 1. Download Tandem Repeats (via RepeatMasker)
# This identifies SINEs, LINEs, and Simple Repeats for the variant caller.
aws s3 cp --no-sign-request \
  s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed .

# 2. Download Accessibility Mask
# This defines regions of 'low mappability' or 'unreliable' coverage.
aws s3 cp --no-sign-request \
  s3://human-pangenomics/T2T/CHM13/assemblies/annotation/accessibility/hs1.combined_mask.bed .

# 3. Subset both to Chr22 to match Reference/FASTQs
grep "^chr22" chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed > chm13v2.0_chr22_repeats.bed
grep "^chr22" hs1.combined_mask.bed > chm13v2.0_chr22_mask.bed
```
