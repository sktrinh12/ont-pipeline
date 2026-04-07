# ONT WGS Pipeline Concepts & Theory

This document explains the biological and computational rationale for each step in the pipeline, connecting the computational processes to the biological questions they answer.

---

## The Big Picture

**The scientific question:** Across chromosome 22 of this Ashkenazi Jewish family, what genetic variants exist, which ones are inherited versus potentially de novo, how are the two copies of chromosome 22 in each person structured, and how do those variants appear when mapped against not just one reference genome but the full spectrum of human genetic diversity?

**Why Chr22:** It's the second smallest autosome (~51 Mb), gene-dense, contains the immunoglobulin lambda locus (IGL), and has historically been one of the hardest chromosomes to assemble due to repetitive content near the centromere.

**Why ONT:** Oxford Nanopore Technology produces reads that are:
- **Long** — routinely spanning 10–100 kb, capable of crossing repetitive regions that stymie short-read sequencers
- **Directly informative** — the electrical signal contains methylation information (5mC, 5hmC) without bisulfite conversion

---

## INGESTION_QC Subworkflow

Validates input data and performs pre-alignment quality control.

### POD5_INSPECT — POD5 File Integrity Check

**What it does technically:** Opens each POD5 file and verifies it's not truncated or corrupted. Extracts read count, signal statistics (mean signal, dwell time), and channel health metrics.

**What it means biologically:** POD5 is the raw electrical signal file from the ONT sequencer — it's what exists before basecalling converts squiggles into DNA letters. Checking POD5 integrity before spending GPU hours basecalling prevents mid-job failures. Similar to checking a film reel isn't broken before developing it.

**Why it may be skipped:** If the input is pre-basecalled FASTQ/CRAM (as in this run), the POD5 files never existed. The pipeline gracefully handles this by filtering input type — POD5 inspection only runs when `meta.input_type == 'pod5'`.

---

### PYCOQC — Sequencing Run QC

**What it does technically:** Reads the `sequencing_summary.txt` file produced by MinKNOW/Dorado during sequencing. Generates per-channel yield plots, read-length distributions over time, and pore health metrics.

**What it would show:** Whether certain pores on the flow cell were underperforming, when quality dropped during the run (flow cells degrade over 72 hours), and which barcodes captured the most reads in a multiplexed run.

**Note:** The `sequencing_summary.txt` file only exists if the sequencer was run locally or the data provider included it. Public datasets (EBI/SRA) almost never include it.

---

### NANOPLOT_RAW — Pre-Alignment Read QC

**What it does technically:** Reads every entry in the FASTQ and computes read length and quality score distributions. Outputs histograms, scatter plots, and NanoStats.txt (N50, median Q-score, total bases).

**Biological significance:** ONT reads range from hundreds of bases to hundreds of kilobases. For chromosome 22 (~51 Mb), longer reads are better because they can span repetitive regions. **N50** is the key metric — if N50 is 15 kb, half the sequenced bases come from reads longer than 15 kb. These reads routinely span multiple genes and bridge structural variants.

**For the trio:** You'd expect similar N50 values across HG002, HG003, HG004 since they were sequenced as part of the same project with the same protocol.

---

### CHOPPER_FILTER — Quality and Length Filtering

**What it does technically:** Streams FASTQ reads through quality and length filters using the Chopper tool (a Rust rewrite of NanoFilt). Removes reads shorter than the configured threshold (default 1000 bp) and below Q10.

**Biological significance:** Very short reads (<1 kb) are often chimeric artifacts — two unrelated DNA fragments ligated together during library prep, or adapter dimers. Including them in alignment creates spurious mappings. Low-quality reads (Q<10) introduce noise into variant calling.

**The trade-off:** Some coverage is lost, but the coverage retained is trustworthy. For fresh off-instrument data, this may remove 10–20% of reads. For published high-quality ONT data (like GIAB), it typically removes <2%.

---

## ALIGNMENT_QC Subworkflow

Maps reads to the reference genome and performs post-alignment QC.

### SEQKIT_VALIDATE — FASTQ Format Validation

**What it does:** Checks that FASTQ files are not truncated, that every read has a matching quality string of the same length, and that no malformed records are present.

**Why it matters:** A truncated FASTQ causes cryptic failures deep in minimap2 or Clair3 that are very hard to debug. Validating at the start costs seconds and saves hours.

---

### MINIMAP2_ALIGN — Read Mapping

**What it does technically:** Maps every read against the T2T-CHM13v2.0 reference using minimap2 with the `lr:hq` preset (long-read, high-quality). Outputs a BAM file where each read is placed at its best matching position with a CIGAR string encoding matches, mismatches, insertions, and deletions.

**Why `lr:hq`:** This preset was introduced in minimap2 2.26 and is specifically tuned for R10.4.1 SUP reads (>Q20). It provides better sensitivity for SNPs compared to the older `map-ont` preset.

**Why T2T-CHM13v2.0:** 
- GRCh38 has gaps and unresolved regions in chromosome 22, particularly near the centromere and in pericentromeric repeats
- CHM13 fills these gaps — it includes complete centromere sequences that are critical for structural variant detection
- CHM13 has no alternate contigs, simplifying the alignment model
- Better representation of African and other non-European haplotypes reduces reference bias

**Read group tags (-R @RG):** Each BAM gets a header tag identifying the sample (SM:HG002, etc.). Downstream tools like Clair3 require this to label VCF columns correctly.

**Biological significance:** This step transforms raw DNA sequences into genomic positions. Every read is anchored to a specific location on chr22. The accuracy of this step determines the accuracy of everything downstream — if a read maps to the wrong place, any variants called from it are wrong.

---

### SAMTOOLS_SORT — Coordinate Sorting

**What it does:** Reorders the BAM so reads are sorted by chromosome position (coordinate sort).

**Why it matters:** Almost every downstream tool — variant callers, coverage tools, genome browsers — requires coordinate-sorted BAM. Random access to any genomic region takes microseconds via the index rather than requiring a full file scan.

---

### SAMTOOLS_MARKDUP — Duplicate Marking

**What it does:** Identifies pairs of reads that are likely PCR duplicates — reads with identical start and end coordinates that were amplified from the same original DNA molecule. Marks them with a flag in the BAM but does not remove them.

**Why it matters for ONT:** PCR duplicates are less common in ONT than Illumina because ONT can sequence native DNA without amplification. For published GIAB data, duplicate rates should be under 1%. Marking allows variant callers to optionally ignore them, preventing artificial inflation of variant confidence at duplicate positions.

---

### SAMTOOLS_INDEX — BAM Indexing

**What it does:** Creates a `.bai` index file alongside the BAM — a lookup table that enables jumping directly to any chromosome:position without reading from the beginning.

**Why it matters:** Without this, every tool that needs to access a specific genomic region would have to scan the entire BAM. With it, accessing any 1 Mb window takes milliseconds.

---

### SAMTOOLS_FLAGSTAT — Alignment Statistics

**What it does:** Counts reads in each alignment category: total reads, mapped reads, properly paired, secondary/supplementary alignments, duplicates. Outputs a plain text summary.

**What the statistics show biologically:**
- **Mapping rate:** What fraction of reads mapped to chr22. For data extracted from a chr22-specific CRAM, expect >99%
- **Supplementary alignments:** Reads that map to two locations, indicating they span a structural variant breakpoint. A high supplementary rate is actually informative — many reads are long enough to bridge SV breakpoints, which is exactly what Sniffles2 requires

---

### NANOPLOT_BAM — Post-Alignment QC

**What it does:** Runs NanoPlot on the aligned BAM rather than the raw FASTQ. Computes alignment-specific metrics: alignment identity, read length vs. quality scatter, coverage uniformity.

**What it adds over raw NanoPlot:** Post-alignment metrics reveal how reads actually performed against the reference. A read that looked Q15 in the FASTQ might have alignment identity of 99.2% — indicating the Q-score model is conservative.

---

### MOSDEPTH — Coverage Analysis

**What it does:** Calculates sequencing depth at every position. Reports mean, median, and per-window depth with cumulative distributions.

**Biological significance:** Coverage depth directly determines the ability to call variants confidently:
- At 10× coverage, both alleles of a heterozygous variant are expected ~50% of the time each
- At 30×, statistical confidence is much higher
- For structural variants, Sniffles2 requires a minimum of 5 reads supporting a variant — so regions with <5× coverage will have no SV calls

**The coverage gate:** The pipeline checks whether each sample meets the minimum coverage threshold (default 15×, bypassed in test mode). Samples failing this gate are excluded from downstream variant calling.

---

## VARIANT_CALLING Subworkflow

Detects genetic variants — SNPs, indels, and structural variants.

### CLAIR3 — SNP and Indel Calling

**What it does technically:** A deep learning model trained specifically on ONT reads that identifies SNPs and small indels (up to ~50 bp) from the aligned BAM. Runs in GVCF mode, outputting a record for every position.

**What it means biologically:** This is a per-sample catalogue of every position on chr22 where each person differs from the CHM13 reference.

**Why the trio matters:** A variant in HG002 (son) that is also in HG003 (father) or HG004 (mother) is likely inherited. A variant in HG002 that appears in neither parent is a candidate de novo mutation — the biologically most interesting calls.

**GVCF mode:** Outputs a record for every position, even positions with no variant. This is needed for joint genotyping — when HG002 has a variant at position X but HG003 has a GVCF record showing 0/0, the analysis can confidently state HG003 does not have that variant.

**Model selection:** Clair3 uses models matched to Dorado basecalling accuracy (SUP/HAC/FAST). The pipeline auto-derives the appropriate model from the Dorado model string.

---

### BCFTOOLS_NORM (snp_indel) — VCF Normalization

**What it does:** Left-aligns indels and splits multiallelic records. For example, a deletion represented as `ATCG→A` at position 100 gets normalized to a canonical form.

**Why it matters:** Without normalization, the same biological deletion in HG002 and HG003 might be written differently and appear as non-matching when comparing. Every downstream comparison tool requires normalized VCFs.

---

### BCFTOOLS_FILTER (snp_indel) — Quality Filtering

**What it does:** Keeps only variants where `FILTER=PASS` and `QUAL≥20`. Removes low-confidence calls.

**Biological significance:** Clair3 assigns every variant a quality score and filter status. Variants in repetitive regions or with low read support get flagged. This step removes noise, keeping only trustworthy variants.

---

### BCFTOOLS_STATS (snp_indel) — Variant Statistics

**What it does:** Generates statistics — SNP count, indel count, Ti/Tv ratio (transitions vs. transversions), heterozygous/homozygous ratio.

**The Ti/Tv ratio:** A key QC metric. In human genomes, true variants have Ti/Tv ~2.0–2.1 for whole genome (higher for exomes due to CpG enrichment). If the ratio is 1.5, there are too many false positive transversions. If it's 3.0+, something is wrong.

---

### SNIFFLES2 — Structural Variant Calling

**What it does technically:** Unlike Clair3 which looks at individual base mismatches, Sniffles2 looks at patterns across entire reads to find large structural variants — deletions (DEL), insertions (INS), inversions (INV), duplications (DUP), and breakend translocations (BND). Uses CIGAR strings and SA tags (supplementary alignments) to identify reads that cross SV breakpoints.

**Biological significance for chr22:**
- **Immunoglobulin lambda (IGL) locus** (~900 kb) is a tandemly duplicated gene cluster where CNV is common and functionally important for antibody diversity. Standard short reads cannot reliably genotype CNVs here — ONT reads can.
- **22q11.2 region** contains one of the most common microdeletion syndromes (DiGeorge syndrome, 1 in 4,000 births). SVs here are clinically significant.
- **Centromeric/pericentromeric SVs** that were invisible in GRCh38 are now resolvable with CHM13.

**The .snf output:** Each sample produces a binary `.snf` file encoding all read-level SV evidence. These can be merged post-pipeline with `sniffles --input *.snf` for joint genotyping.

**The tandem-repeat flag:** Without a TRF BED file, Sniffles2 over-calls SVs in satellite and microsatellite regions. The pipeline uses a CHM13-specific TRF annotation — GRCh38 TRF tracks are NOT compatible with CHM13 coordinates.

---

### BCFTOOLS_NORM_SV, BGZIP_TABIX, SV_SUMMARY, BCFTOOLS_STATS_SV

**What they do collectively:** Normalize the SV VCF, compress and index it (so IGV can display it), generate a per-sample table of SV counts by type and size, and compute VCF-level statistics.

**SV_SUMMARY table:** Particularly useful for the trio — comparing SV type/size distributions across HG002, HG003, HG004 reveals whether the son's SV landscape looks like a combination of the parents or has outlier counts worth investigating.

---

## PHASING Subworkflow

Resolves which variants are on the same chromosome copy (haplotype).

### LONGPHASE_PHASE — Haplotype Phasing

**What it does technically:** Takes the SNP VCF from Clair3 and the aligned BAM, determines which variants are on the same physical chromosome copy (haplotype). A single ONT read spanning two heterozygous positions indicates those alleles are in cis (same chromosome). By chaining these connections across thousands of variants, longphase builds "phase blocks."

**Algorithm:** Longphase uses read-backed phasing:
1. For each heterozygous variant, reads spanning multiple het sites connect variants into phase blocks (PS tags)
2. With long ONT reads (N50 >10 kb), individual reads routinely span dozens of het sites, producing very long phase blocks (N50 >10 Mb)
3. SVs can be included via `--sv-file`; junction-spanning reads link SVs to SNPs, dramatically improving phase block continuity in SV-dense regions

**Biological significance:** Every human has two copies of chr22 — one from each parent. Standard variant calling indicates *that* a person is heterozygous at a position but not *which* copy carries the alternate allele. Phasing answers that question.

**Phase block N50:** The key output metric. If N50 is 8 Mb on chr22 (~51 Mb), half the phased bases are in blocks longer than 8 Mb — most of chr22 is resolved into two continuous haplotypes. Short-read phasing achieves N50 ~100 kb. ONT long reads achieve N50 >10 Mb routinely.

---

### PHASE_STATS — Phase Block Statistics

**What it does:** Parses the phased VCF and computes N50 phase block length, number of phase blocks, total phased bases, and median block size.

**Biological interpretation for chr22 at 30× ONT coverage:**
- N50 > 5 Mb = excellent, chr22 mostly resolved into two haplotypes
- N50 1–5 Mb = good, some fragmentation at repetitive regions
- N50 < 1 Mb = concerning, suggests low heterozygosity or coverage issues

---

### LONGPHASE_HAPLOTAG — Read Haplotype Assignment

**What it does:** Goes through every read in the BAM and assigns it to a haplotype — `HP:i:1`, `HP:i:2`, or `HP:i:0` (unphased). Writes these tags into the BAM itself.

**What this enables biologically:** The haplotagged BAM is the gateway to allele-specific analysis:
- **Allele-specific methylation (ASM):** Run modkit with `--partition-tag HP` to see whether haplotype 1 and haplotype 2 have different methylation patterns at the same CpG site. Imprinted genes (where only one parental copy is expressed) show this.
- **Allele-specific expression:** Matching phased variants to RNA-seq data shows which allele is transcribed
- **Visual inspection in IGV:** Coloring by HP tag shows two separate streams of reads — one per chromosome copy — at any heterozygous position

**Expected haplotagging rate:** For ONT WGS at 30× with R10.4.1 reads, ~85–92% of mapped reads are assigned to a haplotype. Unassigned reads (HP:i:0) are typically in homozygous regions or repetitive elements where SNP density is insufficient for phasing.

---

## PANGENOME Subworkflow

Maps reads to a variation graph representing human genetic diversity.

### VG_GIRAFFE — Pangenome Mapping

**What it does technically:** Maps reads to the HPRC v1.1 pangenome variation graph using the haplotype-aware giraffe algorithm. This graph encodes 94 haplotypes from 47 diverse humans as a network of nodes and edges.

**Why pangenome?** The linear CHM13 reference represents one person's genome. The HPRC graph represents the known genetic diversity of humanity at every locus. When reads are mapped to the graph, they naturally follow paths corresponding to the sample's specific genotype — including alleles that don't exist in CHM13 but exist in other HPRC individuals.

**Reference bias:** Reads from individuals whose haplotypes diverge significantly from CHM13 (e.g., African populations) are forced to align to a foreign sequence in linear alignment, introducing systematic bias. Pangenome mapping reduces this bias.

**Index files required:**
- `.gbz` — GBZ-format graph (sequence + topology + haplotypes)
- `.min` — minimizer index (k-mer seeds for fast mapping)
- `.dist` — distance index (for gap costs in seed extension)
- `.zip` — zipcode index (for named-coordinate reporting)

**Memory requirements:** The full HPRC v1.1 graph requires ~300 GB RAM, making it impractical for most laptop or desktop environments. However, the pipeline supports chromosome-level subgraphs. By extracting only chromosome 22 from the full graph, the memory requirements drop dramatically — enabling pangenome analysis on consumer hardware. This approach was used in the test run, where only the chr22.pangenome.vg file was downloaded and converted to gbz, min, dist, and zip formats.

**For non-European samples:** Giraffe pangenome mapping typically recovers 1–3% more reads than minimap2 to CHM13 alone.

---

### VG_STATS — Pangenome Mapping Statistics

**What it does:** Reports mapping statistics from the GAM file — total reads, percentage mapped, mapping quality distribution, perfect alignments.

**What to look for:** Compare mapping rate here to minimap2 mapping rate. If vg giraffe maps 1–3% more reads, that's the pangenome advantage — reads in regions of high diversity that couldn't find a clean home in CHM13 found their correct path in the graph.

---

### VG_SURJECT — Graph-to-Linear Projection

**What it does:** Projects each graph alignment back onto the linear CHM13 coordinate system, producing a standard BAM.

**Why this step exists:** The bioinformatics ecosystem (IGV, bcftools, most variant callers) speaks linear-reference BAM, not graph GAM. Surjection provides the best of both worlds — reads are mapped with pangenome awareness, but the output is compatible with standard tools.

**Comparison:** The surjected BAM can be compared with the minimap2 BAM side by side in IGV — surjected BAM performs better in regions with high haplotype diversity.

---

### VG_CALL — Graph-Based Variant Calling

**What it does:** Calls variants directly from the graph alignments. Identifies which paths through the pangenome graph are supported by each sample's reads — including paths that encode known population-specific alleles.

**What this gives that Clair3 doesn't:** VG_CALL can identify variants pre-encoded in the pangenome graph from the other 46 HPRC individuals. If a 2 kb insertion exists in 10 of the 94 HPRC haplotypes but not in CHM13, Clair3 struggles to call it (it appears as soft-clipped reads with no reference). VG_CALL identifies it cleanly because the insertion already exists as a graph path.

**For the trio:** Comparing VG_CALL results to Clair3 results reveals which variants are "novel to CHM13" vs. "already known in the pangenome." Variants only in VG_CALL are pangenome-aware calls; variants in both are robustly confirmed.

---

## REPORTING:MULTIQC — Results Aggregation

**What it does:** Aggregates all QC outputs into a single interactive HTML report showing:
- Read quality and length distributions (NanoPlot)
- Alignment rates and read category breakdown (samtools flagstat)
- Coverage depth across chr22 (mosdepth)
- Variant counts, Ti/Tv ratio, indel size distribution (bcftools stats)
- SV type/size breakdown per sample (SV_SUMMARY)
- Phase block N50 per sample (PHASE_STATS)
- Pangenome mapping rate vs. linear mapping rate (VG_STATS)

**The biological story in one page:** From top to bottom — the reads were good quality, they mapped well, coverage was sufficient, each person has ~3–4 million SNPs and ~10,000 SVs on chr22, the chromosome is largely phased into two continuous haplotypes, and the pangenome recovered additional reads beyond linear mapping.

---

## What Is Missing Without POD5 (Methylation)

The only subworkflow that does not run without POD5 files is **methylation**. Conceptually, what is lost:

Chromosome 22 carries several imprinted regions where only one parental copy is methylated (silenced) and the other is active. With POD5 → Dorado → modkit, per-CpG methylation calls are obtained for all four samples, enabling:

- **Imprinting analysis:** Does HG002's methylation pattern at imprinted loci reflect the expected maternal vs. paternal silencing?
- **Allele-specific methylation:** Is there ASM at non-imprinted sites where the two haplotypes have different methylation?
- **Epigenetic inheritance:** Do HG003 and HG004 show consistent methylation at the same sites that could be passed to HG002?

That analysis requires the raw signal and cannot be recovered from FASTQ. Everything else in the pipeline captures the DNA sequence layer completely — the methylation layer is simply absent from the FASTQ data format.

---

## Pipeline Flow Summary

| Stage | Biological Question | Key Tools |
|-------|---------------------|-----------|
| INGESTION_QC | Is the sequencing data high quality? | pod5_inspect, NanoPlot, Chopper |
| ALIGNMENT_QC | Where do the reads map? | minimap2, samtools, mosdepth |
| VARIANT_CALLING | What genetic variants exist? | Clair3 (SNPs/indels), Sniffles2 (SVs) |
| PHASING | Which variants are on the same chromosome copy? | longphase |
| PANGENOME | What population-specific variants are missed by linear reference? | vg giraffe, vg call |
| REPORTING | Is everything working correctly? | MultiQC |
