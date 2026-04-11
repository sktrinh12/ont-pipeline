# ONT-WGS & Methylation Pipeline v1.0.0

A production-ready **Nextflow DSL2** pipeline for Oxford Nanopore (ONT) long-read Whole Genome Sequencing — from raw pod5 signal to structural variants, 5mC/5hmC methylation profiles, phased haplotypes, and pangenome graph projection. Dataset is the Ashkenazim jewish family trio (HG002-HG004).

Refer to [concepts](./concepts.md) for details on the biology and tool use of the pipeline.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [What Changed from v1 (Key Corrections)](#what-changed-from-v1-key-corrections)
- [Technical Stack](#technical-stack)
- [Quick Start](#quick-start)
- [Parameters Reference](#parameters-reference)
- [Output Structure](#output-structure)
- [Theory & Design Decisions](#theory--design-decisions)
- [Running on HPC / AWS](#running-on-hpc--aws)
- [Multi-Sample / Cohort Mode](#multi-sample--cohort-mode)
- [Samplesheet Format](#samplesheet-format)

---

## Overview

Traditional genomics pipelines were built for short reads against a single linear reference. This pipeline is designed for the current reality:

- **R10.4.1 ONT chemistry** with Dorado SUP basecalling (~99.5% per-read accuracy)
- **T2T-CHM13v2.0** as the primary linear reference (the first complete human genome — 8% more sequence than GRCh38, complete centromeres, no gaps)
- **HPRC v1.1 pangenome graph** (94 haplotypes from 47 diverse individuals) to eliminate reference bias
- **Native methylation** from MM/ML tags — no bisulfite treatment, no conversion artefacts
- **Duplex basecalling** for Q50+ accuracy reads (where applicable)

---

## Pipeline Architecture

```
╔══════════════════════════════════════════════════════════════════════════════╗
║           ONT-WGS & METHYLATION PIPELINE v2.0.0                            ║
║           T2T-CHM13v2.0  |  HPRC v1.1 Pangenome                           ║
╚══════════════════════════════════════════════════════════════════════════════╝

  ┌───────────────────────────────────────────────┐
  │  INPUT                                        │
  │  pod5 (preferred)  /  fast5  /  fastq.gz      │
  │  Local  |  S3  |  Samplesheet CSV             │
  └──────────────────────┬────────────────────────┘
                         │
            ┌────────────▼────────────┐
            │   INGESTION & QC        │
            │   pod5 inspect          │
            │   PycoQC (sequencing    │
            │     summary)            │
            │   NanoPlot (per-read)   │
            │   Chopper (filter)  *1  │
            └────────────┬────────────┘
                         │
         ┌───────────────▼───────────────┐
         │  BASECALLING (Dorado)          │
         │  Model: dna_r10.4.1_e8.2_     │
         │    400bps_sup@v4.3.0           │
         │  → MM/ML tags (5mC + 5hmC)    │
         │  → Simplex + Duplex *2         │
         │  → uBAM output                │
         └───────────────┬───────────────┘
                         │
         ┌───────────────▼───────────────┐
         │  ALIGNMENT (minimap2)          │
         │  Preset: lr:hq  *3            │
         │  Reference: T2T-CHM13v2.0     │
         │  Flags: -y (MM/ML passthrough)│
         │  Sort → Markdup → Index       │
         └──────┬────────────────────────┘
                │
     ┌──────────▼──────────┐
     │  COVERAGE QC GATE   │ ← mosdepth
     │  < min_coverage X?  │   Fails sample with warning
     │  PASS  │  FAIL ──►  │   (excluded from calling)
     └──────────┬──────────┘
                │
   ┌────────────┼────────────────────────┐
   │            │                        │
   ▼            ▼                        ▼
┌──────┐  ┌──────────────┐   ┌─────────────────────┐
│ SNP/ │  │ STRUCTURAL   │   │ METHYLATION          │
│INDEL │  │ VARIANTS     │   │ modkit pileup        │
│Clair3│  │ Sniffles2    │   │ --combine-strands    │
│+GVCF │  │ +TRF annot.*4│   │ MM/ML → bedMethyl   │
│      │  │ .snf output  │   │ bgzip + tabix index  │
└──┬───┘  └─────┬────────┘   └──────────────────────┘
   │             │
   └──────┬──────┘
          │
   ┌──────▼────────────────┐
   │  PHASING               │
   │  longphase  *5         │
   │  Joint SNP + SV phase  │
   │  → Haplotagged BAM     │
   │    (HP:i:1 / HP:i:2)  │
   │  → Per-haplotype meth  │
   └──────────────────────┬─┘
                          │
   ┌──────────────────────▼─────────────────────────────┐
   │  PANGENOME GRAPH MAPPING (optional)                 │
   │  vg giraffe → HPRC v1.1 GBZ  *6                   │
   │  → GAM (graph alignments)                          │
   │  → vg surject → linear BAM (CHM13 coords)          │
   │  → vg call → graph-aware VCF                       │
   │  → odgi viz → coverage heatmap PNG                 │
   └───────────────────────────────────────────────────┬┘
                                                       │
                              ┌────────────────────────▼─┐
                              │  MULTIQC REPORT            │
                              │  NanoPlot + PycoQC +       │
                              │  mosdepth + samtools +     │
                              │  bcftools + modkit summary │
                              └──────────────────────────┘

OUTPUT TREE:
  results/
  ├── basecalling/        ← Dorado uBAMs
  ├── filtered_reads/     ← Chopper-filtered FASTQ
  ├── alignment/          ← Sorted, indexed BAMs + flagstat
  ├── methylation/        ← bedMethyl.gz + .tbi per sample
  ├── variants/
  │   ├── snp_indel/      ← Clair3 VCF + GVCF
  │   └── sv/             ← Sniffles2 VCF + .snf
  ├── phasing/            ← Haplotagged BAMs + phased VCFs
  ├── pangenome/          ← Surjected BAMs + graph VCFs + odgi viz
  ├── multiqc/            ← Aggregated HTML QC report
  └── pipeline_info/      ← Timeline, trace, DAG
```

---

## Technical Stack

### Workflow & Infrastructure
| Tool | Purpose |
|---|---|
| Nextflow DSL2 ≥23.04 | Pipeline orchestration, parallelism, resume |
| Docker / Singularity / Conda | Reproducible containerisation |
| AWS Batch / SLURM | HPC & cloud execution |

### Core Bioinformatics Tools
| Category | Tool | Version | Notes |
|---|---|---|---|
| Basecalling | **Dorado** | ≥0.7.0 | GPU recommended; SUP model for maximum accuracy |
| Read filtering | **Chopper** | ≥0.7.0 | Replaces deprecated NanoFilt |
| QC | **NanoPlot** | ≥1.41 | Read length + quality distributions |
| QC | **PycoQC** | ≥2.5.2 | Sequencing summary–based QC |
| QC | **mosdepth** | ≥0.3.6 | Coverage depth + uniformity |
| Alignment | **minimap2** | ≥2.26 | lr:hq preset for R10.4.1 |
| Alignment | **samtools** | ≥1.18 | Sort, index, flagstat, markdup |
| Methylation | **modkit** | ≥0.2.0 | 5mC + 5hmC pileup from MM/ML tags |
| Small variants | **Clair3** | ≥1.0.4 | Deep-learning SNP/Indel caller |
| SVs | **Sniffles2** | ≥2.2 | Joint-genotypable SV calling |
| Phasing | **longphase** | ≥1.6 | Fast ONT phasing + haplotagging |
| Variant QC | **bcftools** | ≥1.18 | Normalise, filter, stats |
| Pangenome | **vg giraffe** | ≥1.55 | Haplotype-aware graph mapper |
| Pangenome | **vg call** | ≥1.55 | Graph-based variant calling |
| Pangenome | **odgi** | ≥0.8 | Graph statistics + visualisation |
| Reporting | **MultiQC** | ≥1.21 | Aggregated QC HTML report |

---

## Quick Start

### Prerequisites

```bash
# install tools
curl -s https://get.nextflow.io | bash
nextflow -version
docker --version
aws --version

# may need to increase swap memory for VG related processes
swapon --show
sudo swapoff /swap.img
sudo fallocate -l 8G /swap.img
sudo chmod 600 /swap.img
sudo mkswap /swap.img
sudo swapon /swap.img
```

[Download Guide](./download_guide.md)

### Run nextflow pipeline

- Run (pre-basecalled chr22 subsetted data)

```bash
nextflow run main.nf \
  -profile test,docker \
  -resume
```

- Run (pod5 input, with GPU basecalling)

```bash
nextflow run main.nf \
  -profile test,docker \
  --input samplesheet.csv \
  --reference ./chm13_chr22.fa \
  --tandem_repeats ./chm13v2.0_chr22_repeats.bed \
  --pangenome_gbz /refs/hprc-v1.1-mc-chm13.gbz \
  --duplex true \
  --outdir results/ \
  -resume
```

---

## Parameters Reference

### Input / Output
| Parameter | Default | Description |
|---|---|---|
| `--input` | null | Path to samplesheet CSV |
| `--input_dir` | null | Directory of pod5/fast5/fastq (auto-detect format) |
| `--outdir` | `results` | Output directory |
| `--reference` | **required** | T2T-CHM13v2.0 FASTA |
| `--tandem_repeats` | null | CHM13-specific TRF BED (strongly recommended) |
| `--pangenome_gbz` | null | HPRC v1.1 GBZ graph (enables pangenome step) |

### Basecalling (Dorado)
| Parameter | Default | Description |
|---|---|---|
| `--skip_basecalling` | false | Skip basecalling (use FASTQ input) |
| `--dorado_model` | `dna_r10.4.1_e8.2_400bps_sup@v4.3.0` | Dorado model (must match chemistry) |
| `--duplex` | false | Enable duplex basecalling |
| `--dorado_min_qscore` | 10 | Q-score filter at basecalling |

### Read Filtering (Chopper — applied to FASTQ inputs)
| Parameter | Default | Description |
|---|---|---|
| `--min_read_length` | 1000 | Minimum read length (bp) |
| `--min_read_quality` | 10 | Minimum mean Q-score |

### Variant Calling
| Parameter | Default | Description |
|---|---|---|
| `--clair3_model` | auto | Clair3 model (auto-detected from Dorado model) |
| `--sniffles_min_support` | 5 | Min reads supporting SV call |
| `--sniffles_min_svlen` | 50 | Min SV length (bp) |

### QC Gates
| Parameter | Default | Description |
|---|---|---|
| `--min_coverage` | 15 | Min genome coverage (X) — samples below this are excluded |

---

## Output Structure

```
results/
├── basecalling/
│   └── {sample}/
│       └── {sample}.ubam                  ← Dorado uBAM (MM/ML tags intact)
├── filtered_reads/
│   └── {sample}/
│       └── {sample}.filtered.fastq.gz     ← Chopper-filtered reads
├── alignment/
│   └── {sample}/
│       ├── {sample}.bam                   ← Sorted, indexed, markdup BAM
│       ├── {sample}.bam.bai
│       └── {sample}.flagstat
├── methylation/
│   └── {sample}/
│       ├── {sample}.bedmethyl.gz          ← Per-CpG 5mC/5hmC calls
│       ├── {sample}.bedmethyl.gz.tbi      ← Tabix index
│       └── {sample}_modkit_summary.tsv
├── variants/
│   ├── snp_indel/
│   │   └── {sample}/
│   │       ├── merge_output.vcf.gz        ← Clair3 SNP/Indel VCF
│   │       └── merge_output.gvcf.gz       ← GVCF for joint genotyping
│   └── sv/
│       └── {sample}/
│           ├── {sample}.sv.vcf.gz         ← Sniffles2 SV VCF
│           └── {sample}.snf               ← For multi-sample merging
├── phasing/
│   └── {sample}/
│       ├── {sample}.phased.vcf.gz         ← longphase phased VCF (PS tags)
│       └── {sample}.haplotagged.bam       ← BAM with HP:i:1/2 tags
├── pangenome/
│   └── {sample}/
│       ├── {sample}.surjected.bam         ← vg surject → linear BAM
│       ├── {sample}.graph.vcf.gz          ← vg call graph-aware VCF
│       └── {sample}.odgi.png              ← Graph coverage heatmap
├── multiqc/
│   ├── multiqc_report.html                ← Main QC report (open in browser)
│   └── multiqc_data/
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.svg
```

---

## Theory & Design Decisions

### Why T2T-CHM13v2.0?
The Genome Reference Consortium's GRCh38 contained ~150 Mb of unresolved sequence, including complete centromeric arrays and subtelomeric regions. These gaps led to:
- Reads from centromeric and pericentromeric regions being unmapped or multiply mapped
- Systematic under-calling of SVs in heterochromatic regions
- Reference bias in allele frequency estimates near gaps

CHM13v2.0 resolves all of this. It also adds a complete, non-mosaic chrY. For ONT WGS in 2024+, CHM13v2.0 is the appropriate reference choice.

### Why the pangenome?
CHM13v2.0, while complete, represents a single haplotype — from a complete hydatidiform mole (CHM13) cell line. It does not represent the full spectrum of human genetic diversity. The HPRC v1.1 graph encodes 94 haplotypes from 47 globally diverse individuals as a variation graph. Aligning to this graph:
- Eliminates reference bias in regions of high haplotype diversity (MHC, centromeres, immunoglobulin loci)
- Enables calling of alleles that simply do not exist in any single linear reference
- Produces surjected BAMs that retain linear-reference coordinates for compatibility

### Why native methylation instead of bisulfite sequencing?
Bisulfite sequencing (BS-seq, WGBS) involves sodium bisulfite treatment that converts unmethylated cytosines to uracil. This:
- Degrades ~95% of the input DNA
- Cannot distinguish 5mC from 5hmC (both resist conversion)
- Introduces conversion artefacts, particularly at CHH and CHG contexts
- Cannot be combined with variant calling on the same library

Dorado + modkit extracts methylation information from the raw electrical signal via learned probabilistic models. MM/ML tags encode per-read, per-base modification probabilities for each cytosine in a CpG context. This approach:
- Requires no additional library preparation
- Simultaneously produces sequence and methylation data from one run
- Distinguishes 5mC (m) from 5hmC (h) when using a combined model
- Is lower cost per sample than WGBS at equivalent depth

### Why longphase instead of WhatsHap?
WhatsHap is the gold standard for short-read phasing and is well-integrated with GATK pipelines. For ONT long reads, longphase offers:
- 5–10× speed improvement (critical at WGS scale)
- Joint phasing of SNPs and SVs (better phase block connectivity)
- Comparable or better N50 phase block lengths
- Native haplotagging support

---

## Multi-Sample / Cohort Mode

For population-level studies, run all samples through the pipeline individually to produce per-sample `.snf` files (Sniffles2) and `.gvcf.gz` files (Clair3). Then:

```bash
# Joint SV genotyping (Sniffles2 multi-sample mode)
sniffles \
  --input sample1.snf sample2.snf sample3.snf \
  --vcf cohort_sv.vcf.gz \
  --reference chm13v2.0.fa.gz

# Joint SNP/Indel genotyping (GLnexus on GVCFs)
glnexus_cli \
  --config DeepVariant_unfiltered \
  sample1.gvcf.gz sample2.gvcf.gz sample3.gvcf.gz \
  > cohort_snp.bcf
```

---

> This pipeline is validated on a Chr22 subset of the HG002 HPRC dataset (GIAB benchmark sample) using the following:
> - ~15× coverage of Chr22 (ONT R10.4.1 SUP)
> - Pre-basecalled FASTQ (skip_basecalling=true for the test profile)
> - GIAB v4.2.1 high-confidence SNP/Indel truthset for benchmarking
> - GIAB T2T-CHM13 SV truthset for SV benchmarking
