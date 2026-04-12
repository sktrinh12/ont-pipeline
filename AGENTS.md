# AGENTS.md — ONT-WGS Pipeline

## Overview

Nextflow DSL2 pipeline for Oxford Nanopore long-read WGS: basecalling, alignment, methylation profiling (5mC/5hmC), variant calling (SNPs/Indels/SVs), phasing, optional pangenome graph projection.

## Build/Lint/Test Commands

### Running the Pipeline

```bash
# Test profile (Chr22 HG00X subset, no GPU)
nextflow run main.nf -profile test,docker -resume
```

### Validation & Linting

```bash
# Check syntax and logic (always use -stub, NOT -validate)
nextflow run main.nf -stub

# Lint (check for issues)
nextflow lint .
```

---

## Testing Changes

When modifying modules, subworkflows, or processes:

```bash
# 1. Syntax check (always use -stub)
nextflow run main.nf -stub

# 2. Lint specific files
nextflow lint modules/local/multiqc/main.nf
nextflow lint subworkflows/local/reporting.nf

# 3. Run single subworkflow for targeted testing
nextflow run main.nf -profile test,docker -resume -entry INGESTION_QC

# 4. Full test profile run
nextflow run main.nf -profile test,docker -resume
```

**Important:**
- **NEVER** use `nextflow validate` or `-validate` — it is not a valid command
- **ALWAYS** use `nextflow run main.nf -stub` to check syntax and logic
- Execute Nextflow commands from `/home/spencer-trinh/Downloads/ont-wgs-pipeline`

### Debugging

```bash
# View timeline/report/trace
nextflow run main.nf -with-timeline timeline.html
nextflow run main.nf -with-report report.html
nextflow run main.nf -with-trace trace.txt

# Dry run (no execution)
nextflow run main.nf -dryrun
```

Never delete a file without asking permission; always prompt before deleting.

---

## Code Style Guidelines

### Language & File Structure
- **Nextflow DSL2** — `nextflow.enable.dsl = 2` at top of main.nf
- All `.nf` files follow DSL2 syntax

```
main.nf                    # Main workflow entry point
nextflow.config           # Configuration, profiles, params
nextflow_schema.json      # JSON schema for params validation
modules/local/<tool>/main.nf   # One process per file
subworkflows/local/*.nf   # Composite workflows
assets/                   # MultiQC config, templates
test/                     # Test data (samplesheet_test.csv, refs/)
```

### Process Naming
- Use UpperCamelCase: `SAMTOOLS_SORT`, `DORADO_BASECALLER`, `MINIMAP2_ALIGN`
- One process per file (nf-test compatibility)
- All processes must be in `modules/local/`.
- File name matches process: `modules/local/samtools/main.nf` contains `process SAMTOOLS_SORT`

### Subworkflow Naming
- Use UpperCamelCase: `INGESTION_QC`, `BASECALLING`, `ALIGNMENT_QC`
- Include in main.nf: `include { WORKFLOW_NAME } from './subworkflows/local/workflow_name'`

### Parameters
- Lowercase with underscores: `params.input`, `params.dorado_model`, `params.skip_methylation`
- Default values in `nextflow.config` under `params {}` block
- Required params validated in `validateParams()` in main.nf

### Input/Output Declarations

```nextflow
process PROCESS_NAME {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/subdir", mode: 'copy'
    container 'image:tag'

    input:
    tuple val(meta), path(bam)
    val config_value
    path reference

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    command --option ${prefix}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version 2>&1 | head -1 | awk '{print $2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: "stub"
    END_VERSIONS
    """
}
```

### Channel Naming
- Prefix with `ch_`: `ch_samples`, `ch_bam`, `ch_reports`
- Use `.mix()` to combine, `.filter()` to subset, `.map()` to transform

### Metadata Structure
- Standard: `[ id: sample_id, condition: condition_label, input_type: 'pod5|fast5|fastq' ]`
- Access in script: `$meta.id`, `$meta.condition`

### Error Handling
- In `nextflow.config`: retry on OOM/kill signals (104, 134, 137, 139, 143, 247)
- `maxRetries = 2` by default
- Use `errorStrategy = 'terminate'` for critical processes

### Version Tracking
- Every process MUST emit `versions.yml` with tool version
```yaml
"PROCESS_NAME":
    toolname: "1.2.3"
```
- Use stub block for testing (empty output files with `touch`)

### Labels (Resource Classes)
- `process_low`: 4 cpus, 5GB, 2h (QC, indexing)
- `process_medium`: 8 cpus, 32GB, 8h (sorting, filtering)
- `process_high`: 32 cpus, 128GB, 48h (variant calling)

### Code Style
- Use `def` for variable declarations
- Groovy interpolation: `${variable}` or `"${meta.id}.bam"`
- Triple-backslash `\\` for command continuation
- Prefer `task.ext.prefix` over hardcoded naming
- Use `stripIndent()` for multi-line log messages

### Samplesheet Format
```csv
sample_id,input_path,condition,family_id,barcode
HG002,/home/spencer-trinh/Downloads/ont-wgs-pipeline/data/HG002_chr22_extracted.fastq.gz,proband,AJ_trio,PAO83395
```

### Testing Strategy
- Test profile uses Chr22 subset, relaxed thresholds (min_coverage=5)
- Test profile skips basecalling (pre-basecalled FASTQ) and methylation since those steps require the original pod5 ONT raw data files
- Test data in `test/` directory
- Use stub blocks for process validation
