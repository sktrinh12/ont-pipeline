/*
========================================================================================
    SUBWORKFLOW: ALIGNMENT_QC
    Purpose : Align basecalled reads to T2T-CHM13v2.0, sort/index BAM, run
              post-alignment QC, and enforce a coverage depth gate.
    Tools   : minimap2, samtools, mosdepth, NanoPlot (on BAM)

    Key design decisions:
    ─────────────────────────────────────────────────────────────────────────────────
    1. Preset 'lr:hq' (long-read, high-quality) is used instead of the older
       'map-ont' preset. lr:hq was introduced in minimap2 2.26 and is specifically
       tuned for R10.4.1 SUP reads (>Q20), giving better sensitivity for SNPs.

    2. Methylation-aware alignment: when input is a uBAM (MM/ML tags present),
       minimap2 is called with -y to copy tags from the uBAM into the output BAM.
       This is critical — modkit requires MM/ML tags in the aligned BAM.

    3. Mark duplicates is optional (--markdup). For WGS from a single flow cell,
       duplicate rates should be <1%, but marking is best practice.

    4. Coverage gate: mosdepth is run per sample. If coverage < params.min_coverage,
       the sample is flagged and excluded from downstream variant calling with a
       warning. This prevents false-negative SV calls in under-sequenced samples.

    5. T2T-CHM13v2.0 (hs1) replaces GRCh38 as the primary reference. CHM13:
         - Fills 8% of the genome that was previously unresolved in GRCh38
         - Includes complete centromere sequences (important for SV calling)
         - Has no alternate contigs (simpler alignment model)
         - Better representation of African haplotypes
========================================================================================
*/

include { MINIMAP2_ALIGN     } from '../../modules/local/misc/main'
include { SEQKIT_VALIDATE    } from '../../modules/local/misc/main'
include { SAMTOOLS_SORT      } from '../../modules/local/samtools/main'
include { SAMTOOLS_INDEX     } from '../../modules/local/samtools/main'
include { SAMTOOLS_MARKDUP   } from '../../modules/local/samtools/main'
include { SAMTOOLS_FLAGSTAT  } from '../../modules/local/samtools/main'
include { MOSDEPTH           } from '../../modules/local/mosdepth/main'
include { NANOPLOT_BAM       } from '../../modules/local/nanoplot/main'

workflow ALIGNMENT_QC {

    take:
    ch_reads         // [ meta, path(fastq.gz or ubam) ]
    ch_reference     // path(reference.fa)
    ch_reference_fai // path(reference.fa.fai)

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // ── Validate FASTQ reads ───────────────────────────────────────────────────
    // Remove malformed reads where SEQ/QUAL lengths differ to prevent minimap2
    // from failing with "SEQ and QUAL are of different length" errors.
    // Only applies to FASTQ input; skip for uBAM input.
    ch_validated_reads = ch_reads.map { meta, reads ->
        if (reads.toString().endsWith('.bam')) {
            [meta, reads]  // Skip validation for uBAM input
        } else {
            [meta, reads]  // Will be validated below
        }
    }.branch { meta, reads ->
        fastq: !reads.toString().endsWith('.bam')
        ubam:  reads.toString().endsWith('.bam')
    }

    ch_fastq_to_validate = ch_validated_reads.fastq
    ch_ubam = ch_validated_reads.ubam

    SEQKIT_VALIDATE ( ch_fastq_to_validate )
    ch_versions = ch_versions.mix(SEQKIT_VALIDATE.out.versions)

    ch_reads_for_alignment = ch_ubam.mix(SEQKIT_VALIDATE.out.reads)

    // ── Alignment ─────────────────────────────────────────────────────────────
    // -a            : output SAM (samtools converts to BAM inline via pipe)
    // -y            : copy MM/ML tags from uBAM input
    // --secondary=no: suppress secondary alignments for cleaner variant calling
    // -R @RG        : read group tags (required by Clair3)
    MINIMAP2_ALIGN (
        ch_reads_for_alignment,
        ch_reference,
        ch_reference_fai,
        params.minimap2_preset,
        ''
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    // ── Sort → Mark duplicates → Index ───────────────────────────────────────
    SAMTOOLS_SORT ( MINIMAP2_ALIGN.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    if (params.markdup) {
        SAMTOOLS_MARKDUP ( SAMTOOLS_SORT.out.bam )
        ch_sorted_bam = SAMTOOLS_MARKDUP.out.bam
        ch_versions   = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)
    } else {
        ch_sorted_bam = SAMTOOLS_SORT.out.bam
    }

    SAMTOOLS_INDEX ( ch_sorted_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam_bai = ch_sorted_bam.join(SAMTOOLS_INDEX.out.bai, by: 0)

    // ── Post-alignment flagstat ───────────────────────────────────────────────
    SAMTOOLS_FLAGSTAT ( ch_bam_bai.map { meta, bam, bai -> [ meta, bam ] } )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
    ch_reports  = ch_reports.mix(SAMTOOLS_FLAGSTAT.out.stats.map { it[1] })

    // ── NanoPlot on aligned BAM (alignment-specific metrics) ──────────────────
    NANOPLOT_BAM ( ch_bam_bai.map { meta, bam, bai -> [ meta, bam ] }, '' )
    ch_versions = ch_versions.mix(NANOPLOT_BAM.out.versions)
    ch_reports  = ch_reports.mix(NANOPLOT_BAM.out.report.map { it[1] })

    // ── Mosdepth coverage QC ─────────────────────────────────────────────────
    // mosdepth is the fastest coverage tool; it also produces per-window depth
    // BED files useful for identifying low-coverage regions.
    MOSDEPTH (
        ch_bam_bai,
        ch_reference_fai,
        200   // window size in kb
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_reports  = ch_reports.mix(MOSDEPTH.out.summary.map { it[1] })

    // ── Coverage gate ─────────────────────────────────────────────────────────
    // Parse mosdepth summary to get mean genome coverage.
    // Samples below threshold are flagged and removed from downstream steps.
    // Bypass in stub mode or when --skip_coverage_gate is set.
    ch_bam_with_coverage = ch_bam_bai
        .join(MOSDEPTH.out.summary, by: 0)
        .map { meta, bam, bai, summary ->
            def cov = parseMosdepthCoverage(summary)

            // Stub mode bypass: empty files return 0.0, bypass gate
            if (workflow.stubRun) {
                log.warn "[ALIGNMENT_QC] Stub mode: bypassing coverage gate for sample '${meta.id}'"
                return [ meta + [coverage: -1.0], bam, bai ]
            }

            // User-specified bypass
            if (params.skip_coverage_gate) {
                log.warn "[ALIGNMENT_QC] --skip_coverage_gate: bypassing coverage gate for sample '${meta.id}'"
                return [ meta + [coverage: cov], bam, bai ]
            }

            if (cov < params.min_coverage) {
                if (cov == 0.0) {
                    log.warn "[ALIGNMENT_QC] Sample '${meta.id}': coverage parsed as 0.0×. " +
                             "This may indicate an empty mosdepth file or parsing issue. " +
                             "Including sample. To exclude, use --skip_coverage_gate to override."
                    return [ meta + [coverage: cov], bam, bai ]
                }
                log.error "[ALIGNMENT_QC] Sample '${meta.id}': coverage ${cov}× is below " +
                         "minimum threshold of ${params.min_coverage}×. Excluding from " +
                         "variant calling. Check library prep and sequencing yield. " +
                         "Override with --skip_coverage_gate if this is expected."
                return null
            }
            log.info "[ALIGNMENT_QC] Sample '${meta.id}': coverage ${cov}× ✓"
            return [ meta + [coverage: cov], bam, bai ]
        }
        .filter { it != null }

    emit:
    bam      = ch_bam_with_coverage.map { meta, bam, bai -> [ meta, bam ] }
    bai      = ch_bam_with_coverage.map { meta, bam, bai -> [ meta, bai ] }
    bam_bai  = ch_bam_with_coverage
    reports  = ch_reports
    versions = ch_versions
}

// ── Helper: parse mosdepth summary TSV for mean total coverage ───────────────
def parseMosdepthCoverage(summaryFile) {
    try {
        def lines = summaryFile.readLines()
        if (lines.isEmpty()) {
            return 0.0  // Empty file (stub mode)
        }
        def cov = 0.0
        lines.each { line ->
            if (line.startsWith('total\t')) {
                def cols = line.split('\t')
                if (cols.size() >= 4) {
                    cov = cols[3].toFloat()
                }
            }
        }
        return cov
    } catch (Exception e) {
        log.warn "[ALIGNMENT_QC] Failed to parse mosdepth summary '${summaryFile.name}': ${e.message}. Returning 0.0."
        return 0.0
    }
}
