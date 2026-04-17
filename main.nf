#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
=======================================================================================
    SCALABLE ONT-WGS PIPELINE (T2T + Pangenome)
=======================================================================================
    Author      : ONT Genomics Pipeline
    Version     : 1.0.0
    Description : End-to-end Oxford Nanopore WGS pipeline with SV calling,
                  5mC/5hmC methylation profiling, and pangenome graph projection.
    Reference   : T2T-CHM13v2.0 (linear) + HPRC v1.1 (pangenome)
=======================================================================================
*/

// ── Import subworkflows ──────────────────────────────────────────────────────
include { INGESTION_QC       } from './subworkflows/local/ingestion_qc'
include { ALIGNMENT_QC       } from './subworkflows/local/alignment_qc'
include { VARIANT_CALLING    } from './subworkflows/local/variant_calling'
include { PHASING            } from './subworkflows/local/phasing'
include { PANGENOME          } from './subworkflows/local/pangenome'
include { REPORTING          } from './subworkflows/local/reporting'

// ── Validate parameters ──────────────────────────────────────────────────────
def validateParams() {
    def errors = []

    // Input validation
    if (!params.input && !params.input_dir) {
        errors << "ERROR: Provide --input (samplesheet CSV) or --input_dir (directory of pod5/fast5/fastq)"
    }
    if (params.input && !file(params.input).exists()) {
        errors << "ERROR: Samplesheet not found: ${params.input}"
    }
    if (!params.reference) {
        errors << "ERROR: --reference (T2T-CHM13v2.0 FASTA) is required"
    }

    // Hard-stop on errors
    if (errors) {
        errors.each { e -> log.error(e) }
        System.exit(1)
    }
}

// ── Print pipeline header ────────────────────────────────────────────────────
def printHeader() {
    log.info """
    ╔═══════════════════════════════════════════════════════════════╗
    ║         ONT-WGS PIPELINE  v${workflow.manifest.version}            ║
    ╠═══════════════════════════════════════════════════════════════╣
    ║  Input          : ${params.input ?: params.input_dir}
    ║  Reference      : ${params.reference}
    ║  Pangenome GBZ  : ${params.pangenome_gbz ?: 'DISABLED'}
    ║  Phasing        : ${!params.skip_phasing}
    ║  Output Dir     : ${params.outdir}
    ╚═══════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

// ── Build sample channel from samplesheet or directory scan ─────────────────
def buildSampleChannel() {
    if (params.input) {
        // Samplesheet: sample_id,input_path,condition,family_id,barcode
        return channel
            .fromPath(params.input)
            .splitCsv(header: true, strip: true)
            .map { row ->
                def inputPath = row.input_path
                if (!inputPath) return null
                def meta = [
                    id        : row.sample_id,
                    condition : row.condition ?: 'unknown',
                    family_id : row.family_id ?: 'unknown',
                    barcode   : row.barcode ?: 'unknown',
                    input_type: detectInputType(inputPath)
                ]
                [ meta, file(inputPath) ]
            }
            .filter { it -> it != null }
    } else {
        // Auto-detect pod5/fast5/fastq in directory
        def exts = params.skip_basecalling
            ? "**/*.{fastq,fastq.gz,fq,fq.gz}"
            : "**/*.{pod5,fast5}"
        return channel
            .fromPath("${params.input_dir}/${exts}")
            .map { f ->
                def meta = [
                    id        : f.simpleName,
                    condition : 'unknown',
                    input_type: detectInputType(f.toString())
                ]
                [ meta, f ]
            }
            .groupTuple(by: 0) // group by meta to handle multi-file samples
    }
}

def detectInputType(path) {
    if (!path) return 'unknown'
    if (path.endsWith('.pod5'))                              return 'pod5'
    if (path.endsWith('.fast5'))                            return 'fast5'
    if (path.endsWith('.fastq') || path.endsWith('.fq'))    return 'fastq'
    if (path.endsWith('.fastq.gz') || path.endsWith('.fq.gz')) return 'fastq_gz'
    return 'unknown'
}

// ── MAIN WORKFLOW ────────────────────────────────────────────────────────────
workflow {

    validateParams()
    printHeader()

    ch_samples = buildSampleChannel()

    ch_reference    = file(params.reference)
    ch_reference_fai = file("${params.reference}.fai")

    // ── STEP 1: Pre-alignment QC (runs on raw pod5/fast5 or pre-basecalled fastq)
    INGESTION_QC ( ch_samples )
    ch_qc_reports = INGESTION_QC.out.reports

    // ── STEP 2: Basecalling (Dorado simplex + optional duplex)
    // Gate: skip if input is already FASTQ
    ch_for_basecalling = params.skip_basecalling
        ? channel.empty()
        : INGESTION_QC.out.validated_reads.filter { meta, _f -> meta.input_type in ['pod5','fast5'] }

    // Use pre-basecalled FASTQ directly
    ch_basecalled = INGESTION_QC.out.validated_reads
        .filter { meta, _f -> meta.input_type =~ /fastq/ }

    // ── STEP 3: Alignment to T2T-CHM13v2.0 + post-alignment QC
    ALIGNMENT_QC (
        ch_basecalled,
        ch_reference,
        ch_reference_fai
    )
    ch_bam       = ALIGNMENT_QC.out.bam
    ch_bai       = ALIGNMENT_QC.out.bai
    ch_qc_reports = ch_qc_reports.mix(ALIGNMENT_QC.out.reports)

    // ── STEP 5: Variant calling (SNPs/Indels via Clair3, SVs via Sniffles2)
    VARIANT_CALLING (
        ch_bam,
        ch_bai,
        ch_reference,
        ch_reference_fai
    )
    ch_qc_reports = ch_qc_reports.mix(VARIANT_CALLING.out.reports)

    // ── STEP 6: Phasing with longphase + haplotagging
    if (!params.skip_phasing) {
        PHASING (
            ch_bam,
            ch_bai,
            VARIANT_CALLING.out.snp_vcf,
            ch_reference,
            ch_reference_fai
        )
        ch_qc_reports = ch_qc_reports.mix(PHASING.out.reports)
    }

    // ── STEP 7: Pangenome graph mapping (vg giraffe → HPRC v1.1 GBZ)
    if (params.pangenome_gbz) {
        PANGENOME ( ch_basecalled, file(params.pangenome_gbz) )
        ch_qc_reports = ch_qc_reports.mix(PANGENOME.out.reports)
    }

    // ── STEP 8: Aggregate MultiQC report
    REPORTING ( ch_qc_reports )
}

// ── On completion ────────────────────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    log.info """
    ═══════════════════════════════════════════════════════
     Pipeline completed  : ${status}
     Duration            : ${workflow.duration}
     Output directory    : ${params.outdir}
     Trace report        : ${params.outdir}/pipeline_info/trace.txt
    ═══════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}
