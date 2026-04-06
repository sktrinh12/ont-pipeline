/*
========================================================================================
    SUBWORKFLOW: INGESTION_QC
    Purpose : Validate inputs, run pre-alignment QC on raw ONT data
    Tools   : pod5 inspect, NanoPlot, PycoQC, Chopper (filtering)
========================================================================================
*/

// ── Process imports ───────────────────────────────────────────────────────────
include { POD5_INSPECT      } from '../../modules/local/pod5_inspect/main'
include { NANOPLOT_RAW      } from '../../modules/local/nanoplot/main'
include { PYCOQC            } from '../../modules/local/pycoqc/main'
include { CHOPPER_FILTER    } from '../../modules/local/chopper/main'

workflow INGESTION_QC {

    take:
    ch_samples  // [ meta, path(raw_file_or_dir) ]

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // ── Separate by input type ────────────────────────────────────────────────
    ch_pod5   = ch_samples.filter { meta, f -> meta.input_type == 'pod5' }
    ch_fast5  = ch_samples.filter { meta, f -> meta.input_type == 'fast5' }
    ch_fastq  = ch_samples.filter { meta, f -> meta.input_type =~ /fastq/ }

    // ── Pod5 integrity check + summary statistics ─────────────────────────────
    // pod5 inspect provides read count, signal stats, and detects truncated files
    POD5_INSPECT ( ch_pod5 )
    ch_versions = ch_versions.mix(POD5_INSPECT.out.versions)
    ch_reports  = ch_reports.mix(POD5_INSPECT.out.summary.map { it[1] })

    // ── PycoQC on fast5/sequencing_summary (if available) ────────────────────
    // PycoQC reads sequencing_summary.txt produced by MinKNOW / Dorado
    // It generates per-channel and per-read-length QC plots
    ch_fast5_summary = ch_fast5.map { meta, f ->
        def summary = file("${f.parent}/sequencing_summary.txt")
        summary.exists() ? [ meta, f, summary ] : [ meta, f, [] ]
    }
    PYCOQC ( ch_fast5_summary )
    ch_versions = ch_versions.mix(PYCOQC.out.versions)
    ch_reports  = ch_reports.mix(PYCOQC.out.report.map { it[1] })

    // ── NanoPlot on pre-basecalled FASTQ (if input is already FASTQ) ──────────
    // NanoPlot gives N50, mean Q, read length distributions
    NANOPLOT_RAW ( ch_fastq )
    ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions)
    ch_reports  = ch_reports.mix(NANOPLOT_RAW.out.report.map { it[1] })

    // ── Chopper: length + quality filter on pre-basecalled FASTQ ─────────────
    // Chopper is a Rust rewrite of NanoFilt — much faster, actively maintained.
    // Note: We do NOT filter pod5/fast5 here; that is handled at the Dorado
    //       --min-qscore stage during basecalling.
    CHOPPER_FILTER ( ch_fastq )
    ch_versions = ch_versions.mix(CHOPPER_FILTER.out.versions)

    // ── Validated reads: merge pod5, fast5, and filtered fastq ───────────────
    ch_validated = ch_pod5
        .mix(ch_fast5)
        .mix(CHOPPER_FILTER.out.reads)

    emit:
    validated_reads = ch_validated   // [ meta, path ]
    reports         = ch_reports     // for MultiQC
    versions        = ch_versions
}
