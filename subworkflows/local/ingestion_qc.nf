/*
========================================================================================
    SUBWORKFLOW: INGESTION_QC
    Purpose : Validate inputs, run pre-alignment QC on raw ONT data
    Tools   : pod5 inspect, NanoPlot, PycoQC, Chopper (filtering)
========================================================================================
*/

// ── Process imports ───────────────────────────────────────────────────────────
include { NANOPLOT_RAW      } from '../../modules/local/nanoplot/main'
include { CHOPPER_FILTER    } from '../../modules/local/chopper/main'

workflow INGESTION_QC {

    take:
    ch_samples  // [ meta, path(raw_file_or_dir) ]

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    // ── Separate by input type ────────────────────────────────────────────────
    ch_pod5   = ch_samples.filter { meta, _f -> meta.input_type == 'pod5' }
    ch_fast5  = ch_samples.filter { meta, _f -> meta.input_type == 'fast5' }
    ch_fastq  = ch_samples.filter { meta, _f -> meta.input_type =~ /fastq/ }

    // ── NanoPlot on pre-basecalled FASTQ (if input is already FASTQ) ──────────
    // NanoPlot gives N50, mean Q, read length distributions
    NANOPLOT_RAW ( ch_fastq )
    ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions)
    ch_reports  = ch_reports.mix(NANOPLOT_RAW.out.report.map { it -> it[1] })

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
