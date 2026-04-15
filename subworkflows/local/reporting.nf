/*
========================================================================================
    SUBWORKFLOW: REPORTING
    Purpose : Aggregate all QC reports into a single MultiQC HTML report.
    Tools   : MultiQC
========================================================================================
*/

include { MULTIQC } from '../../modules/local/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/local/misc/main'

workflow REPORTING {

    take:
    ch_reports   // Mixed channel of all QC files from upstream subworkflows

    main:
    ch_versions = channel.empty()

    ch_report_files = ch_reports.collect()

    MULTIQC (
        ch_report_files,
        params.multiqc_title ?: 'ONT-WGS Pipeline'
    )

    emit:
    report   = MULTIQC.out.report
    versions = ch_versions
}
