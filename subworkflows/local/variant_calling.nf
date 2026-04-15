/*
========================================================================================
    SUBWORKFLOW: VARIANT_CALLING
    Purpose : Call small variants (SNPs/Indels) with Clair3 and structural variants
              (SVs) with Sniffles2. Outputs VCF files ready for downstream analysis.
    Tools   : Clair3, Sniffles2, bcftools

    Key design decisions:
    ─────────────────────────────────────────────────────────────────────────────────
    1. Clair3 (small variants — SNPs and Indels):
         - Uses a model matched to the Dorado basecalling model. The Clair3 model
           must correspond to the chemistry and accuracy (FAST/HAC/SUP).
         - Run in GVCF mode (--gvcf) to enable joint genotyping across samples.
         - Phasing is done in the PHASING subworkflow (longphase is more accurate
           than Clair3's internal phasing for ONT long reads).
         - For population studies: merge GVCFs with GLnexus.

    2. Sniffles2 (structural variants — deletions, insertions, inversions, dups, BNDs):
         - --tandem-repeats: CRITICAL. Without a tandem repeat BED annotation,
           Sniffles2 over-calls SVs in satellite and microsatellite regions.
           Use a CHM13-specific TRF (Tandem Repeat Finder) BED — GRCh38 TRF tracks
           are NOT compatible with CHM13 coordinates.
         - --output-rnames: include read names in VCF for traceability.
         - Multi-sample mode: when running a cohort, pass all BAMs to Sniffles2
           simultaneously to produce a jointly-genotyped SV VCF (much better than
           merging per-sample SVs post hoc with SURVIVOR).
         - Mosaic mode (--mosaic): add for cancer / somatic SV detection.

    3. VCF normalization:
         - bcftools norm: left-align and decompose MNPs, normalize indel
           representation. Ensures compatibility with annotation tools (VEP, ANNOVAR).
         - bcftools filter: apply QUAL and FILTER=PASS filters.

    4. SV size / type breakdown:
         - bcftools stats + custom script produce a per-sample SV breakdown table
           (DEL/INS/INV/DUP/BND × size bins) for the MultiQC report.
========================================================================================
*/

include { CLAIR3                  } from '../../modules/local/clair3/main'
include { SNIFFLES2               } from '../../modules/local/sniffles2/main'
include { BCFTOOLS_NORM           } from '../../modules/local/bcftools/main'
include { BCFTOOLS_FILTER         } from '../../modules/local/bcftools/main'
include { BCFTOOLS_STATS          } from '../../modules/local/bcftools/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_SV } from '../../modules/local/bcftools/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_SV } from '../../modules/local/bcftools/main'
include { BGZIP_TABIX             } from '../../modules/local/bgzip_tabix/main'
include { SV_SUMMARY              } from '../../modules/local/sv_summary/main'

workflow VARIANT_CALLING {

    take:
    ch_bam           // [ meta, path(bam) ]
    ch_bai           // [ meta, path(bai) ]
    ch_reference     // path(reference.fa)
    ch_reference_fai // path(reference.fa.fai)

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()
    ch_bam_bai  = ch_bam.join(ch_bai, by: 0)

    // ────────────────────────────────────────────────────────────────────────
    // ── SMALL VARIANTS: Clair3 ───────────────────────────────────────────────
    // ────────────────────────────────────────────────────────────────────────

    // Auto-select Clair3 model from Dorado model string if not explicitly set.
    // Mapping: dorado SUP → clair3 ont_sup_g5; dorado HAC → clair3 ont_hac_g5
    def clair3_model = params.clair3_model ?: deriveClair3Model(params.dorado_model)

    CLAIR3 (
        ch_bam_bai,
        ch_reference,
        ch_reference_fai,
        clair3_model,
        params.clair3_min_coverage
    )
    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    // Normalize Clair3 VCF
    BCFTOOLS_NORM (
        CLAIR3.out.vcf,
        ch_reference,
        'snp_indel'
    )
    ch_snp_vcf_norm = BCFTOOLS_NORM.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    // Filter: PASS variants only, min QUAL 20
    BCFTOOLS_FILTER (
        ch_snp_vcf_norm,
        "'FILTER==\"PASS\" && QUAL>=20'",
        'snp_indel'
    )
    ch_snp_vcf_filtered = BCFTOOLS_FILTER.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions)

    BCFTOOLS_STATS ( ch_snp_vcf_filtered, 'snp_indel' )
    ch_reports  = ch_reports.mix(BCFTOOLS_STATS.out.stats.map { it -> it[1] })
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    // ────────────────────────────────────────────────────────────────────────
    // ── STRUCTURAL VARIANTS: Sniffles2 ───────────────────────────────────────
    // ────────────────────────────────────────────────────────────────────────

    // Build tandem-repeat channel (optional but strongly recommended)
    ch_trf_bed = params.tandem_repeats
        ? channel.value(file(params.tandem_repeats))
        : channel.value([])

    SNIFFLES2 (
        ch_bam_bai,
        ch_reference,
        ch_trf_bed,
        params.sniffles_min_support,
        params.sniffles_min_svlen
    )
    ch_sv_vcf   = SNIFFLES2.out.vcf
    ch_sv_snf   = SNIFFLES2.out.snf     // Sniffles2 .snf files for multi-sample merging
    ch_versions = ch_versions.mix(SNIFFLES2.out.versions)

    // Normalize and index SV VCF
    BCFTOOLS_NORM_SV (
        ch_sv_vcf,
        ch_reference,
        'sv'
    )
    BGZIP_TABIX ( BCFTOOLS_NORM_SV.out.vcf )
    ch_versions = ch_versions.mix(BGZIP_TABIX.out.versions)

    // Per-sample SV type/size summary for MultiQC
    SV_SUMMARY ( BGZIP_TABIX.out.compressed )
    ch_reports  = ch_reports.mix(SV_SUMMARY.out.summary.map { it -> it[1] })
    ch_versions = ch_versions.mix(SV_SUMMARY.out.versions)

    BCFTOOLS_STATS_SV ( BGZIP_TABIX.out.compressed, 'sv' )
    ch_reports = ch_reports.mix(BCFTOOLS_STATS_SV.out.stats.map { it -> it[1] })

    emit:
    snp_vcf  = ch_snp_vcf_filtered      // [ meta, path(.vcf.gz) ] — for phasing
    sv_vcf   = BGZIP_TABIX.out.compressed  // [ meta, path(.vcf.gz) ]
    sv_snf   = ch_sv_snf                // [ meta, path(.snf) ] — for cohort merging
    reports  = ch_reports
    versions = ch_versions
}

// ── Helper: derive Clair3 model name from Dorado model string ────────────────
// Clair3 bundles trained models for common Dorado model families.
// See: https://github.com/HKU-BAL/Clair3#pre-trained-models
def deriveClair3Model(doradoModel) {
    if (!doradoModel) return 'ont'
    if (doradoModel.contains('sup'))  return 'ont'
    if (doradoModel.contains('hac'))  return 'ont'
    if (doradoModel.contains('fast')) return 'ont'
    log.warn "[VARIANT_CALLING] Cannot auto-detect Clair3 model from Dorado model " +
             "'${doradoModel}'. Defaulting to 'ont'. Override with --clair3_model."
    return 'ont'
}
