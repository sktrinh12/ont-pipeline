/*
========================================================================================
    SUBWORKFLOW: PHASING
    Purpose : Phase heterozygous SNPs/Indels and SVs with longphase, then haplotag
              the BAM so that reads are assigned to haplotype 1 or 2 (HP:Z:1/2).
    Tools   : longphase, samtools (haplotag)

    Key design decisions:
    ─────────────────────────────────────────────────────────────────────────────────
    1. longphase vs WhatsHap:
         - WhatsHap is widely used but is ~5–10× slower than longphase for ONT.
         - longphase is specifically designed for long-read phasing and achieves
           comparable or better phase block lengths (N50 often >10 Mb on WGS).
         - longphase also jointly phases SVs alongside SNPs, improving haplotype
           block connectivity in SV-dense regions.
         - Use WhatsHap if you need: VCF-based phasing with pre-existing haplotype
           blocks, or strict GATK compatibility.

    2. Haplotagging: after phasing, we use `longphase haplotag` (or samtools phase
       alternative) to write HP:i:1 / HP:i:2 / HP:i:0 (unphased) tags into the BAM.
       The haplotagged BAM is the key output — it enables:
         - Per-haplotype methylation calling (modkit --partition-tag HP)
         - Allele-specific expression analysis
         - De novo assembly of individual haplotypes

    3. Phase block statistics:
         - We compute N50 phase block length using a helper script.
         - Typical values: >10 Mb on 30×+ ONT WGS with R10.4.1 SUP.
         - Values below 1 Mb suggest insufficient heterozygosity or coverage.

    4. Homozygous samples: longphase will still run but phase blocks are trivial.
       Expected for highly inbred lines or when using a reference-matched sample.
========================================================================================
*/

include { LONGPHASE_PHASE      } from '../../modules/local/longphase/main'
include { LONGPHASE_HAPLOTAG   } from '../../modules/local/longphase/main'
include { SAMTOOLS_INDEX        } from '../../modules/local/samtools/main'
include { PHASE_STATS           } from '../../modules/local/phase_stats/main'

workflow PHASING {

    take:
    ch_bam           // [ meta, path(bam) ]
    ch_bai           // [ meta, path(bai) ]
    ch_snp_vcf       // [ meta, path(snp.vcf.gz) ] — from Clair3 (filtered, PASS only)
    ch_reference     // path(reference.fa)
    ch_reference_fai // path(reference.fa.fai)

    main:
    ch_versions = channel.empty()
    ch_reports  = channel.empty()

    ch_bam_bai = ch_bam.join(ch_bai, by: 0)
    ch_phase_input = ch_bam_bai.join(ch_snp_vcf, by: 0)

    // ── longphase phase ───────────────────────────────────────────────────────
    // longphase reads the BAM, the SNP VCF (must be phased-ready: het sites only),
    // and outputs a phased VCF with PS (phase set) tags and HP haplotype blocks.
    //
    // Key flags:
    //   --ont          : use ONT long-read phasing parameters
    //   --sv-file      : (optional) include SV VCF for joint SV+SNP phasing
    //                    This dramatically improves phase blocks in SV-rich regions
    //   --indels       : include indels in phasing (default: SNPs only)
    //   --min-snp-af   : minimum allele frequency to include a site (default 0.15)
    LONGPHASE_PHASE (
        ch_phase_input,
        ch_reference,
        ch_reference_fai,
        params.longphase_threads
    )
    ch_phased_vcf = LONGPHASE_PHASE.out.phased_vcf
    ch_versions   = ch_versions.mix(LONGPHASE_PHASE.out.versions)

    // ── Phase block statistics ────────────────────────────────────────────────
    PHASE_STATS ( ch_phased_vcf )
    ch_reports  = ch_reports.mix(PHASE_STATS.out.stats.map { it -> it[1] })
    ch_versions = ch_versions.mix(PHASE_STATS.out.versions)

    // ── longphase haplotag ────────────────────────────────────────────────────
    // Write HP:i:1, HP:i:2 (or HP:i:0 for unphased) tags into the BAM.
    // Output is a haplotagged BAM, sorted and indexed.
    // This BAM is the key deliverable for:
    //   - Per-haplotype methylation (modkit --partition-tag HP)
    //   - Allele-specific methylation (ASM) analysis
    //   - Trio-aware phasing validation
    LONGPHASE_HAPLOTAG (
        ch_bam_bai.join(ch_phased_vcf, by: 0),
        ch_reference,
        ch_reference_fai
    )
    ch_haplotagged_bam = LONGPHASE_HAPLOTAG.out.bam
    ch_versions        = ch_versions.mix(LONGPHASE_HAPLOTAG.out.versions)

    SAMTOOLS_INDEX ( ch_haplotagged_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    phased_vcf      = ch_phased_vcf          // [ meta, path(phased.vcf.gz) ]
    haplotagged_bam = ch_haplotagged_bam     // [ meta, path(haplotagged.bam) ]
    haplotagged_bai = SAMTOOLS_INDEX.out.bai // [ meta, path(haplotagged.bam.bai) ]
    reports         = ch_reports
    versions        = ch_versions
}
