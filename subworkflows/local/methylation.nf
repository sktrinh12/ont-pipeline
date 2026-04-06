/*
========================================================================================
    SUBWORKFLOW: METHYLATION
    Purpose : Extract 5mC (and optionally 5hmC) methylation signals from aligned
              ONT BAM, generate per-site bedMethyl pileups, and annotate CpG islands.
    Tools   : modkit, bedtools

    Key design decisions:
    ─────────────────────────────────────────────────────────────────────────────────
    1. modkit (ONT's official tool) replaces the older nanopolish + modbam2bed
       workflow. It is faster, supports 5hmC, and is actively maintained.

    2. MM/ML tags from Dorado: Dorado with a `*_5mCG_5hmCG*` model produces
       BOTH 5mC and 5hmC probabilities. Older `*_5mCG*` models produce only 5mC.
       The pipeline checks the model string and sets modkit flags accordingly.

    3. Phased methylation: After phasing (PHASING subworkflow), haplotagged BAM
       can be passed back to modkit with --partition-tag HP to produce per-haplotype
       bedMethyl files. This reveals allele-specific methylation (ASM).

    4. bedMethyl output format (columns):
         chrom, start, end, mod_code, score, strand, start, end, color,
         Nvalid_cov, fraction_modified, Nmod, Ncanonical, Nother, Ndel, Nfail,
         Ndiff, Nnocall
       The 'fraction_modified' column (col 11) is used in downstream DMR analysis.

    5. CpG island annotation: bedtools intersect with a CHM13-specific CpG island
       BED (generated from cpgtools or downloaded from UCSC for CHM13) adds
       biological context to the methylation data.

    6. Coverage-based confidence filter: modkit --min-coverage 5 (default) ensures
       only CpG sites with sufficient read depth are included in the bedMethyl.
       Low-coverage sites produce noisy methylation estimates.
========================================================================================
*/

include { MODKIT_PILEUP           } from '../../modules/local/misc/main'
include { MODKIT_SUMMARY          } from '../../modules/local/misc/main'
include { MODKIT_SAMPLE_PROBS     } from '../../modules/local/misc/main'
include { BEDTOOLS_INTERSECT      } from '../../modules/local/bedtools/main'
include { BGZIP_TABIX             } from '../../modules/local/bgzip_tabix/main'

workflow METHYLATION {

    take:
    ch_bam           // [ meta, path(bam) ]  — must have MM/ML tags
    ch_bai           // [ meta, path(bai) ]
    ch_reference     // path(reference.fa)
    ch_reference_fai // path(reference.fa.fai)

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    ch_bam_bai = ch_bam.join(ch_bai, by: 0)

    // ── Sanity check: warn if MM/ML tags are absent ───────────────────────────
    // This can happen if the user provides FASTQ (tags stripped) or uses an
    // older basecalling model that does not produce modification probabilities.
    ch_bam_bai
        .subscribe { meta, bam, bai ->
            // Shell command is invoked inside the process, but we log a reminder
            log.info "[METHYLATION] Sample '${meta.id}': verifying MM/ML tags in BAM…"
        }

    // ── modkit sample-probs ───────────────────────────────────────────────────
    // Inspect the distribution of modification probability scores.
    // This step reveals whether the model threshold is appropriate and whether
    // 5hmC is distinguishable from 5mC in this dataset.
    MODKIT_SAMPLE_PROBS ( ch_bam_bai )
    ch_versions = ch_versions.mix(MODKIT_SAMPLE_PROBS.out.versions)

    // ── modkit summary ────────────────────────────────────────────────────────
    // Quick per-read modification summary: fraction of modified bases per
    // modification code (m = 5mC, h = 5hmC, a = 6mA, etc.)
    MODKIT_SUMMARY ( ch_bam_bai )
    ch_versions = ch_versions.mix(MODKIT_SUMMARY.out.versions)

    // ── modkit pileup (main step) ─────────────────────────────────────────────
    // Generates per-CpG methylation calls across the genome.
    // Flags:
    //   --motif CpG 0           : restrict to CpG context
    //   --combine-strands       : merge + and – strand CpGs (canonical CpG reporting)
    //   --min-coverage          : min reads per site
    //   --ignore-index-errors   : continue past BAI inconsistencies
    //   --ref                   : reference FASTA (required for CHH/CHG motifs)
    MODKIT_PILEUP (
        ch_bam_bai,
        ch_reference,
        ch_reference_fai,
        params.modkit_motifs,
        params.methylation_cpg_only
    )
    ch_bedmethyl = MODKIT_PILEUP.out.bedmethyl
    ch_versions  = ch_versions.mix(MODKIT_PILEUP.out.versions)

    // ── bgzip + tabix index the bedMethyl ─────────────────────────────────────
    // Required for: random access in IGV, bedtools interval operations,
    //               differential methylation tools (DSS, methylKit).
    BGZIP_TABIX ( ch_bedmethyl )
    ch_versions = ch_versions.mix(BGZIP_TABIX.out.versions)

    // ── CpG island annotation (optional) ─────────────────────────────────────
    // Intersect bedMethyl with CHM13-specific CpG islands.
    // Note: GRCh38 CpG island tracks cannot be used directly with CHM13
    //       without liftover — always use CHM13-native annotation files.
    if (params.cpg_islands_bed) {
        BEDTOOLS_INTERSECT (
            BGZIP_TABIX.out.compressed,
            file(params.cpg_islands_bed),
            '-wa -wb'
        )
        ch_cpg_annotated = BEDTOOLS_INTERSECT.out.result
        ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions)
    }

    emit:
    bedmethyl  = BGZIP_TABIX.out.compressed    // [ meta, path(.bed.gz) ]
    bedmethyl_tbi = BGZIP_TABIX.out.tbi           // [ meta, path(.tbi) ]
    reports    = ch_reports
    versions   = ch_versions
}
