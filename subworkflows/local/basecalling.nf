/*
========================================================================================
    SUBWORKFLOW: BASECALLING
    Purpose : Convert raw pod5/fast5 signal → basecalled FASTQ with MM/ML methylation
              tags embedded (required for modkit downstream).
    Tools   : Dorado (simplex + duplex)

    Key design decisions:
    ─────────────────────────────────────────────────────────────────────────────────
    1. Dorado SUP model is used by default for highest accuracy (~99.5% per-read).
       FAST/HAC are available for speed-vs-accuracy tradeoffs.

    2. --emit-moves is set for duplex pairing (required by dorado duplex).

    3. Dorado natively outputs unmapped BAM (uBAM) with MM/ML methylation tags
       embedded. We convert to FASTQ only when methylation is DISABLED, because
       keeping the BAM format preserves the epigenetic tags needed by modkit.

    4. If duplex mode is enabled, we run:
         a) dorado duplex → duplex BAM (higher-accuracy duplex reads)
         b) dorado basecaller → simplex BAM (all reads)
         c) merge + deduplicate (duplex reads supersede their simplex counterparts)

    5. fast5 inputs are supported but the pipeline emits a deprecation warning.
       ONT has officially deprecated fast5; all new R10.4.1+ flow cells use pod5.
========================================================================================
*/

include { DORADO_BASECALLER  } from '../../modules/local/misc/main'
include { DORADO_DUPLEX      } from '../../modules/local/misc/main'
include { DORADO_MERGE       } from '../../modules/local/misc/main'
include { SAMTOOLS_BAM2FASTQ } from '../../modules/local/samtools/main'

workflow BASECALLING {

    take:
    ch_raw_reads   // [ meta, path(pod5 | fast5) ]

    main:
    ch_versions = Channel.empty()

    // ── Deprecation warning for fast5 ─────────────────────────────────────────
    ch_raw_reads
        .filter { meta, f -> meta.input_type == 'fast5' }
        .subscribe { meta, f ->
            log.warn "[BASECALLING] Sample '${meta.id}': fast5 format is deprecated by ONT. " +
                     "Convert to pod5 using 'pod5 convert fast5' for future runs."
        }

    // ── Simplex basecalling ───────────────────────────────────────────────────
    // Always performed. Produces uBAM with:
    //   - Per-read quality scores
    //   - MM/ML tags (methylation) when a model that supports modbase calling is used
    //   - MV (move tables) tags if --emit-moves is set (required for duplex)
    DORADO_BASECALLER (
        ch_raw_reads,
        params.dorado_model,
        params.dorado_min_qscore
    )
    ch_simplex_ubam = DORADO_BASECALLER.out.ubam
    ch_versions     = ch_versions.mix(DORADO_BASECALLER.out.versions)

    // ── Duplex basecalling (optional) ─────────────────────────────────────────
    // Dorado duplex pairs template + complement strands to produce reads with
    // Q50+ accuracy — roughly 2–3× more accurate than simplex SUP.
    // Requires: R10.4.1 chemistry (not compatible with R9.4.1).
    // Duplex yield is typically 10–30% of total reads (chemistry + flow cell dependent).
    if (params.duplex) {
        DORADO_DUPLEX (
            ch_raw_reads,
            params.dorado_model
        )
        ch_duplex_ubam = DORADO_DUPLEX.out.ubam
        ch_versions    = ch_versions.mix(DORADO_DUPLEX.out.versions)

        // Merge simplex + duplex; duplex reads are tagged DX:i:1 in the BAM header.
        // Simplex reads that formed a duplex pair (DX:i:-1) are removed to avoid
        // double-counting during alignment and variant calling.
        DORADO_MERGE (
            ch_simplex_ubam.join(ch_duplex_ubam, by: 0)
        )
        ch_final_ubam = DORADO_MERGE.out.ubam
        ch_versions   = ch_versions.mix(DORADO_MERGE.out.versions)
    } else {
        ch_final_ubam = ch_simplex_ubam
    }

    // ── Convert uBAM → FASTQ for aligners that need FASTQ input ──────────────
    // vg giraffe requires FASTQ; minimap2 accepts uBAM directly.
    // For methylation-aware runs: keep uBAM as the primary output (MM/ML intact).
    // For pangenome mapping: also produce FASTQ.
    SAMTOOLS_BAM2FASTQ ( ch_final_ubam )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FASTQ.out.versions)

    emit:
    reads   = SAMTOOLS_BAM2FASTQ.out.reads   // [ meta, path(fastq.gz) ]  — for vg
    ubam    = ch_final_ubam                  // [ meta, path(ubam) ]      — for minimap2 + modkit
    versions = ch_versions
}
