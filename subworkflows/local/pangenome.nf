/*
========================================================================================
    SUBWORKFLOW: PANGENOME
    Purpose : Map reads to the HPRC v1.1 pangenome variation graph using vg giraffe,
              call variants directly on the graph, and visualise graph topology.
    Tools   : vg giraffe, vg call, vg stats, odgi (viz)

    Key design decisions:
    ─────────────────────────────────────────────────────────────────────────────────
    1. Why pangenome in addition to linear alignment?
         - The T2T-CHM13v2.0 linear reference, while complete, represents a single
           haplotype. Reads from individuals whose haplotypes diverge significantly
           from CHM13 (e.g., African populations) are forced to align to a foreign
           sequence, introducing systematic reference bias.
         - The HPRC v1.1 pangenome graph encodes 94 haplotypes from 47 diverse
           individuals, covering ~10% additional genetic variation vs. GRCh38.
         - Mapping to the graph first, then projecting onto the linear reference,
           reduces bias in regions of high structural diversity (e.g., MHC, centromeres).

    2. vg giraffe vs. vg map:
         - vg map: classic full-alignment to graph, accurate but slow.
         - vg giraffe: haplotype-aware seed-extension mapper, 5–10× faster.
           Recommended for all production pangenome mapping.
         - giraffe requires three pre-built index files from the GBZ graph:
             .gbz  — GBZ-format graph (sequence + topology + haplotypes)
             .min  — minimizer index (k-mer seeds)
             .dist — distance index (for optimal alignment scoring)
           These are available for HPRC v1.1 from the HPRC S3 bucket:
             s3://human-pangenomics/pangenomes/freeze/freeze1/graph/

    3. Graph-based variant calling:
         - vg call operates directly on the graph paths to produce a VCF against
           a chosen reference path (hs1 = CHM13v2.0 in the HPRC graph).
         - These VCFs capture pangenome-specific variants — alleles that are absent
           from the CHM13 linear reference but present in the HPRC graph.
         - Complement rather than replace Clair3/Sniffles2 results.

    4. odgi visualisation:
         - odgi viz generates a 1D linearized heat map of the graph showing path
           coverage, useful for QC and publication figures.
         - Output is a PNG/GIF, not a publication-quality figure, but it confirms
           that mapping succeeded and coverage is uniform.

    5. Chromosome-level parallelism:
         - vg giraffe is memory-intensive (~300 GB RAM for the full HPRC v1.1 graph).
         - For HPC runs with limited RAM, use the per-chromosome subgraph strategy:
           extract chr-specific subgraphs with `vg chunk` and map in parallel.
         - The pipeline defaults to full-graph mapping but supports chunked mode
           via params.vg_chunk_mode.
========================================================================================
*/

include { VG_GIRAFFE      } from '../../modules/local/vg/main'
include { VG_SURJECT       } from '../../modules/local/vg/main'
include { VG_CALL          } from '../../modules/local/vg/main'
include { VG_STATS         } from '../../modules/local/vg/main'
include { ODGI_BUILD       } from '../../modules/local/odgi/build/main'
include { ODGI_VIZ         } from '../../modules/local/odgi/viz/main'
include { SAMTOOLS_SORT    } from '../../modules/local/samtools/main'
include { SAMTOOLS_INDEX   } from '../../modules/local/samtools/main'

workflow PANGENOME {

    take:
    ch_reads       // [ meta, path(reads.fastq.gz) ]
    ch_gbz         // path(hprc_v1.1.gbz)

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Validate index files - use separate params if provided, otherwise derive from gbz
    ch_min  = params.pangenome_min_index ? file(params.pangenome_min_index) : file("${ch_gbz}.min")
    ch_dist = params.pangenome_dist ? file(params.pangenome_dist) : file("${ch_gbz}.dist")
    ch_zip  = params.pangenome_zipcode ? file(params.pangenome_zipcode) : null
    if (!ch_min.exists()) {
        log.error "[PANGENOME] Minimizer index not found: ${ch_min}. " +
                  "Run: vg minimizer -p -o graph.min graph.gbz"
        System.exit(1)
    }
    if (!ch_dist.exists()) {
        log.error "[PANGENOME] Distance index not found: ${ch_dist}. " +
                  "Run: vg index -j graph.dist graph.gbz"
        System.exit(1)
    }

    // ── vg giraffe: map to pangenome graph ────────────────────────────────────
    // Outputs: GAF (Graph Alignment Format) or GAM (Graph Alignment/Map binary)
    // We use GAM for downstream vg call; GAF is human-readable for debugging.
    //
    // Key flags:
    //   -Z  : GBZ graph input
    //   -m  : minimizer index
    //   -d  : distance index
    //   -x  : also output as BAM via surjection (--output-format BAM)
    //   -f  : input FASTQ
    //   -t  : threads
    //   --zipcode-name: use zipcode index for named-coordinate reporting
    VG_GIRAFFE (
        ch_reads,
        ch_gbz,
        ch_min,
        ch_dist,
        ch_zip,
        params.vg_threads
    )
    ch_gam      = VG_GIRAFFE.out.gam
    ch_versions = ch_versions.mix(VG_GIRAFFE.out.versions)

    // ── vg stats: mapping quality report ─────────────────────────────────────
    // Reports: total reads, mapped %, perfect mappings, mapping quality distribution
    VG_STATS ( ch_gam )
    ch_reports  = ch_reports.mix(VG_STATS.out.stats.map { it[1] })
    ch_versions = ch_versions.mix(VG_STATS.out.versions)

    // ── vg surject: project graph alignments → linear BAM ────────────────────
    // Surjection maps each graph alignment back to the closest reference path
    // (CHM13 = CHM13v2.0). This produces a standard BAM against the linear reference,
    // enabling comparison with minimap2 alignments and use with linear-reference tools.
    //
    // The surjected BAM is complementary to the minimap2 BAM:
    //   - minimap2 BAM: better for regions well-represented in CHM13
    //   - surjected BAM: better for regions with high haplotype diversity
    VG_SURJECT (
        ch_gam,
        ch_gbz,
        params.vg_threads
    )
    ch_surject_bam = VG_SURJECT.out.bam
    ch_versions    = ch_versions.mix(VG_SURJECT.out.versions)

    SAMTOOLS_SORT  ( ch_surject_bam )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

    // ── vg call: variant calling on the graph ─────────────────────────────────
    // Detects variants that are encoded as alternate paths in the pangenome graph.
    // These include population-specific alleles absent from the CHM13 linear ref.
    VG_CALL (
        ch_gam,
        ch_gbz,
        'CHM13'
    )
    ch_versions = ch_versions.mix(VG_CALL.out.versions)

    // ── odgi viz: 1D pangenome graph visualisation ────────────────────────────
    // Requires ODGI format graph built from GBZ.
    // odgi viz produces a PNG strip showing node coverage depth across the graph.
    // This is used to confirm mapping uniformity and identify under-mapped regions.
    //
    // Note: odgi viz on the full genome graph is very large. We restrict to
    // Chr22 in the test profile, or a user-specified region via params.odgi_region.
    // Region format: PATH:start-end (e.g., "chr22:0-51324926")
    //
    // TODO: Enable alignment overlay (VG issue #4240)
    // To re-enable: (1) Add GAM to ODGI_BUILD input, (2) Uncomment vg augment step
    // in odgi/build/main.nf, (3) Add -A ${prefix} -S flags to odgi viz command
    if (params.odgi_region) {
        ODGI_BUILD(
            ch_gam.map { meta, _gam -> [meta, ch_gbz] }
        )
        ch_versions = ch_versions.mix(ODGI_BUILD.out.versions)

        ODGI_VIZ(
            ODGI_BUILD.out.odgi,
            params.odgi_region
        )
        ch_reports  = ch_reports.mix(ODGI_VIZ.out.png.map { it[1] })
        ch_versions = ch_versions.mix(ODGI_VIZ.out.versions)
    }

    emit:
    surjected_bam = SAMTOOLS_SORT.out.bam          // [ meta, path(surjected.bam) ]
    surjected_bai = SAMTOOLS_INDEX.out.bai
    graph_vcf     = VG_CALL.out.vcf                // [ meta, path(graph.vcf.gz) ]
    reports       = ch_reports
    versions      = ch_versions
}
