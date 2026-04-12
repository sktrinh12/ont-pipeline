/*
=======================================================================================
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
             .min  — LR minimizer index (from vg autoindex --workflow lr-giraffe)
             .dist — LR distance index (from vg autoindex)
             .zipcodes — LR zipcode index (from vg autoindex)

    3. Read renaming (Part 3 of visualization workflow):
         - Zero-pad read names so lexicographic sort matches numeric sort
         - Format: @SAMPLENAME_BARCODE.000001
         - Critical for proper ordering in odgi visualizations

    4. Filtering (Part 5 of visualization workflow):
         - vg filter -q 1 removes MAPQ=0 (unmapped) reads
         - Unmapped reads produce empty paths that crash odgi build

    5. Augmentation (Part 6 of visualization workflow):
         - vg augment --label-paths adds read names as graph paths
         - Essential for sample reads to appear in odgi viz

    6. odgi visualisation (Part 9 of visualization workflow):
         - heatmap: compressed coverage across all paths
         - depth: per-read depth coloring (black→blue)
         - clustered: Jaccard-similarity ordering

    7. Chromosome-level parallelism:
         - vg giraffe is memory-intensive (~300 GB RAM for the full HPRC v1.1 graph).
         - For HPC runs with limited RAM, use the per-chromosome subgraph strategy:
           extract chr-specific subgraphs with `vg chunk` and map in parallel.
         - The pipeline defaults to full-graph mapping but supports chunked mode
           via params.vg_chunk_mode.
=======================================================================================
*/

include { READ_RENAME        } from '../../modules/local/vg/main'
include { VG_GIRAFFE         } from '../../modules/local/vg/main'
include { VG_FILTER_MAPQ     } from '../../modules/local/vg/main'
include { VG_SURJECT         } from '../../modules/local/vg/main'
include { VG_CALL            } from '../../modules/local/vg/main'
include { VG_STATS           } from '../../modules/local/vg/main'
include { VG_AUGMENT         } from '../../modules/local/vg/main'
include { VG_COMBINE         } from '../../modules/local/vg/main'
include { ADD_SAMPLE_PATHS   } from '../../modules/local/vg/main'
include { ODGI_BUILD         } from '../../modules/local/odgi/main'
include { ODGI_PATHS_FILTER  } from '../../modules/local/odgi/main'
include { ODGI_VIZ_HEATMAP   } from '../../modules/local/odgi/main'
include { ODGI_VIZ_DEPTH     } from '../../modules/local/odgi/main'
include { ODGI_VIZ_CLUSTERED } from '../../modules/local/odgi/main'
include { SAMTOOLS_SORT      } from '../../modules/local/samtools/main'
include { SAMTOOLS_INDEX     } from '../../modules/local/samtools/main'

workflow PANGENOME {

    take:
    ch_reads       // [ meta, path(reads.fastq.gz) ]
    ch_gbz         // path(hprc_v1.1.gbz)

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Validate LR index files - use pangenome_prefix to derive paths
    def prefix = params.pangenome_prefix ?: 'chr22_lr'
    ch_min  = params.pangenome_min_index ? file(params.pangenome_min_index) : file("${prefix}.longread.withzip.min")
    ch_dist = params.pangenome_dist ? file(params.pangenome_dist) : file("${prefix}.dist")
    ch_zip  = params.pangenome_zipcode ? file(params.pangenome_zipcode) : file("${prefix}.longread.zipcodes")
    
    if (!ch_min.exists()) {
        log.error "[PANGENOME] LR minimizer index not found: ${ch_min}. " +
                  "Run: vg autoindex --workflow lr-giraffe --prefix ${prefix} --gbz ${ch_gbz}"
        System.exit(1)
    }
    if (!ch_dist.exists()) {
        log.error "[PANGENOME] LR distance index not found: ${ch_dist}. " +
                  "Run: vg autoindex --workflow lr-giraffe --prefix ${prefix} --gbz ${ch_gbz}"
        System.exit(1)
    }
    if (!ch_zip.exists()) {
        log.warn "[PANGENOME] LR zipcode index not found: ${ch_zip}. Proceeding without zipcodes."
        ch_zip = null
    }

    // ── Part 3: Read renaming ─────────────────────────────────────────────────────
    // Zero-pad read names so lexicographic sort matches numeric sort
    // Format: @SAMPLEID_BARCODE.000001
    READ_RENAME( ch_reads )
    ch_renamed_reads = READ_RENAME.out.fastq
    ch_versions = ch_versions.mix(READ_RENAME.out.versions)

    // ── Part 4: vg giraffe: map to pangenome graph ────────────────────────────────
    // Outputs: GAF (Graph Alignment Format) or GAM (Graph Alignment/Map binary)
    // We use GAM for downstream vg call; GAF is human-readable for debugging.
    //
    // Key flags:
    //   -b r10   : ONT R10.4.1 (use -b hifi for PacBio HiFi)
    //   -Z  : GBZ graph input
    //   -m  : LR minimizer index
    //   -d  : LR distance index
    //   -z  : LR zipcode index
    //   -f  : input FASTQ
    //   -t  : threads
    //   -p  : progressive output
    VG_GIRAFFE (
        ch_renamed_reads,
        ch_gbz,
        ch_min,
        ch_dist,
        ch_zip,
        params.vg_threads
    )
    ch_gam      = VG_GIRAFFE.out.gam
    ch_versions = ch_versions.mix(VG_GIRAFFE.out.versions)

    // ── Part 5: vg filter - Remove unmapped reads (MAPQ=0) ────────────────────────
    // Unmapped reads produce empty paths that crash odgi build
    VG_FILTER_MAPQ( ch_gam )
    ch_filtered_gam = VG_FILTER_MAPQ.out.gam
    ch_versions     = ch_versions.mix(VG_FILTER_MAPQ.out.versions)

    // ── vg stats: mapping quality report ─────────────────────────────────────────
    // Reports: total reads, mapped %, perfect mappings, mapping quality distribution
    VG_STATS ( ch_filtered_gam )
    ch_reports  = ch_reports.mix(VG_STATS.out.stats.map { it[1] })
    ch_versions = ch_versions.mix(VG_STATS.out.versions)

    // ── vg surject: project graph alignments → linear BAM ────────────────────────
    // Surjection maps each graph alignment back to the closest reference path
    // (CHM13 = CHM13v2.0). This produces a standard BAM against the linear reference,
    // enabling comparison with minimap2 alignments and use with linear-reference tools.
    //
    // The surjected BAM is complementary to the minimap2 BAM:
    //   - minimap2 BAM: better for regions well-represented in CHM13
    //   - surjected BAM: better for regions with high haplotype diversity
    VG_SURJECT (
        ch_filtered_gam,
        ch_gbz,
        params.vg_threads
    )
    ch_surject_bam = VG_SURJECT.out.bam
    ch_versions    = ch_versions.mix(VG_SURJECT.out.versions)

    SAMTOOLS_SORT  ( ch_surject_bam )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

    // ── vg call: variant calling on the graph ──────────────────────────────────────
    // Detects variants that are encoded as alternate paths in the pangenome graph.
    // These include population-specific alleles absent from the CHM13 linear ref.
    VG_CALL (
        ch_filtered_gam,
        ch_gbz,
        'CHM13'
    )
    ch_versions = ch_versions.mix(VG_CALL.out.versions)

    // ── Part 6-9: vg augment + odgi viz: consolidated family visualization ───────
    // Instead of per-sample processing, we:
    // 1. Augment each sample individually to get VGs with sample paths (--label-paths)
    // 2. Combine all augmented VGs into one graph with VG_COMBINE
    // 3. Build one ODGI from the combined graph
    // 4. Generate 3 visualizations: heatmap, depth, clustered
    //
    // This produces consolidated PNGs instead of one per sample.
    if (params.odgi_region) {
        // Augment each sample individually (with --label-paths to preserve read paths)
        VG_AUGMENT(
            ch_filtered_gam,
            ch_gbz
        )
        ch_single_augmented = VG_AUGMENT.out.vg
        ch_versions         = ch_versions.mix(VG_AUGMENT.out.versions)

        // Collect all augmented VGs into a single combined graph
        // Group all augmented VGs together and create family meta
        ch_augmented_list = ch_single_augmented
            .map { [ [id: 'family', family_id: it[0].family_id], it[1] ] }
            .collect()

        // Combine all augmented VGs into one
        VG_COMBINE(
            ch_augmented_list
        )
        ch_combined_vg = VG_COMBINE.out.vg
        ch_versions    = ch_versions.mix(VG_COMBINE.out.versions)

        // Extract sample IDs from the input reads channel
        // Format: sampleid_barcode for matching renamed reads
        ch_sample_ids = ch_reads.map { "${it[0].id}_${it[0].barcode}" }.collect()

        // Build ODGI from combined graph (sample read paths already added via --label-paths in VG_AUGMENT)
        ODGI_BUILD(
            ch_combined_vg
        )
        ch_odgi = ODGI_BUILD.out.odgi
        ch_versions = ch_versions.mix(ODGI_BUILD.out.versions)

        // ── Part 8: Prepare display path list ─────────────────────────────────────
        // Restrict visualizations to CHM13 reference + sample reads
        ODGI_PATHS_FILTER(
            ch_odgi,
            ch_sample_ids
        )
        ch_path_list = ODGI_PATHS_FILTER.out.path_list
        ch_versions  = ch_versions.mix(ODGI_PATHS_FILTER.out.versions)

        // ── Part 9: Generate visualizations ────────────────────────────────────────
        
        // 9a: Compressed Coverage Heatmap
        ODGI_VIZ_HEATMAP(
            ch_odgi,
            params.odgi_region,
            params.odgi_ignore_prefix
        )
        ch_reports  = ch_reports.mix(ODGI_VIZ_HEATMAP.out.heatmap_png)
        ch_versions = ch_versions.mix(ODGI_VIZ_HEATMAP.out.versions)

        // 9b: Per-Read Depth Coloring
        ODGI_VIZ_DEPTH(
            ch_odgi,
            ch_path_list,
            params.odgi_region
        )
        ch_reports  = ch_reports.mix(ODGI_VIZ_DEPTH.out.depth_png)
        ch_versions = ch_versions.mix(ODGI_VIZ_DEPTH.out.versions)

        // 9c: Clustered View (Jaccard similarity ordering)
        ODGI_VIZ_CLUSTERED(
            ch_odgi,
            params.odgi_region,
            params.odgi_ignore_prefix
        )
        ch_reports  = ch_reports.mix(ODGI_VIZ_CLUSTERED.out.clustered_png)
        ch_versions = ch_versions.mix(ODGI_VIZ_CLUSTERED.out.versions)
    }

    emit:
    surjected_bam = SAMTOOLS_SORT.out.bam          // [ meta, path(surjected.bam) ]
    surjected_bai = SAMTOOLS_INDEX.out.bai
    graph_vcf     = VG_CALL.out.vcf                // [ meta, path(graph.vcf.gz) ]
    reports       = ch_reports
    versions      = ch_versions
}