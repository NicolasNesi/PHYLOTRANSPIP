#!/usr/bin/env nextflow

/*
 * Évaluation de complétude BUSCO sur l'assemblage brut et filtré
 */

process BUSCO {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/busco_Summaries", mode: 'copy'

    input:
    tuple val(sample_id), path(raw_assembly), path(filtered_assembly)

    output:
    path "${sample_id}_raw_short_summary.txt", emit: raw_summary
    path "${sample_id}_raw_full_table.tsv", emit: raw_table
    path "${sample_id}_raw_single_copy_busco_sequences", emit: raw_sequences
    path "${sample_id}_filtered_short_summary.txt", emit: filtered_summary
    path "${sample_id}_filtered_full_table.tsv", emit: filtered_table
    path "${sample_id}_filtered_single_copy_busco_sequences", emit: filtered_sequences
    path "${sample_id}_busco_figure.png", emit: figure

    script:
    """
    busco \\
        --in ${raw_assembly} \\
        --out BUSCO_raw_out \\
        --lineage_dataset ${params.busco_lineage} \\
        --mode transcriptome \\
        --cpu ${task.cpus} \\
        --force \\
        --metaeuk \\
        --datasets_version odb12

    busco \\
        --in ${filtered_assembly} \\
        --out BUSCO_filtered_out \\
        --lineage_dataset ${params.busco_lineage} \\
        --mode transcriptome \\
        --cpu ${task.cpus} \\
        --force \\
        --metaeuk \\
        --datasets_version odb12

    mv BUSCO_raw_out/run_${params.busco_lineage}/short_summary.txt ${sample_id}_raw_short_summary.txt
    mv BUSCO_raw_out/run_${params.busco_lineage}/full_table.tsv ${sample_id}_raw_full_table.tsv
    cp -r BUSCO_raw_out/run_${params.busco_lineage}/busco_sequences/single_copy_busco_sequences ${sample_id}_raw_single_copy_busco_sequences

    mv BUSCO_filtered_out/run_${params.busco_lineage}/short_summary.txt ${sample_id}_filtered_short_summary.txt
    mv BUSCO_filtered_out/run_${params.busco_lineage}/full_table.tsv ${sample_id}_filtered_full_table.tsv
    cp -r BUSCO_filtered_out/run_${params.busco_lineage}/busco_sequences/single_copy_busco_sequences ${sample_id}_filtered_single_copy_busco_sequences

    mkdir busco_summaries_tmp
    cp ${sample_id}_raw_short_summary.txt busco_summaries_tmp/
    cp ${sample_id}_filtered_short_summary.txt busco_summaries_tmp/
    busco --plot busco_summaries_tmp
    mv busco_summaries_tmp/busco_figure.png ${sample_id}_busco_figure.png
    """
}