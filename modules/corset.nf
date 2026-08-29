#!/usr/bin/env nextflow

/*
 * Clustering de transcrits avec Salmon + Corset, puis sélection du transcrit représentatif par cluster
 */

process CORSET {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/salmon_corset", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(reads_r1), path(reads_r2)

    output:
    tuple val(sample_id), path("${sample_id}.largest_cluster_transcripts.fa"), emit: representative
    path "${sample_id}.redundant_cluster_transcripts.fa", emit: redundant
    path "${sample_id}.largest_cluster.csv",   emit: largest_csv
    path "${sample_id}.redundant_cluster.csv", emit: redundant_csv
    path "${sample_id}_salmon-clusters.txt",   emit: clusters

    script:
    """
    salmon index \\
        --transcripts ${assembly} \\
        --index ${sample_id}_salmon_index \\
        --threads ${task.cpus}

    salmon quant \\
        --index ${sample_id}_salmon_index \\
        --libType A \\
        --dumpEq \\
        --hardFilter \\
        --skipQuant \\
        --threads ${task.cpus} \\
        --mates1 ${reads_r1} \\
        --mates2 ${reads_r2} \\
        --output ${sample_id}_salmon_quant

    gunzip ${sample_id}_salmon_quant/aux_info/eq_classes.txt.gz

    corset \\
        -i salmon_eq_classes ${sample_id}_salmon_quant/aux_info/eq_classes.txt \\
        -m 5 \\
        -p ${sample_id}_salmon

    conda deactivate
    conda activate python2_env

    filter_corset_output.py ${assembly} ${sample_id}_salmon-clusters.txt ./

    rm -rf ${sample_id}_salmon_quant ${sample_id}_salmon_index
    """
}