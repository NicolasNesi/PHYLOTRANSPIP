#!/usr/bin/env nextflow

/*
 * Évaluation qualité de l'assemblage avec Transrate + filtrage des transcrits mal supportés
 */
 
process TRANSRATE {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/transrate", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(reads_r1), path(reads_r2)

    output:
    tuple val(sample_id), path("${sample_id}.good_transcripts.short_name.fa"), emit: filtered
    path "${sample_id}.good_transcripts.csv", emit: good_csv
    path "${sample_id}.bad_transcripts.csv",  emit: bad_csv
    path "${sample_id}_all_contigs.csv",      emit: contigs_csv

    script:
    """
    transrate \\
        --assembly ${assembly} \\
        --left ${reads_r1} \\
        --right ${reads_r2} \\
        --threads ${task.cpus} \\
        --output transrate_out

    conda deactivate
    conda activate python2_env

    echo "Start filtering bad transcripts"
    filter_transcripts_transrate.py ${assembly} transrate_out/*/contigs.csv ./

    mv transrate_out/*/contigs.csv ${sample_id}_all_contigs.csv
    """
}