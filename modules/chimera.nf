#!/usr/bin/env nextflow

/*
 * Détection de chimères par blastx (Yang & Smith 2013), sur les deux assemblages filtrés
 */
 
process CHIMERA {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/chimera", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path blast_db

    output:
    tuple val(sample_id), path("chimera_out/${sample_id}.filtered_transcripts.fa"), emit: filtered
    path "chimera_out/${sample_id}.chimera_transcripts.fa", emit: chimeras

    script:
    """
    mkdir chimera_out
    run_chimera_detection.py \\
        ${assembly} \\
        ${blast_db} \\
        ${task.cpus} \\
        ./chimera_out/

    rm -f chimera_out/Hg38_*
    """
}