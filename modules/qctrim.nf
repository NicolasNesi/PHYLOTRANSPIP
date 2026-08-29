#!/usr/bin/env nextflow
/*
 * FastQC on cleaned and decomtamined reads
 */
process QCTRIM {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/fastqc/${stage}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)
    val stage

    output:
	tuple val(sample_id), path("*.zip"), emit: zip
    path "*.html", emit: html

    script:
    """
    fastqc \\
        -o ./ \\
        -t ${task.cpus} \\
        --noextract \\
        --nogroup \\
        ${reads_r1} ${reads_r2}
    """
}
 