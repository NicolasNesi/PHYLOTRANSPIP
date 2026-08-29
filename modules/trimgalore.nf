#!/usr/bin/env nextflow

/*
 * Cleaning of raw reads with Trim galore
 */
 
process TRIMGALORE {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/trimgalore", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fixed_val_1.fq.gz"), path("${sample_id}_R2.fixed_val_2.fq.gz"), emit: trimmed
    path "*_unpaired_*.fq.gz",      emit: unpaired, optional: true
    path "*trimming_report.txt",    emit: report

    script:
    """
    trim_galore \\
        --phred33 \\
        --gzip \\
        --quality 20 \\
        --length 50 \\
        --cores ${task.cpus} \\
        --retain_unpaired \\
        --paired \\
        ${reads_r1} ${reads_r2}
    """
}