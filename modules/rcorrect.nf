#!/usr/bin/env nextflow

/*
 * RNA-seq error corrector
 */


process RCORRECT {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/trimgalore", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fixed.fq.gz"), path("${sample_id}_R2.fixed.fq.gz"), emit: corrected

    script:
    """
    run_rcorrector.pl \\
        -t ${task.cpus} \\
        -verbose \\
        -1 ${reads_r1} \\
        -2 ${reads_r2} > /dev/null

    conda deactivate
    conda activate python2_env

    python fixrcorrector.py -o rcorrected -1 *R1.cor.fq.gz -2 *R2.cor.fq.gz

    mv rcorrected_${sample_id}_R1.cor.fq ${sample_id}_R1.fixed.fq
    mv rcorrected_${sample_id}_R2.cor.fq ${sample_id}_R2.fixed.fq

    gzip *.fixed.fq
    rm -f *.cor.fq.gz
    """
}