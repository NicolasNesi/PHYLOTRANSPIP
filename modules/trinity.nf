#!/usr/bin/env nextflow

/*
 * Assemblage de novo avec Trinity
 */
 
process TRINITY {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/trinity", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)

    output:
    tuple val(sample_id), path("${sample_id}.Trinity.fasta"), emit: assembly

    script:
    """
    ulimit -s unlimited

    Trinity \\
        --max_memory ${task.memory.toGiga()}G \\
        --CPU ${task.cpus} \\
        --output trinity_out \\
        --seqType fq \\
        --SS_lib_type FR \\
        --full_cleanup \\
        --no_normalize_reads \\
        --left ${reads_r1} \\
        --right ${reads_r2} \\
        --path_reinforcement_distance 75 \\
        --bflyHeapSpaceMax 100G \\
        --bflyCalculateCPU

    mv trinity_out.Trinity.fasta ${sample_id}.Trinity.fasta
    """
}