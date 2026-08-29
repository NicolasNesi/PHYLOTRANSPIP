#!/usr/bin/env nextflow

/*
 * Retire les contigs Trinity qui matchent le rRNA
 */
 
process REBOWTIE {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/rebowtie", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path rrna_index_files

    output:
    tuple val(sample_id), path("${sample_id}.Trinity.Cleaned.fasta"), emit: cleaned
    path "${sample_id}_rebowtie.log", emit: log

    script:
    def rrna_index = rrna_index_files[0].name.replaceAll(/\.\d\.bt2l?$/, '')
    """
    bowtie2 \\
        -k 1 \\
        --time \\
        --threads ${task.cpus} \\
        -f \\
        --al matches_rRNA.fasta \\
        --un ${sample_id}.Trinity.Cleaned.fasta \\
        -x ${rrna_index} \\
        -U ${assembly} \\
        > /dev/null 2> ${sample_id}_rebowtie.log
    """
}