#!/usr/bin/env nextflow

/*
 * Map reads against contaminant genomes
 */

process DECONTAM_MAP {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/decontam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_R1.clean.fq.gz"), path("${sample_id}_R2.clean.fq.gz"), emit: clean
    path "${sample_id}_decontam.log",       emit: log
    path "${sample_id}_ref_counts.tsv",     emit: ref_counts

    script:
    def index_basename = index_files[0].name.replaceAll(/\.\d\.bt2l?$/, '')
    """
    bowtie2 \\
        -x ${index_basename} \\
        -1 ${reads_r1} -2 ${reads_r2} \\
        -p ${task.cpus} \\
        --un-conc-gz ${sample_id}_R%.clean.fq.gz \\
        2> ${sample_id}_decontam.log \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -

    samtools index ${sample_id}.sorted.bam

    echo -e "reference\\tn_reads_mapped" > ${sample_id}_ref_counts.tsv
    samtools idxstats ${sample_id}.sorted.bam | awk -F'\\t' '\$3>0 {print \$1"\\t"\$3}' >> ${sample_id}_ref_counts.tsv

    rm ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
    """
}