#!/usr/bin/env nextflow

/*
 * Get SRA RNA reads
 */


process WGETSRA {
    tag "$sample_id"
	
	publishDir "${params.outdir}/CombinedSRR/${sample_id}/sra", mode: 'copy'

    input:
    tuple val(sample_id), val(srr)

    output:
    tuple val(sample_id), path("${sample_id}_${srr}_R1.fastq"), path("${sample_id}_${srr}_R2.fastq"), emit: reads

    script:
    """
    # récupère les URLs SRA et télécharge
    wget \$(ffq --ncbi ${srr} | jq -r '.[] | .url' | tr '\\n' ' ')

    # dump SRA -> fastq
    fasterq-dump \\
        --split-files \\
        --threads ${task.cpus} \\
        --outdir ./ \\
        --seq-defline '@\$sn/\$ri' \\
        --outfile ${srr} \\
        --progress \\
        ${srr}*lite.1

    mv ${srr}_1.fastq ${sample_id}_${srr}_R1.fastq
    mv ${srr}_2.fastq ${sample_id}_${srr}_R2.fastq
    """
}