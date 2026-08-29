#!/usr/bin/env nextflow

/*
 * Retrait des séquences nucléotidiques redondantes (cd-hit-est)
 */
 
process CDHIT {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/cd-hit-est", mode: 'copy'

    input:
    tuple val(sample_id), path(cds)

    output:
    tuple val(sample_id), path("${sample_id}_cdhitest.fa"), emit: unique
    path "${sample_id}_cdhitest.fa.clstr", emit: clusters
	path "${sample_id}_unique_transcripts.tsv", emit: summary

    script:
    """
    cd-hit-est \\
        -c 0.95 \\
        -n 10 \\
        -r 0 \\
        -M ${task.memory.toMega()} \\
        -T ${task.cpus} \\
        -i ${cds} \\
        -o ${sample_id}_cdhitest.fa
		
    n_input=\$(grep -c '^>' ${cds})
    n_unique=\$(grep -c '^>' ${sample_id}_cdhitest.fa)

    echo -e "sample_id\\tn_input_transcripts\\tn_unique_transcripts" > ${sample_id}_unique_transcripts.tsv
    echo -e "${sample_id}\\t\${n_input}\\t\${n_unique}" >> ${sample_id}_unique_transcripts.tsv
    """
}