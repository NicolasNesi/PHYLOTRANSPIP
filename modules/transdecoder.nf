#!/usr/bin/env nextflow

/*
 * Prédiction de CDS avec TransDecoder (LongOrfs + blastp + Predict)
 */

process TRANSDECODER {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/transdecoder", mode: 'copy'

    input:
    tuple val(sample_id), path(transcripts)
    path blastp_db_files

    output:
    tuple val(sample_id), path("${sample_id}.pep.fa"), emit: proteins
    tuple val(sample_id), path("${sample_id}.cds.fa"), emit: cds
    path "${sample_id}.blastp.outfmt6.gz", emit: blastp_out
    path "${sample_id}.log", emit: log

    script:
    def blastp_db = blastp_db_files[0].name.replaceAll(/\.(phr|pin|psq)$/, '')
    """
    transdecoder_wrapper.py ${transcripts} ${task.cpus} stranded ./ ${blastp_db}
    """
}