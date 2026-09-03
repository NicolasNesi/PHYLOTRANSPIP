#!/usr/bin/env nextflow

/*
 * Alignement MAFFT du génome mitochondrial entier, tous échantillons confondus
 */
 
process MITO_WHOLE_ALIGN {
    publishDir "${params.outdir}/mito_multifasta", mode: 'copy'

    input:
    path whole_multifasta

    output:
    path "MitoGenome_Whole_aligned.fasta", emit: aligned
    path "MitoGenome_Whole_stats.tsv", emit: stats

    script:
    """
    mafft --auto --thread ${task.cpus} ${whole_multifasta} > MitoGenome_Whole_aligned.fasta

    seqkit stats MitoGenome_Whole_aligned.fasta -T | csvtk cut -t -f "file,num_seqs,sum_len" > MitoGenome_Whole_stats.tsv
    """
}