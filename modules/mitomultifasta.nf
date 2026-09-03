#!/usr/bin/env nextflow

/*
 * Regroupe les gènes mitochondriaux et les génomes entiers de tous les échantillons
 * en multifasta (un par gène + un pour le génome entier)
 */
 
process MITO_MULTIFASTA {
    publishDir "${params.outdir}/mito_multifasta", mode: 'copy'

    input:
    path all_gene_fastas
    path all_whole_mito
    path gene_list

    output:
    path "*_MitoMultiFasta.fasta", emit: gene_multifasta
    path "MitoGenome_Whole_MultiFasta.fasta", emit: whole_multifasta
    path "Statistic_MitoMultiFasta.tsv", emit: stats

    script:
    """
    while read gene; do
        cat *_\${gene}.fasta >> \${gene}_MitoMultiFasta.fasta 2>/dev/null || true
    done < ${gene_list}

    find . -maxdepth 1 -name "*_MitoMultiFasta.fasta" -empty -delete

    cat ${all_whole_mito} > MitoGenome_Whole_MultiFasta.fasta

    seqkit stats *_MitoMultiFasta.fasta MitoGenome_Whole_MultiFasta.fasta -T | csvtk cut -t -f "file,num_seqs,sum_len" > Statistic_MitoMultiFasta.tsv
    """
}