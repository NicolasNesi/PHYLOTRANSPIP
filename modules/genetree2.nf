#!/usr/bin/env nextflow

/*
 * Arbre de gène (2e passe) avec IQ-TREE sur les matrices post-TreeShrink, header édité
 */
 
process GENETREE2 {
    tag "$gene"

    publishDir "${params.outdir}/Genetree2", mode: 'copy'

    input:
    tuple val(gene), path(matrix)

    output:
    tuple val(gene), path("${gene}.treefile"), emit: tree, optional: true
    path "${gene}.iqtree", optional: true
    path "${gene}.log", optional: true

    script:
    """
    iqtree2 \\
        -s ${matrix} \\
        -m MFP \\
        --ufboot 1000 \\
        --abayes \\
        --alrt 1000 \\
        -bnni \\
        --prefix ${gene} \\
        -T AUTO
    """
}