#!/usr/bin/env nextflow

/*
 * Arbre de gène avec IQ-TREE, uniquement sur les alignements avec assez d'espèces (QC)
 */
 
process GENETREE {
    tag "$gene"

    publishDir "${params.outdir}/Genetree", mode: 'copy'

    input:
    tuple val(gene), path(alignment)

    output:
    tuple val(gene), path("${gene}.treefile"), emit: tree, optional: true
    path "${gene}.iqtree", optional: true
    path "${gene}.log", optional: true

    script:
    """
    n_seq=\$(grep -c '^>' ${alignment})
    if [ "\$n_seq" -lt ${params.min_species_genetree} ]; then
        echo "Skipping ${gene}: only \$n_seq species (<${params.min_species_genetree})"
        exit 0
    fi

    iqtree \\
        -s ${alignment} \\
        -m MPF \\
        --ufboot 1000 \\
        -bnni \\
        --prefix ${gene} \\
        -T AUTO
    """
}