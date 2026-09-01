#!/usr/bin/env nextflow
/*
 * QC des arbres de gènes avec TreeShrink : détecte et retire les taxons outliers (branches longues)
 */
process QCTREE {
    publishDir "${params.outdir}/TreeShrink", mode: 'copy'

    input:
    path alignments
    path trees

    output:
    path "MatrixOK/*_MatrixOK.fasta", emit: filtered_alignments
    path "indir/*/output.tree", emit: shrunk_trees
    path "List_seq_removed_TreeShrink.txt", emit: removed_seqs

    script:
    """
    mkdir -p indir MatrixOK

    for t in *.treefile; do
        gene=\$(basename "\$t" .treefile)
        mkdir -p indir/\${gene}
        cp "\$t" indir/\${gene}/input.tree
        cp "\${gene}_GuidancePRANK_algn_cleaned.fasta" indir/\${gene}/input.fasta
    done

    run_treeshrink.py -i indir -t input.tree -a input.fasta > treeshrink.log

    for t in *.treefile; do
        gene=\$(basename "\$t" .treefile)
        cp indir/\${gene}/output.fasta MatrixOK/\${gene}_MatrixOK.fasta
    done

    cat indir/*/output.txt > List_seq_removed_TreeShrink.txt
    
    # Retire |ENSG_ENST_GENE des en-têtes pour préparer la concaténation / GENETREE2
    for f in MatrixOK/*_MatrixOK.fasta; do
        sed -E 's/\\|ENSG.*//g' "\$f" > "\${f%.fasta}_ForConcat.fasta"
    done
    """
}