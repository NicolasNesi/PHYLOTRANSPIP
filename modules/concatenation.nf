#!/usr/bin/env nextflow

/*
 * Concaténation des alignements filtrés (TreeShrink) en une matrice unique + fichiers de partition
 */
 
process CONCATENATION {
    publishDir "${params.outdir}/MatrixOK", mode: 'copy'

    input:
    path matrix_ok_alignments

    output:
    path "*_MatrixOK*.svg"
    path "Concatenated_Matrix.fasta", emit: concat_matrix
    path "Concatenated_Partition_byCodonPos.nex", emit: partition_nexus
    path "Partition_byGene.raxml", emit: partition_raxml_gene
    path "Partition_byGeneCodon.raxml", emit: partition_raxml_genecodon
    path "Statistics_MatrixOK.tsv", emit: stats
    path "Statistic_Concatenated_Matrix.tsv", emit: concat_stats
    path "Statistic_Concatenated_Matrix_PerTaxa.tsv", emit: concat_stats_taxa

    script:
    """
    for f in *_MatrixOK.fasta; do
        gene=\$(basename "\$f" _MatrixOK.fasta)
        CIAlign --infile "\$f" --outfile "\${gene}_MatrixOK" --keep_gaponly --plot_format svg --plot_output
    done
    rm -f *_log.txt

    seqkit stat --all *_MatrixOK.fasta -T > Statistics_MatrixOK.tsv

    AMAS.py concat \\
        --in-files *_ForConcat.fasta \\
        --data-type dna \\
        --in-format fasta \\
        --out-format fasta \\
        --cores ${task.cpus} \\
        --check-align \\
        --concat-out Concatenated_Matrix.fasta \\
        --concat-part Concatenated_Partition_byCodonPos.nex \\
        --part-format nexus \\
        --codons 123

    AMAS.py summary \\
        --in-files Concatenated_Matrix.fasta \\
        --data-type dna \\
        --in-format fasta \\
        --cores ${task.cpus} \\
        --by-taxon \\
        --summary-out Statistic_Concatenated_Matrix.txt

    mv Statistic_Concatenated_Matrix.txt Statistic_Concatenated_Matrix.tsv
    
	mv Concatenated_Matrix.fasta-seq-summary.txt Statistic_Concatenated_Matrix_PerTaxa.tsv

    AMAS.py concat \\
        --in-files *_ForConcat.fasta \\
        --data-type dna \\
        --in-format fasta \\
        --out-format fasta \\
        --cores ${task.cpus} \\
        --check-align \\
        --concat-out temp_concat1.fasta \\
        --concat-part Partition_byGene.raxml \\
        --part-format raxml
    rm -f temp_concat1.fasta

    AMAS.py concat \\
        --in-files *_ForConcat.fasta \\
        --data-type dna \\
        --in-format fasta \\
        --out-format fasta \\
        --cores ${task.cpus} \\
        --check-align \\
        --concat-out temp_concat2.fasta \\
        --concat-part Partition_byGeneCodon.raxml \\
        --part-format raxml \\
        --codons 123
    rm -f temp_concat2.fasta
    """
}