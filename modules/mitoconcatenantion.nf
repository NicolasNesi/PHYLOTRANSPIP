#!/usr/bin/env nextflow

/*
 * Concaténation des alignements mitochondriaux :
 * - matrice complète (PCG+tRNA+rRNA), partition par gène et par gène x codon
 * - matrice PCG seule, partition par gène x codon
 * - versions positions 1+2 uniquement (PCG et complète)
 */
 
process MITO_CONCAT {
    publishDir "${params.outdir}/mito_multifasta/Concatenated", mode: 'copy'

    input:
    path all_cleaned_alignments
    path all_pos12_alignments
    path pcg_gene_list

    output:
    path "MitoGenes_Concatenated.fasta", emit: full_matrix
    path "Partition_byGene.raxml", emit: partition_full_bygene
    path "MitoPCG_Concatenated.fasta", emit: pcg_matrix
    path "Partition_byGeneCodon.raxml", emit: partition_pcg_codon
    path "MitoPCG_pos12_Concatenated.fasta", emit: pcg_pos12
    path "MitoGenes_pos12_Concatenated.fasta", emit: full_pos12
    path "Statistics_MitoAlignments.tsv", emit: stats

    script:
    """
    seqkit stat --all *_GuidancePRANK_algn_cleaned.fasta -T > Statistics_MitoAlignments.tsv

    # -- Matrice complète (PCG+tRNA+rRNA), partition par gène --
    AMAS.py concat \\
        --in-files *_GuidancePRANK_algn_cleaned.fasta \\
        --data-type dna \\
        --in-format fasta \\
        --out-format fasta \\
        --cores ${task.cpus} \\
        --check-align \\
        --concat-out MitoGenes_Concatenated.fasta \\
        --concat-part Partition_byGene.raxml \\
        --part-format raxml

    # -- Sélectionne les fichiers PCG parmi les alignements nettoyés --
    pcg_files=""
    while read gene; do
        f="\${gene}_GuidancePRANK_algn_cleaned.fasta"
        [ -f "\$f" ] && pcg_files="\$pcg_files \$f"
    done < ${pcg_gene_list}

    # -- Matrice PCG seule, partition par gène x codon --
    AMAS.py concat \\
        --in-files \$pcg_files \\
        --data-type dna \\
        --in-format fasta \\
        --out-format fasta \\
        --cores ${task.cpus} \\
        --check-align \\
        --concat-out MitoPCG_Concatenated.fasta \\
        --concat-part Partition_byGeneCodon.raxml \\
        --part-format raxml \\
        --codons 123

    # -- Matrice PCG positions 1+2 (fichiers déjà produits par FILTERALIGNMENTS -m c3) --
    pcg_pos12_files=""
    while read gene; do
        f="\${gene}_GuidancePRANK_algn_pos12.fasta"
        [ -f "\$f" ] && pcg_pos12_files="\$pcg_pos12_files \$f"
    done < ${pcg_gene_list}

    AMAS.py concat \\
        --in-files \$pcg_pos12_files \\
        --data-type dna \\
        --in-format fasta \\
        --out-format fasta \\
        --cores ${task.cpus} \\
        --check-align \\
        --concat-out MitoPCG_pos12_Concatenated.fasta \\
        --concat-part Partition_byGene_PCGpos12.raxml \\
        --part-format raxml

    # -- Matrice complète positions 1+2 : PCG en pos12 + tRNA/rRNA en entier --
    non_pcg_files=""
    for f in *_GuidancePRANK_algn_cleaned.fasta; do
        base=\$(basename "\$f" _GuidancePRANK_algn_cleaned.fasta)
        if ! grep -qx "\$base" ${pcg_gene_list}; then
            non_pcg_files="\$non_pcg_files \$f"
        fi
    done

    AMAS.py concat \\
        --in-files \$pcg_pos12_files \$non_pcg_files \\
        --data-type dna \\
        --in-format fasta \\
        --out-format fasta \\
        --cores ${task.cpus} \\
        --check-align \\
        --concat-out MitoGenes_pos12_Concatenated.fasta \\
        --concat-part Partition_byGene_Fullpos12.raxml \\
        --part-format raxml
    """
}