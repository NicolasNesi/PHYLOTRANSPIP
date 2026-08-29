#!/usr/bin/env nextflow

/*
 * Construit un multifasta par gène orthologue, à partir des CDS 1-to-1 de tous les échantillons
 */
 
process MAKEMULTIFASTA {
    publishDir "${params.outdir}/multifasta", mode: 'copy'

    input:
    path all_orthologs
    path gene_list

    output:
    path "*_MultiFasta_forGuidance_GeneName.fasta", emit: multifasta
    path "Statistic_MultiFasta.tsv", emit: stats

    script:
    """
    while read l; do
        cat *_\${l}_CDS_OnetoOne.fasta >> \${l}_MultiFasta_forGuidance.fasta 2>/dev/null || true
    done < ${gene_list}

    find . -maxdepth 1 -name "*_MultiFasta_forGuidance.fasta" -empty -delete

    for f in *_MultiFasta_forGuidance.fasta; do
        [ -e "\$f" ] || continue
        ensg_code=\$(basename "\$f" | cut -d'_' -f1)
        gene_symbol=\$(basename "\$f" | cut -d'_' -f2)
        awk -v g="\$gene_symbol" -v e="\$ensg_code" '
            /^>/ {
                split(\$0, parts, "|")
                sample = substr(parts[1], 2)
                \$0 = ">"g"_"e"_"sample
            } 1
        ' "\$f" > "\${f%.fasta}_GeneName.fasta"
    done

    seqkit stats *_MultiFasta_forGuidance_GeneName.fasta -T | csvtk cut -t -f "file,num_seqs,sum_len" > Statistic_MultiFasta.tsv
    """
}