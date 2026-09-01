#!/usr/bin/env nextflow

/*
 * Filtrage des CDS orthologues : header, longueur min, codons stop prématurés, sélection 1-to-1
 */
 
process FILTERFASTA {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/orthologs", mode: 'copy'

    input:
    tuple val(sample_id), path(orthologs)
    path ortholog_1to1_list

    output:
    tuple val(sample_id), path("*_CDS_OnetoOne.fasta"), emit: onetoone
    path "${sample_id}_statistic_CDSs.tsv", emit: stats

    script:
    """
    awk '{split(\$0,a,"|")} /^>/ {\$0=a[1]"|"a[2]"_"a[3]}1' ${orthologs} > temp.file
    sed 's/Hg38@//g' temp.file > ${sample_id}_multifasta_EditedHeader.fasta
    rm temp.file

    seqkit seq --min-len 150 ${sample_id}_multifasta_EditedHeader.fasta > ${sample_id}_noShort150.fasta

    premature_stops.pl -i ${sample_id}_noShort150.fasta -o ${sample_id}_noShort150_noPrematureSTOP.fasta

    while read -r line; do
        p=\$(grep -h -A 1 "\$line" ${sample_id}_noShort150_noPrematureSTOP.fasta)
        if [ -n "\$p" ]; then
            echo "\$p" > ${sample_id}_\${line}_CDS_OnetoOne.fasta
        fi
    done < ${ortholog_1to1_list}

    # Garde la séquence la plus longue en cas d'ID dupliqué, puis retire les doublons exacts
    for f in *_CDS_OnetoOne.fasta; do
        keep_longest.py "\$f" > "\${f%.fasta}_tmp.fasta"
        seqkit rmdup --by-seq -w 0 "\${f%.fasta}_tmp.fasta" -o "\$f"
        rm "\${f%.fasta}_tmp.fasta"
    done

    seqkit stats *_CDS_OnetoOne.fasta -T | csvtk cut -t -f "file,num_seqs,sum_len" > ${sample_id}_statistic_CDSs.tsv
    """
}