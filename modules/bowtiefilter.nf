#!/usr/bin/env nextflow

/*
 * Filter rRNA + mtDNA of your taxa of interest (e.g. Phyllostomidae)
 */
 
process BOWTIEFILTER {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/bowtie", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)
    path rrna_index_files
    path mito_index_files

    output:
    tuple val(sample_id), path("${sample_id}_cleaned_edited.1.fq.gz"), path("${sample_id}_cleaned_edited.2.fq.gz"), emit: cleaned
    path "${sample_id}_rRNA.log", emit: rrna_log
    path "${sample_id}_mito.log", emit: mito_log
    path "${sample_id}_contam_summary.tsv", emit: summary

    script:
    def rrna_index = rrna_index_files[0].name.replaceAll(/\.\d\.bt2l?$/, '')
    def mito_index = mito_index_files[0].name.replaceAll(/\.\d\.bt2l?$/, '')
    """
    bowtie2 \\
        -k 1 \\
        --time \\
        --threads ${task.cpus} \\
        --al-conc matches_rRNA.fq \\
        --un-conc ${sample_id}_rRNAcleaned.fq \\
        -x ${rrna_index} \\
        -1 ${reads_r1} \\
        -2 ${reads_r2} \\
        > /dev/null 2> ${sample_id}_rRNA.log

    bowtie2 \\
        -k 1 \\
        --time \\
        --threads ${task.cpus} \\
        --al-conc matches_mito.fq \\
        --un-conc ${sample_id}_Cleaned.fq \\
        -x ${mito_index} \\
        -1 ${sample_id}_rRNAcleaned.1.fq \\
        -2 ${sample_id}_rRNAcleaned.2.fq \\
        > /dev/null 2> ${sample_id}_mito.log

    sed 's/\\/1//g' ${sample_id}_Cleaned.1.fq | gzip > ${sample_id}_cleaned_edited.1.fq.gz
    sed 's/\\/2//g' ${sample_id}_Cleaned.2.fq | gzip > ${sample_id}_cleaned_edited.2.fq.gz

    rrna_pct=\$(grep "overall alignment rate" ${sample_id}_rRNA.log | awk '{print \$1}' | tr -d '%')
    mito_pct_of_remaining=\$(grep "overall alignment rate" ${sample_id}_mito.log | awk '{print \$1}' | tr -d '%')
    mito_pct_of_total=\$(awk -v m="\$mito_pct_of_remaining" -v r="\$rrna_pct" 'BEGIN { printf "%.2f", m * (100 - r) / 100 }')

    echo -e "sample_id\\trRNA_pct\\tmtDNA_pct_of_remaining\\tmtDNA_pct_of_total" > ${sample_id}_contam_summary.tsv
    echo -e "${sample_id}\\t\${rrna_pct}\\t\${mito_pct_of_remaining}\\t\${mito_pct_of_total}" >> ${sample_id}_contam_summary.tsv
    """
}