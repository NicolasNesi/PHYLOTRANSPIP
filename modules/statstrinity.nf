#!/usr/bin/env nextflow

/*
 * Statistiques d'assemblage Trinity
 */
 
process STATSTRINITY {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/trinity", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}_Assembly-Stats.txt", emit: full_stats
    path "${sample_id}_stats_summary.tsv",  emit: summary

    script:
    """
    perl TrinityStats.pl ${assembly} > ${sample_id}_Assembly-Stats.txt

    total_genes=\$(grep "Total trinity 'genes'" ${sample_id}_Assembly-Stats.txt | cut -d: -f2 | tr -d ' ')
    total_transcripts=\$(grep "Total trinity transcripts" ${sample_id}_Assembly-Stats.txt | cut -d: -f2 | tr -d ' ')
    percent_gc=\$(grep "Percent GC" ${sample_id}_Assembly-Stats.txt | cut -d: -f2 | tr -d ' ')
    contig_n50=\$(grep "Contig N50" ${sample_id}_Assembly-Stats.txt | cut -d: -f2 | tr -d ' ')
    median_len=\$(grep "Median contig length" ${sample_id}_Assembly-Stats.txt | cut -d: -f2 | tr -d ' ')
    average_len=\$(grep "Average contig" ${sample_id}_Assembly-Stats.txt | cut -d: -f2 | tr -d ' ')
    total_bases=\$(grep "Total assembled bases" ${sample_id}_Assembly-Stats.txt | cut -d: -f2 | tr -d ' ')

    echo -e "sample_id\\tTotal_trinity_genes\\tTotal_Trinity_Transcripts\\tPercentGC\\tContig_N50\\tMedian_contig_Length\\tAverage_contig\\tTotal_assembled_bases" > ${sample_id}_stats_summary.tsv
    echo -e "${sample_id}\\t\${total_genes}\\t\${total_transcripts}\\t\${percent_gc}\\t\${contig_n50}\\t\${median_len}\\t\${average_len}\\t\${total_bases}" >> ${sample_id}_stats_summary.tsv
    """
}