#!/usr/bin/env nextflow

/*
 * Statistiques par gène : saturation, treeness/RCV, long branch score, DVMC, taux évolutif,
 * outliers taxonomiques, biais de composition, stop codons internes, recombinaison (PhiPack)
 */
 
process STATSALGTREE {
    tag "$gene"

    publishDir "${params.outdir}/StatsAlgTree/${gene}", mode: 'copy'

    input:
    tuple val(gene), path(tree), path(alignment)

    output:
    path "${gene}_stats_summary.tsv", emit: stats
    path "${gene}_long_branch_score.txt", emit: lb_score_detail
    path "${gene}_outlier_taxa.txt", emit: outlier_taxa
    path "${gene}_noPrematureSTOP.fasta", emit: no_stop_codons
    path "${gene}_recombination.txt", emit: recombination

    script:
    """
    saturation_raw=\$(phykit saturation -a ${alignment} -t ${tree})
    saturation=\$(echo "\$saturation_raw" | awk '{print \$1}')
    saturation_dev_from_1=\$(echo "\$saturation_raw" | awk '{print \$2}')

    toverr_raw=\$(phykit treeness_over_rcv -a ${alignment} -t ${tree})
    treeness=\$(echo "\$toverr_raw" | awk '{print \$1}')
    rcv_toverr=\$(echo "\$toverr_raw" | awk '{print \$2}')
    toverr=\$(echo "\$toverr_raw" | awk '{print \$3}')

    phykit long_branch_score ${tree} > ${gene}_long_branch_score.txt
    lbscore_median=\$(grep '^median:' ${gene}_long_branch_score.txt | awk '{print \$2}')

    dvmc=\$(phykit degree_of_violation_of_a_molecular_clock ${tree})
    evorate=\$(phykit evolutionary_rate ${tree})
    rcv=\$(phykit relative_composition_variability ${alignment})

    phykit alignment_outlier_taxa ${alignment} > ${gene}_outlier_taxa.txt
    n_outliers=\$(tail -n +3 ${gene}_outlier_taxa.txt | wc -l)

    echo -e "gene\\tsaturation\\tsaturation_dev_from_1\\ttreeness\\trcv_toverr\\ttreeness_over_rcv\\tlong_branch_score_median\\tdvmc\\tevolutionary_rate\\trcv\\tn_outlier_taxa" > ${gene}_stats_summary.tsv
    echo -e "${gene}\\t\${saturation}\\t\${saturation_dev_from_1}\\t\${treeness}\\t\${rcv_toverr}\\t\${toverr}\\t\${lbscore_median}\\t\${dvmc}\\t\${evorate}\\t\${rcv}\\t\${n_outliers}" >> ${gene}_stats_summary.tsv

    premature_stops.pl -i ${alignment} -o ${gene}_noPrematureSTOP.fasta

    n_in=\$(grep -c '^>' ${alignment})
    n_out=\$(grep -c '^>' ${gene}_noPrematureSTOP.fasta)
    echo -e "${gene}\\tn_seq_input\\t\${n_in}\\tn_seq_no_internal_stop\\t\${n_out}" >> ${gene}_stats_summary.tsv

    Phi -r ${alignment} -p 1000 -o -v > ${gene}_recombination.txt
    """
}