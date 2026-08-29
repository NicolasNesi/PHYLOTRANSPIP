#!/usr/bin/env nextflow

/*
 * Reciprocal Best Hit BLAST pour identifier les CDS orthologues
 */

process RCPBLAST {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/ReciprocalBLAST", mode: 'copy'

    input:
    tuple val(sample_id), path(cds)
    path proteome_ref

    output:
    tuple val(sample_id), path("${sample_id}_candidate_orthologs.fa"), emit: orthologs
    path "${sample_id}_candidate_orthologs.log", emit: log
    path "blastx_${sample_id}_cds_vs_Hg38Proteome.outfmt6.txt", emit: blastx_out
    path "tblastn_Hg38Proteome_vs_${sample_id}_cds.outfmt6.txt", emit: tblastn_out

    script:
    """
    makeblastdb -in ${cds} -dbtype 'nucl'

    blastx \\
        -query ${cds} \\
        -db ${proteome_ref} \\
        -out blastx_${sample_id}_cds_vs_Hg38Proteome.outfmt6.txt \\
        -evalue 1e-06 \\
        -max_target_seqs 1 \\
        -outfmt 6 \\
        -num_threads ${task.cpus}

    tblastn \\
        -db ${cds} \\
        -query ${proteome_ref} \\
        -out tblastn_Hg38Proteome_vs_${sample_id}_cds.outfmt6.txt \\
        -evalue 1e-06 \\
        -max_target_seqs 1 \\
        -outfmt 6 \\
        -num_threads ${task.cpus}

    analyze_blastPlus_topHit_coverage.pl \\
        blastx_${sample_id}_cds_vs_Hg38Proteome.outfmt6.txt \\
        ${cds} \\
        ${proteome_ref}
		
    conda deactivate
    conda activate ptp_env

    filter_reciprocal_blast_matches.pl \\
        -i blastx_${sample_id}_cds_vs_Hg38Proteome.outfmt6.txt.w_pct_hit_length \\
        -r tblastn_Hg38Proteome_vs_${sample_id}_cds.outfmt6.txt \\
        -d ${cds} \\
        -s ${sample_id} \\
        -o ${sample_id}_candidate_orthologs
    """
}