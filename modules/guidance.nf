#!/usr/bin/env nextflow

/*
 * Alignement Guidance3 (PRANK), par gène orthologue
 */
 
process GUIDANCE {
    tag "$gene"

    publishDir "${params.outdir}/Guidance-Prank/MACSE", mode: 'copy', pattern: "*stats_TrimNonHomologous.txt"

    input:
    path multifasta

    output:
    tuple val(gene), path("${gene}_GuidancePRANK_algn.fasta"), emit: aligned, optional: true
    path "${gene}_stats_TrimNonHomologous.txt", optional: true

    script:
    gene = multifasta.name.replaceAll(/_MultiFasta_forGuidance_GeneName\.fasta$/, '')
    """
    n_seq=\$(grep -c '^>' ${multifasta})
    if [ "\$n_seq" -lt 6 ]; then
        echo "Skipping ${gene}: only \$n_seq sequences (<6)"
        exit 0
    fi

    macse -prog trimNonHomologousFragments \\
        -seq ${multifasta} \\
        -out_trim_info ${gene}_stats_TrimNonHomologous.txt \\
        -min_MEM_length 6 \\
        -min_trim_in 40 \\
        -min_trim_ext 20 \\
        -min_homology_to_keep_seq 0.1 \\
        -min_internal_homology_to_keep_seq 0.5

    guidance3 \\
        --seqFile ${gene}_MultiFasta_forGuidance_GeneName_NT.fasta \\
        --msaProgram PRANK \\
        --seqType codon \\
        --outDir ./${gene}_guidance_out \\
        --dataset ${gene} \\
        --seqCutoff 0.6 \\
        --colCutoff 0.93 \\
        --genCode 1 \\
        --bootstraps 100 \\
        --proc_num ${task.cpus}

    cp ${gene}_guidance_out/${gene}.PRANK.aln.With_Names ${gene}_GuidancePRANK_algn.fasta
    """
}