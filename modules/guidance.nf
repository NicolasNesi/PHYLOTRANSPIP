#!/usr/bin/env nextflow
/*
 * Alignement Guidance3 (PRANK), générique — gènes nucléaires (codon) ou mitochondriaux (codon/nuc)
 */
process GUIDANCE {
    tag "$gene"

    publishDir "${params.outdir}/Guidance-Prank/MACSE", mode: 'copy', pattern: "*stats_TrimNonHomologous.txt"

    input:
    tuple path(multifasta), val(seqType)

    output:
    tuple val(gene), path("${gene}_GuidancePRANK_algn.fasta"), val(seqType), emit: aligned, optional: true
    path "${gene}_stats_TrimNonHomologous.txt", optional: true

    script:
    gene = multifasta.name.replaceAll(/(_MultiFasta_forGuidance_GeneName|_MitoMultiFasta)\.fasta$/, '')
    def extra_args = (seqType == 'codon') ? '--genCode 1' : ''
    """
    n_seq=\$(grep -c '^>' ${multifasta})
    if [ "\$n_seq" -lt ${params.min_seq_guidance} ]; then
        echo "Skipping ${gene}: only \$n_seq sequences (<${params.min_seq_guidance})"
        exit 0
    fi

    base=\$(basename ${multifasta} .fasta)
    seqfile=${multifasta}

    if [ "${seqType}" = "codon" ]; then
        macse -prog trimNonHomologousFragments \\
            -seq ${multifasta} \\
            -out_trim_info ${gene}_stats_TrimNonHomologous.txt \\
            -min_MEM_length 6 \\
            -min_trim_in 40 \\
            -min_trim_ext 20 \\
            -min_homology_to_keep_seq 0.1 \\
            -min_internal_homology_to_keep_seq 0.5
        seqfile=\${base}_NT.fasta
    fi

    guidance3 \\
        --seqFile \${seqfile} \\
        --msaProgram PRANK \\
        --seqType ${seqType} \\
        --outDir ./${gene}_guidance_out \\
        --dataset ${gene} \\
        --seqCutoff 0.6 \\
        --colCutoff 0.93 \\
        ${extra_args} \\
        --bootstraps 100 \\
        --proc_num ${task.cpus}

    cp ${gene}_guidance_out/${gene}.PRANK.aln.With_Names ${gene}_GuidancePRANK_algn.fasta
    """
}