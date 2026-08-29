#!/usr/bin/env nextflow

/*
 * Alignement Guidance3 (PRANK) + nettoyage MACSE/clipKit/CIAlign, par gène orthologue
 */
 
process GUIDANCE {
    tag "$gene"

    publishDir "${params.outdir}/Guidance-Prank", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('_raw.fasta'))            "Alignment_Raw/${filename}"
        else if (filename.endsWith('_cleaned.fasta'))    "Alignment_Cleaned/${filename}"
        else if (filename.endsWith('.svg'))              "ImageAln/${filename}"
        else if (filename.contains('stats'))             "MACSE/${filename}"
        else null
    }

    input:
    path multifasta

    output:
    tuple val(gene), path("${gene}_GuidancePRANK_algn_raw.fasta"), emit: raw, optional: true
    tuple val(gene), path("${gene}_GuidancePRANK_algn_cleaned.fasta"), emit: cleaned, optional: true
    path "${gene}_stats_TrimNonHomologous.txt", optional: true
    path "${gene}_stats_trimEnd.txt", optional: true
    path "${gene}_GuidancePRANK_algn_raw_output.svg", optional: true
    path "${gene}_GuidancePRANK_algn_cleaned_output.svg", optional: true

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

    conda deactivate
    conda activate ptp_env

    seqkit sort --by-name ${gene}_GuidancePRANK_algn.fasta > temp.file
    seqkit seq -w0 temp.file > ${gene}_GuidancePRANK_algn_sorted.fasta
    mv ${gene}_GuidancePRANK_algn_sorted.fasta ${gene}_GuidancePRANK_algn.fasta
    rm temp.file

    remove_short_fasta.py ${gene}_GuidancePRANK_algn.fasta ${gene}_GuidancePRANK_algn_trimshort.fasta 0.75

    macse -prog trimAlignment \\
        -align ${gene}_GuidancePRANK_algn_trimshort.fasta \\
        -min_percent_NT_at_ends 0.5 \\
        -out_trim_info ${gene}_stats_trimEnd.txt

    clipkit ${gene}_GuidancePRANK_algn_trimshort_NT.fasta \\
        -m gappy \\
        --gaps 0.9 \\
        --codon \\
        --output ${gene}_GuidancePRANK_algn_trimshort.trimmed.fa

    remove_short_fasta.py ${gene}_GuidancePRANK_algn_trimshort.trimmed.fa ${gene}_GuidancePRANK_algn_trimshort.trimmed2.fa 0.80

    seqkit sort --by-name ${gene}_GuidancePRANK_algn_trimshort.trimmed2.fa > temp2.file
    seqkit seq -w0 temp2.file > ${gene}_GuidancePRANK_algn_cleaned.fasta
    rm temp2.file

    mv ${gene}_GuidancePRANK_algn.fasta ${gene}_GuidancePRANK_algn_raw.fasta

    CIAlign --infile ${gene}_GuidancePRANK_algn_raw.fasta \\
        --outfile ${gene}_GuidancePRANK_algn_raw \\
        --keep_gaponly --plot_format svg --plot_output

    CIAlign --infile ${gene}_GuidancePRANK_algn_cleaned.fasta \\
        --outfile ${gene}_GuidancePRANK_algn_cleaned \\
        --keep_gaponly --plot_format svg --plot_output
    """
}