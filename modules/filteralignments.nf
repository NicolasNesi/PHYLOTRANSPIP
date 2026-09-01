#!/usr/bin/env nextflow

/*
 * Filtering and cleaning of alignments
 */
 
process FILTERALIGNMENTS {
    tag "$gene"

    publishDir "${params.outdir}/Guidance-Prank", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('_raw.fasta'))          "Alignment_Raw/${filename}"
        else if (filename.endsWith('_cleaned.fasta'))  "Alignment_Cleaned/${filename}"
        else if (filename.endsWith('.svg'))            "ImageAln/${filename}"
        else if (filename.contains('stats'))           "MACSE/${filename}"
        else null
    }

    input:
    tuple val(gene), path(alignment)

    output:
    tuple val(gene), path("${gene}_GuidancePRANK_algn_raw.fasta"), emit: raw, optional: true
    tuple val(gene), path("${gene}_GuidancePRANK_algn_cleaned.fasta"), emit: cleaned, optional: true
    path "${gene}_stats_trimEnd.txt", optional: true
    path "${gene}_GuidancePRANK_algn_raw_output.svg", optional: true
    path "${gene}_GuidancePRANK_algn_cleaned_output.svg", optional: true

    script:
    """
    seqkit sort --by-name ${alignment} > temp.file
    seqkit seq -w0 temp.file > ${gene}_GuidancePRANK_algn_sorted.fasta
    rm temp.file

    remove_short_fasta.py ${gene}_GuidancePRANK_algn_sorted.fasta ${gene}_GuidancePRANK_algn_trimshort.fasta 0.75

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

    mv ${gene}_GuidancePRANK_algn_sorted.fasta ${gene}_GuidancePRANK_algn_raw.fasta

    CIAlign --infile ${gene}_GuidancePRANK_algn_raw.fasta \\
        --outfile ${gene}_GuidancePRANK_algn_raw \\
        --keep_gaponly --plot_format svg --plot_output

    CIAlign --infile ${gene}_GuidancePRANK_algn_cleaned.fasta \\
        --outfile ${gene}_GuidancePRANK_algn_cleaned \\
        --keep_gaponly --plot_format svg --plot_output
    """
}