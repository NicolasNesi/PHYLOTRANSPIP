#!/usr/bin/env nextflow

/*
 * Filtre les paires de reads contenant une séquence surreprésentée (>0.1%)
 */
 
process FILTEROVERREP {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/overrep_filter", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2), path(fastqc_zips)

    output:
    tuple val(sample_id), path("${sample_id}_R1.overrep_filtered.fq.gz"), path("${sample_id}_R2.overrep_filtered.fq.gz"), emit: filtered
    path "${sample_id}_overrep_seqs.fasta", emit: overrep_seqs

    script:
    """
    mkdir unzipped
    for z in ${fastqc_zips}; do unzip -o -q \$z -d unzipped; done

    find unzipped -name fastqc_data.txt | xargs awk '
        BEGIN { FS="\\t" }
        /^>>Overrepresented sequences/ { flag=1; next }
        /^>>END_MODULE/ { flag=0 }
        flag && \$0 !~ /^#/ && NF>0 {
            if (\$3+0 > 0.1) print ">overrep_"NR"\\n"\$1
        }
    ' > ${sample_id}_overrep_seqs.fasta

    if [ -s ${sample_id}_overrep_seqs.fasta ]; then
        seqkit grep -s -r -f ${sample_id}_overrep_seqs.fasta ${reads_r1} | seqkit seq -n -i > ids_r1.txt
        seqkit grep -s -r -f ${sample_id}_overrep_seqs.fasta ${reads_r2} | seqkit seq -n -i > ids_r2.txt
        cat ids_r1.txt ids_r2.txt | sort -u > exclude_ids.txt

        seqkit grep -v -f exclude_ids.txt ${reads_r1} -o ${sample_id}_R1.overrep_filtered.fq.gz
        seqkit grep -v -f exclude_ids.txt ${reads_r2} -o ${sample_id}_R2.overrep_filtered.fq.gz
    else
        cp ${reads_r1} ${sample_id}_R1.overrep_filtered.fq.gz
        cp ${reads_r2} ${sample_id}_R2.overrep_filtered.fq.gz
    fi
    """
}