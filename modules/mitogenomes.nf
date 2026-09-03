#!/usr/bin/env nextflow

/*
 * Assemblage du mitogénome (NOVOPlasty) + annotation (MITOS2)
 */
 
process MITOGENOMES {
    tag "$sample_id"

    publishDir "${params.outdir}/CombinedSRR/${sample_id}/novoplasty", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_r1), path(reads_r2)
    path seed

    output:
    tuple val(sample_id), path("Circularized_assembly_1_${sample_id}.fas"), emit: assembly, optional: true
    path "mitos_out", emit: mitos_annotation, optional: true
    path "genes/*.fasta", emit: gene_fastas, optional: true
    tuple val(sample_id), path("${sample_id}_whole_mito.fasta"), emit: whole_mito, optional: true
    path "log_${sample_id}.txt", emit: log

    script:
    """
    cat > config_${sample_id}.txt << EOF
Project:
-----------------------
Project name          = ${sample_id}
Type                  = mito
Genome Range          = ${params.novoplasty_genome_range}
K-mer                 = ${params.novoplasty_kmer}
Max memory            =
Extended log          = 0
Save assembled reads  = yes
Seed Input            = ${seed}
Extend seed directly  = no
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = ${params.novoplasty_read_length}
Insert size           = ${params.novoplasty_insert_size}
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = ${reads_r1}
Reverse reads         = ${reads_r2}
Store Hash             =

Heteroplasmy:
-----------------------
MAF                    =
HP exclude list        =
PCR-free                =

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Reduce ambigious N's  =
Output path           = ./
EOF

    NOVOPlasty4.3.5.pl -c config_${sample_id}.txt > log_${sample_id}.txt 2>&1

    if [ -f "Circularized_assembly_${sample_id}.fasta" ]; then
        cp Circularized_assembly_${sample_id}.fasta Circularized_assembly_1_${sample_id}.fas
        echo "Status: circularized (single contig)" >> log_${sample_id}.txt

    elif ls Option_*_${sample_id}.fasta >/dev/null 2>&1; then
        first_option=\$(ls Option_*_${sample_id}.fasta | sort | head -1)
        cp "\$first_option" Circularized_assembly_1_${sample_id}.fas
        echo "Status: merged from multiple contigs, using \$first_option (check Merged_contigs_${sample_id}.txt manually)" >> log_${sample_id}.txt

    elif [ -f "Contigs_1_${sample_id}.fasta" ]; then
        cp Contigs_1_${sample_id}.fasta Circularized_assembly_1_${sample_id}.fas
        echo "Status: NOT circularized, using raw contig (may be incomplete/partial mitogenome)" >> log_${sample_id}.txt

    else
        echo "Status: NOVOPlasty produced no usable contig for ${sample_id}" >> log_${sample_id}.txt
    fi

    if [ -f "Circularized_assembly_1_${sample_id}.fas" ]; then
        len=\$(grep -v '^>' Circularized_assembly_1_${sample_id}.fas | tr -d '\\n' | wc -c)
        echo "Assembly length: \${len} bp" >> log_${sample_id}.txt
        if [ "\$len" -lt 14000 ]; then
            echo "WARNING: assembly shorter than expected mitogenome size, likely incomplete" >> log_${sample_id}.txt
        fi

        mkdir -p mitos_out genes
        runmitos.py \\
            --input Circularized_assembly_1_${sample_id}.fas \\
            --linear \\
            --code ${params.mitos_code} \\
            --refseqver ${params.mitos_refseqver} \\
            --refdir ${params.mitos_refdir} \\
            --outdir mitos_out

        split_mitos_genes.py mitos_out/result.fas ${sample_id} genes

        awk -v id="${sample_id}" '/^>/{print ">"id; next} {print}' Circularized_assembly_1_${sample_id}.fas > ${sample_id}_whole_mito.fasta
    else
        echo "NOVOPlasty did not produce a circularized assembly for ${sample_id}" >> log_${sample_id}.txt
    fi
    """
}