#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    folder       = null
    illumina     = 'IlluminaOutput'
    sra_list     = null
    samp_list    = null
    download_sra = true
    taxon        = 'phyllostomidae'
	ortholog_taxon = 'Microbat' 
    outdir       = 'results'
    help         = false
    h            = false
}

def helpMessage() {
    log.info """\
    ===========================================================================
    N E X T F L O W - P H Y L O T R A N S C R I P T O M I C S  P I P E L I N E 
    ( P H Y L O T R A N S P I P)
    ===========================================================================

    Usage:
      nextflow run main.nf [options]

    Required:
      --folder               <dir>       Root folder of the project
                                          (must include Scripts/, and if applicable IlluminaOutput/)

    Mode SRA download (--download_sra true):
      --sra_list       <file>     CSV file (in Scripts/): Genus_species_locality_code,SRRcode

    Mode reads déjà présents (--download_sra false):
      --samp_list      <file>      CSV file (in Scripts/): sample_id,...
      --illumina       <dir>       Name of subfolder with fastq.gz [default: IlluminaOutput]

   Optional:
      --outdir          <dir>       Output directory [default: results]
      --taxon           <name>      Prefix for the rRNA/mtDNA bowtie2 indexes of the species of interest
                                     in databases/bowtie/<taxon>_rRNA and <taxon>_mtDNA
                                     [default: phyllostomidae]
     --ortholog_taxon  <name>      Name of the Ensembl reference taxon used in the 1-to-1 ortholog list
                                     (Scripts/List_Ortholog_1to1_<ortholog_taxon>_Hsapiens_Ensembl.tsv)
                                     Format: one ENSG_genename per line
                                     [default: Microbat]

    Other:
      --help / --h                 Show this help message and exit


    Examples:
      # Whole pipeline
      nextflow run main.nf --download_sra true --sra_list list_sra.csv \\
          --taxon phyllostomidae \\
          --folder /dlocal/home/2019013/Data/ViroCaen/PHYLOTRANSPIP

      # No download of sra, use sra already downloaded or Illumina sequenced reads
      nextflow run main.nf --download_sra false --samp_list List_samples.csv \\
          --taxon phyllostomidae \\
          --folder /dlocal/home/2019013/Data/ViroCaen/PHYLOTRANSPIP

    Author       : ${workflow.manifest.author}
    Homepage     : ${workflow.manifest.homePage}
    Version      : ${workflow.manifest.version}
    ==================================================================
    """
    .stripIndent()
}

if (params.help || params.h) {
    helpMessage()
    System.exit(0)
}

log.info """\
    N E X T F L O W - P H Y L O T R A N S C R I P T O M I C S  P I P E L I N E
    ( P H Y L O T R A N S P I P)
    ==================================================================
    Output folder                   : ${params.outdir}
    SRA list                        : ${params.sra_list}
    Sample list                     : ${params.samp_list}
    Taxon (rRNA/mtDNA index)        : ${params.taxon}
	Ortholog reference taxon        : ${params.ortholog_taxon}
	
    Author                          : ${workflow.manifest.author}
    Homepage                        : ${workflow.manifest.homePage}
    Pipeline Version                : ${workflow.manifest.version}
    """
    .stripIndent()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WGETSRA }         from './modules/wgetsra.nf'
include { RCORRECT }        from './modules/rcorrect.nf'
include { TRIMGALORE }      from './modules/trimgalore.nf'
include { DECONTAM_MAP }    from './modules/decontammap.nf'
include { QCTRIM }          from './modules/qctrim.nf'
include { BOWTIEFILTER }    from './modules/bowtiefilter.nf'
include { FILTEROVERREP }   from './modules/filteroverrep.nf'
include { TRINITY }         from './modules/trinity.nf'
include { STATSTRINITY }    from './modules/statstrinity.nf'
include { REBOWTIE }        from './modules/rebowtie.nf'
include { TRANSRATE }       from './modules/transrate.nf'
include { CHIMERA }         from './modules/chimera.nf'
include { CORSET }          from './modules/corset.nf'
include { TRANSDECODER }    from './modules/transdecoder.nf'
include { CDHIT }           from './modules/cdhit.nf'
include { RCPBLAST }        from './modules/rcpblast.nf'
include { FILTERFASTA }     from './modules/filterfasta.nf'
include { MAKEMULTIFASTA }  from './modules/makemultifasta.nf'
include { GUIDANCE }        from './modules/guidance.nf'
include { FILTERALIGNMENTS }from './modules/filteralignments.nf'
include { GENETREE }        from './modules/genetree.nf'
include { TREESHRINK }      from './modules/treeshrink.nf'
include { BUSCO }           from './modules/busco.nf'
include { MITOGENOMES }     from './modules/mitogenomes.nf'
include { GENETREE2 }       from './modules/genetree2.nf'

// Variables set in functions called below
download_sra = true   // false if using Illumina raw reads from sequencing


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    if (!params.folder) {
        error "Must indicate --folder (root folder of the project)"
    }

    def scripts_dir     = "${params.folder}/Scripts"
    def illumina_dir    = "${params.folder}/${params.illumina}"
    def sra_list_path   = "${scripts_dir}/${params.sra_list}"
    def samp_list_path  = "${scripts_dir}/${params.samp_list}"
    def db_bowtie_dir   = "${scripts_dir}/databases/bowtie"
    def db_blast_dir    = "${scripts_dir}/databases/blast"

    // -- validation --
    if (params.download_sra && !params.sra_list) {
        error "Must indicate --sra_list when --download_sra true"
    }
    if (!params.download_sra && !params.samp_list) {
        error "Must indicate --samp_list when --download_sra false"
    }
    if (!params.download_sra && !file(illumina_dir).exists()) {
        error "Illumina folder missing: ${illumina_dir}"
    }
    if (!file("${db_bowtie_dir}/contam_bt2_index.1.bt2").exists()) {
        error "Contaminant bowtie2 index missing: ${db_bowtie_dir}/contam_bt2_index*"
    }
    if (!file("${db_blast_dir}/hg38Proteome.faa").exists()) {
        error "Blast protein reference missing: ${db_blast_dir}/hg38Proteome.faa"
    }

    def rrna_index = "${db_bowtie_dir}/${params.taxon}_rRNA"
    def mito_index = "${db_bowtie_dir}/${params.taxon}_mtDNA"

    if (!file("${rrna_index}.1.bt2").exists()) {
        error "rRNA bowtie2 index missing: ${rrna_index}*"
    }
    if (!file("${mito_index}.1.bt2").exists()) {
        error "mtDNA bowtie2 index missing: ${mito_index}*"
    }

    def contam_index_ch = Channel.fromPath("${db_bowtie_dir}/contam_bt2_index*").collect()
    def rrna_index_ch   = Channel.fromPath("${rrna_index}*").collect()
    def mito_index_ch   = Channel.fromPath("${mito_index}*").collect()
    def blast_db_ch     = Channel.fromPath("${db_blast_dir}/hg38Proteome.faa")
    def transdecoder_blastdb_ch = Channel.fromPath("${db_blast_dir}/hg38Proteome.p*").collect()

    if (params.download_sra) {
        Channel
            .fromPath(sra_list_path)
            .splitCsv(header: false)
            .map { row -> tuple(row[0], row[1]) }
            .set { sra_ch }

        WGETSRA(sra_ch)
        RCORRECT(WGETSRA.out.reads)

    } else {
        Channel
            .fromPath(samp_list_path)
            .splitCsv(header: true)
            .map { row ->
                def r1 = file("${illumina_dir}/${row.sample_id}_R1.fastq.gz")
                def r2 = file("${illumina_dir}/${row.sample_id}_R2.fastq.gz")
                tuple(row.sample_id, r1, r2)
            }
            .set { reads_ch }

        RCORRECT(reads_ch)
    }

    TRIMGALORE(RCORRECT.out.corrected)
    DECONTAM_MAP(TRIMGALORE.out.trimmed, contam_index_ch)

    QCTRIM(TRIMGALORE.out.trimmed, 'post_trimgalore')
    QCTRIM(DECONTAM_MAP.out.clean, 'post_decontam')

    BOWTIEFILTER(DECONTAM_MAP.out.clean, rrna_index_ch, mito_index_ch)
    BOWTIEFILTER.out.summary
        .collectFile(name: 'all_samples_contam_summary.tsv', storeDir: params.outdir, keepHeader: true, skip: 1)

    def qc_bowtiefilter = QCTRIM(BOWTIEFILTER.out.cleaned, 'post_bowtiefilter')
    def overrep_input = BOWTIEFILTER.out.cleaned.join(qc_bowtiefilter.zip)
    FILTEROVERREP(overrep_input)

    TRINITY(FILTEROVERREP.out.filtered)

    STATSTRINITY(TRINITY.out.assembly)
    STATSTRINITY.out.summary
        .collectFile(name: 'all_samples_trinity_stats.tsv', storeDir: params.outdir, keepHeader: true, skip: 1)

    REBOWTIE(TRINITY.out.assembly, rrna_index_ch)

    def transrate_input = REBOWTIE.out.cleaned.join(BOWTIEFILTER.out.cleaned)
    TRANSRATE(transrate_input)

    CHIMERA(TRANSRATE.out.filtered, blast_db_ch)

    def corset_input = CHIMERA.out.filtered.join(BOWTIEFILTER.out.cleaned)
    CORSET(corset_input)

    TRANSDECODER(CORSET.out.representative, transdecoder_blastdb_ch)

    CDHIT(TRANSDECODER.out.cds)
    CDHIT.out.summary
        .collectFile(name: 'all_samples_unique_transcripts.tsv', storeDir: params.outdir, keepHeader: true, skip: 1)

    RCPBLAST(CDHIT.out.unique, blast_db_ch)
	
    def ortholog_1to1_ch = Channel.fromPath("${scripts_dir}/List_Ortholog_1to1_${params.ortholog_taxon}_Hsapiens_Ensembl.tsv")

    FILTERFASTA(RCPBLAST.out.orthologs, ortholog_1to1_ch, genename_1to1_ch)
	
	def all_orthologs_ch = FILTERFASTA.out.onetoone
    .map { it[1] }      // ne garde que les fichiers, retire sample_id
    .flatten()           // aplati la liste de listes (chaque échantillon peut avoir plusieurs fichiers)
    .collect()            // regroupe tous les fichiers de tous les échantillons en un seul channel
	
	// -- convergence multi-échantillons pour la phylogénomique --
    MAKEMULTIFASTA(all_orthologs_ch, ortholog_1to1_ch)
	
	GUIDANCE(MAKEMULTIFASTA.out.multifasta.flatten())
	
}

    FILTERALIGNMENTS(GUIDANCE.out)
    GENETREE(FILTERALIGNMENTS.out)
    TREESHRINK(GENETREE.out)
    GENETREE2(TREESHRINK.out)

    BUSCO(TRINITY.out)                // souvent indépendant, évalue l'assemblage
    MITOGENOMES(TRIMGALORE.out)       // souvent à part, tourne sur les reads bruts
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
