#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.folder               = null
params.illumina             = 'IlluminaOutput'
params.sra_list             = null
params.samp_list            = null
params.download_sra         = true

params.taxon                = 'phyllostomidae'
params.ortholog_taxon       = 'Microbat'
params.busco_lineage        = 'mammalia_odb12'
params.min_species_genetree = 10

params.novoplasty_seed         = null
params.novoplasty_genome_range = '12000-22000'
params.novoplasty_kmer         = 23
params.novoplasty_read_length  = 151
params.novoplasty_insert_size  = 300
params.mitos_refdir            = '/dlocal/home/2019013/Databases/MITOS2'
params.mitos_refseqver         = 'refseq63m'
params.mitos_code              = 2

params.mode             = 'full'   // full | mito_only | per_sample | downstream
params.previous_outdirs = null     // requis pour --mode downstream

params.outdir = 'results'
params.help   = false
params.h      = false

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
    log.info """\
    ===========================================================================
    N E X T F L O W - P H Y L O T R A N S C R I P T O M I C S  P I P E L I N E 
    ( P H Y L O T R A N S P I P)
    ===========================================================================

    Usage:
      nextflow run main.nf [options]

    Modes (--mode):
      full        Whole pipeline: per-sample assembly through phylogenomics (default)
      mito_only   Reads -> trimming -> NOVOPlasty/MITOS2 mitogenome assembly only
      per_sample  Per-sample steps only (assembly through 1-to-1 orthologs), stops
                  before the multi-sample convergence step. Use this to process
                  independent sequencing/SRA batches separately.
      downstream  Multi-sample phylogenomics steps only, starting from the
                  1-to-1 orthologs produced by one or more previous --mode per_sample
                  runs (see --previous_outdirs).

    Required:
      --folder               <dir>       Project root folder
                                          (must contain Scripts/main.nf, Scripts/databases/,
                                          and, if applicable, IlluminaOutput/)

    Mode SRA download (--download_sra true) :
      --sra_list       <file>      CSV file (in Scripts/) : Genus_species_locality_code,SRRcode

    Mode reads already present (--download_sra false) :
      --samp_list      <file>      CSV file (in Scripts/) : sample_id,...
      --illumina       <dir>       Fastq.gz subfolder name [default: IlluminaOutput]

    Databases (in Scripts/databases/) :
      --taxon           <name>     Prefix for the rRNA/mtDNA bowtie2 indexes of the species
                                    of interest, in databases/bowtie/<taxon>_rRNA and
                                    <taxon>_mtDNA [default: phyllostomidae]
      --ortholog_taxon  <name>     Name of the Ensembl reference taxon used in the 1-to-1
                                    ortholog list (Scripts/List_Ortholog_1to1_<ortholog_taxon>_Hsapiens_Ensembl.tsv)
                                    [default: Microbat]
      --busco_lineage   <name>     BUSCO lineage dataset name (online mode, auto-downloaded)
                                    [default: mammalia_odb12]

    Mitogenome assembly (NOVOPlasty + MITOS2) :
      --novoplasty_seed          <fasta>   Mitochondrial seed sequence (e.g. COI), required
      --novoplasty_genome_range  <min-max> Expected mitogenome size range [default: 12000-22000]
      --novoplasty_kmer          <int>     NOVOPlasty k-mer size [default: 23]
      --novoplasty_read_length   <int>     Read length [default: 151]
      --novoplasty_insert_size   <int>     Insert size [default: 300]
      --mitos_refdir             <dir>     MITOS2 reference database directory
      --mitos_refseqver          <name>    MITOS2 reference sequence version [default: refseq63m]
      --mitos_code               <int>     Genetic code for MITOS2 [default: 2]

    Downstream mode only:
      --previous_outdirs   <list>   Comma-separated list of previous --mode per_sample
                                     output directories to combine, e.g.:
                                     results_batch1,results_batch2

    Optional:
      --outdir                <dir>   Output directory [default: results]
      --min_species_genetree   <int>  Minimum number of species required to build a gene tree
                                       [default: 10]

    Other:
      --help / --h                 Show this help message and exit

    Examples:
      # Whole pipeline, single batch
      nextflow run main.nf --mode full --download_sra false --samp_list List_samples.csv \\
          --novoplasty_seed Seed_COI.fasta --folder /dlocal/.../PHYLOTRANSPIP

      # Per-sample only, first sequencing batch
      nextflow run main.nf --mode per_sample --samp_list Batch1_samples.csv \\
          --novoplasty_seed Seed_COI.fasta --outdir results_batch1 --folder /dlocal/.../PHYLOTRANSPIP

      # Per-sample only, SRA batch
      nextflow run main.nf --mode per_sample --download_sra true --sra_list Batch2_sra.csv \\
          --novoplasty_seed Seed_COI.fasta --outdir results_batch2 --folder /dlocal/.../PHYLOTRANSPIP

      # Combine batches for phylogenomics
      nextflow run main.nf --mode downstream --previous_outdirs results_batch1,results_batch2 \\
          --outdir results_final --folder /dlocal/.../PHYLOTRANSPIP

      # Mitogenome only
      nextflow run main.nf --mode mito_only --download_sra false --samp_list List_samples.csv \\
          --novoplasty_seed Seed_COI.fasta --folder /dlocal/.../PHYLOTRANSPIP

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
    Mode                             : ${params.mode}
    Output folder                    : ${params.outdir}
    SRA list                         : ${params.sra_list}
    Sample list                      : ${params.samp_list}
    Taxon (rRNA/mtDNA index)         : ${params.taxon}
    Ortholog reference taxon         : ${params.ortholog_taxon}
    BUSCO lineage                    : ${params.busco_lineage}
    Author                           : ${workflow.manifest.author}
    Homepage                         : ${workflow.manifest.homePage}
    Pipeline Version                 : ${workflow.manifest.version}
    """
    .stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WGETSRA }          from './modules/wgetsra.nf'
include { RCORRECT }         from './modules/rcorrect.nf'
include { TRIMGALORE }       from './modules/trimgalore.nf'
include { DECONTAM_MAP }     from './modules/decontam_map.nf'
include { QCTRIM }           from './modules/qctrim.nf'
include { BOWTIEFILTER }     from './modules/bowtiefilter.nf'
include { FILTEROVERREP }    from './modules/filteroverrep.nf'
include { TRINITY }          from './modules/trinity.nf'
include { STATSTRINITY }     from './modules/statstrinity.nf'
include { REBOWTIE }         from './modules/rebowtie.nf'
include { TRANSRATE }        from './modules/transrate.nf'
include { CHIMERA }          from './modules/chimera.nf'
include { CORSET }           from './modules/corset.nf'
include { TRANSDECODER }     from './modules/transdecoder.nf'
include { CDHIT }            from './modules/cdhit.nf'
include { RCPBLAST }         from './modules/rcpblast.nf'
include { FILTERFASTA }      from './modules/filterfasta.nf'
include { BUSCO }            from './modules/busco.nf'
include { MITOGENOMES }      from './modules/mitogenomes.nf'

include { MAKEMULTIFASTA }   from './modules/makemultifasta.nf'
include { GUIDANCE }         from './modules/guidance.nf'
include { FILTERALIGNMENTS } from './modules/filteralignments.nf'
include { GENETREE }         from './modules/genetree.nf'
include { QCTREE }           from './modules/qctree.nf'
include { CONCATENATION }    from './modules/concatenation.nf'
include { GENETREE2 }        from './modules/genetree2.nf'
include { STATSALGTREE }     from './modules/statsalgtree.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    if (!params.folder) {
        error "Must indicate --folder (root folder of the project)"
    }

    def scripts_dir    = "${params.folder}/Scripts"
    def illumina_dir   = "${params.folder}/${params.illumina}"
    def sra_list_path  = "${scripts_dir}/${params.sra_list}"
    def samp_list_path = "${scripts_dir}/${params.samp_list}"
    def db_bowtie_dir  = "${scripts_dir}/databases/bowtie"
    def db_blast_dir   = "${scripts_dir}/databases/blast"

    /*
     * ==========================================================
     *  MODE: downstream — multi-sample phylogenomics only
     * ==========================================================
     */
    if (params.mode == 'downstream') {

        if (!params.previous_outdirs) {
            error "Must indicate --previous_outdirs when --mode downstream (comma-separated list)"
        }

        def outdirs = params.previous_outdirs.split(',').collect { it.trim() }

        def ortholog_1to1_ch = Channel.fromPath("${scripts_dir}/List_Ortholog_1to1_${params.ortholog_taxon}_Hsapiens_Ensembl.tsv")

        def all_orthologs_ch = Channel
            .fromPath(outdirs.collect { "${it}/CombinedSRR/*/orthologs/*_CDS_OnetoOne.fasta" })
            .flatten()
            .collect()

        MAKEMULTIFASTA(all_orthologs_ch, ortholog_1to1_ch)
        GUIDANCE(MAKEMULTIFASTA.out.multifasta.flatten())
        FILTERALIGNMENTS(GUIDANCE.out.aligned)
        GENETREE(FILTERALIGNMENTS.out.cleaned)

        def qctree_alignments = FILTERALIGNMENTS.out.cleaned.map { it[1] }.collect()
        def qctree_trees      = GENETREE.out.tree.map { it[1] }.collect()

        QCTREE(qctree_alignments, qctree_trees)
        CONCATENATION(QCTREE.out.filtered_alignments, QCTREE.out.matrix_ok_forconcat)

        def genetree2_input = QCTREE.out.matrix_ok_forconcat
            .flatten()
            .map { f -> tuple(f.name.replaceAll(/_MatrixOK_ForConcat\.fasta$/, ''), f) }

        GENETREE2(genetree2_input)

        def statsalgtree_input = GENETREE2.out.tree
            .join(
                QCTREE.out.matrix_ok_forconcat
                    .flatten()
                    .map { f -> tuple(f.name.replaceAll(/_MatrixOK_ForConcat\.fasta$/, ''), f) }
            )

        STATSALGTREE(statsalgtree_input)
        STATSALGTREE.out.stats
            .collectFile(name: 'all_genes_stats_summary.tsv', storeDir: params.outdir, keepHeader: true, skip: 1)

        return
    }

    /*
     * ==========================================================
     *  MODES: full / mito_only / per_sample — validations communes
     * ==========================================================
     */
    if (params.download_sra && !params.sra_list) {
        error "Must indicate --sra_list when --download_sra true"
    }
    if (!params.download_sra && !params.samp_list) {
        error "Must indicate --samp_list when --download_sra false"
    }
    if (!params.download_sra && !file(illumina_dir).exists()) {
        error "Illumina folder missing: ${illumina_dir}"
    }
    if (!params.novoplasty_seed) {
        error "Must indicate --novoplasty_seed (mitochondrial seed fasta, e.g. COI)"
    }

    def seed_ch = Channel.fromPath(params.novoplasty_seed)

    // -- construction des reads (SRA ou locaux), commune à tous les modes --
    def reads_ch
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

    /*
     * ==========================================================
     *  MODE: mito_only — s'arrête ici
     * ==========================================================
     */
    if (params.mode == 'mito_only') {
        MITOGENOMES(TRIMGALORE.out.trimmed, seed_ch)
        return
    }

    /*
     * ==========================================================
     *  MODES: full / per_sample — assemblage complet par échantillon
     * ==========================================================
     */
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
    def ortholog_1to1_ch = Channel.fromPath("${scripts_dir}/List_Ortholog_1to1_${params.ortholog_taxon}_Hsapiens_Ensembl.tsv")

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

    FILTERFASTA(RCPBLAST.out.orthologs, ortholog_1to1_ch)

    // -- QC par échantillon, en parallèle, quel que soit full/per_sample --
    def busco_input = TRINITY.out.assembly.join(CHIMERA.out.filtered)
    BUSCO(busco_input)

    MITOGENOMES(TRIMGALORE.out.trimmed, seed_ch)

    /*
     * ==========================================================
     *  MODE: per_sample — s'arrête ici
     * ==========================================================
     */
    if (params.mode == 'per_sample') {
        return
    }

    /*
     * ==========================================================
     *  MODE: full — enchaîne directement sur la phylogénomique
     * ==========================================================
     */
    def all_orthologs_ch = FILTERFASTA.out.onetoone.map { it[1] }.flatten().collect()

    MAKEMULTIFASTA(all_orthologs_ch, ortholog_1to1_ch)
    GUIDANCE(MAKEMULTIFASTA.out.multifasta.flatten())
    FILTERALIGNMENTS(GUIDANCE.out.aligned)
    GENETREE(FILTERALIGNMENTS.out.cleaned)

    def qctree_alignments = FILTERALIGNMENTS.out.cleaned.map { it[1] }.collect()
    def qctree_trees      = GENETREE.out.tree.map { it[1] }.collect()

    QCTREE(qctree_alignments, qctree_trees)
    CONCATENATION(QCTREE.out.filtered_alignments, QCTREE.out.matrix_ok_forconcat)

    def genetree2_input = QCTREE.out.matrix_ok_forconcat
        .flatten()
        .map { f -> tuple(f.name.replaceAll(/_MatrixOK_ForConcat\.fasta$/, ''), f) }

    GENETREE2(genetree2_input)

    def statsalgtree_input = GENETREE2.out.tree
        .join(
            QCTREE.out.matrix_ok_forconcat
                .flatten()
                .map { f -> tuple(f.name.replaceAll(/_MatrixOK_ForConcat\.fasta$/, ''), f) }
        )

    STATSALGTREE(statsalgtree_input)
    STATSALGTREE.out.stats
        .collectFile(name: 'all_genes_stats_summary.tsv', storeDir: params.outdir, keepHeader: true, skip: 1)
}