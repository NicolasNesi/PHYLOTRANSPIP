# Phylotranscriptomics Pipeline PHYLOTRANSPIP

! VERSION NEXTFLOW IN PROGRESS !

Pipeline: from RNAseq reads (SRA or Ilumina raw reads) to matrices of orthologs coding sequencing genes (CDSs) for phylogenomics and comparative genomics analyses

To add:
 - companion script to create the contaminants bt2 index (human genome, UniVec, 37 bacteria, 7 fungi, 2 animals; list from PhyloProcessR pipeline)
 - companion script to get list of SRA available for your group of interest (e.g. Phyllostomidae) and then create your list_sra.csv file
 - list of softs needed and yml files to create mamba environments
 - create an apptainer/singularity?
 - calculate site saturation per gene alignment (Philippe et al. 2011) -> phykit sat -a gene1_msa.aln -t gene1_genetree.tre
 - calculate Treeness (signal-to-noise among branch lengths) divided by relative composition variability (measure of composition bias) (treeness/RCV) -> phykit toverr -a gene1_msa.aln -t gene1_genetree.tre
 - calulate long branch score -> phykit lb_score -t gene1_genetree.tre
 - pipeline flowchart
 - option pipeline: stop at FilterFasta step. continue to create CDSs matrix using multiple sequencing RUN / SRA download
 - create partition files (raxml format) for iqtree: per gene, per gene and codon position
 - option to run only mitochondrion step
 - Complete README.md to do on final version

