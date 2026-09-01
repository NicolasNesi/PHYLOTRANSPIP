# Phylotranscriptomics Pipeline PHYLOTRANSPIP

! VERSION NEXTFLOW IN PROGRESS !

Pipeline: from RNAseq reads (SRA or Ilumina raw reads) to matrices of orthologs coding sequencing genes (CDSs) for phylogenomics and comparative genomics analyses

To add:
 - companion script to create the contaminants bt2 index (human genome, UniVec, 37 bacteria, 7 fungi, 2 animals; list from PhyloProcessR pipeline)
 - companion script to get list of SRA available for your group of interest (e.g. Phyllostomidae) and then create your list_sra.csv file
 - list of softs needed and yml files to create mamba environments
 - create an apptainer/singularity?
 - calculate saturation for a given tree and alignment (Philippe et al. 2011) -> phykit saturation --alignment <alignment> --tree <tree>
 - calculate Treeness (signal-to-noise among branch lengths) divided by Relative Composition Variability (measure of composition bias) (treeness/RCV) -> phykit treeness_over_rcv --alignment <alignment> --tree <tree>
 - calculate long branch score -> phykit long_branch_score <tree>
 - phykit degree_of_violation_of_a_molecular_clock <tree>
 - phykit evolutionary_rate <tree>
 - pipeline flowchart
 - option pipeline: stop at FilterFasta step. continue to create CDSs matrix using multiple sequencing RUN / SRA download
 - create partition files (raxml format) for iqtree: per gene, per gene and codon position
 - option to run only mitochondrion step
 - QC matrix for Positive Selection Tests (internal STOP), recombination (phipack) etc...
 - Detailed README.md to do on final version

