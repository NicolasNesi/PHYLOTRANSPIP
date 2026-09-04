# PhyloTransPip

**Phylotranscriptomics Pipeline** — a Nextflow pipeline from RNA-seq reads to CDSs for phylogenomics and comparative genomics analyses.

From raw Illumina reads or SRA accessions, the pipeline performs read cleaning and decontamination, *de novo* transcriptome assembly with Trinity, chimera and redundancy filtering, ortholog identification against a human reference proteome, multi-sample gene alignment, gene tree inference, and phylogenomic quality control (saturation, treeness/RCV, recombination, etc.). It also assembles and annotates a mitochondrial genome per sample, and builds concatenated mitogenome matrices (protein-coding genes, tRNA, rRNA) across all samples.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 25.10
- A SLURM cluster (the pipeline is configured for SLURM by default; adjust `nextflow.config` for other executors)
- `mamba`/`conda`

### Conda/mamba environments

The pipeline expects the following environments to already exist (see `nextflow.config` for the exact activation commands used by each step):

| Environment | Used for |
|---|---|
| `ptp_env` | amas, bowtie2, cd-hit, cialign, clipkit, corset, csvtk, fastqc, ffq, IQ-TREE2, MACSE, mafft, NOVOplasty, perl-bioperl, phykit, PhiPack, pysradb, rcorrector, salmon, seqkit, sratoolkit, treeshrink, trim-galore, blast+ |
| `trinity_env` | analyze_blastPlus_topHit_coverage.pl, Trinity assembly, TrinityStats.pl, blast+ (recip. blast) |
| `busco_env` | BUSCO |
| `transdecoder_env` | TransDecoder |
| `python2_env` | Legacy Python 2 scripts (chimera detection), blast+ |
| `guidance3_env` | Guidance3 (alignment), MACSE, PRANK |
| `phylo_env` | MITOS2 |

### Custom scripts (`Scripts/bin/`)

The following custom scripts must be present in `Scripts/bin/`, executable, with a correct shebang. Nextflow automatically adds this folder to the `PATH` of every process:

- `fixrcorrector.py`, `seq.py`
- `run_chimera_detection.py`, `detect_chimera_from_blastx_modifed.py`
- `filter_transcripts_transrate.py`
- `filter_corset_output.py`
- `transdecoder_wrapper.py`
- `filter_reciprocal_blast_matches.pl`
- `premature_stops.pl`
- `remove_short_fasta.py`
- `keep_longest.py`
- `split_mitos_genes.py` — splits a MITOS2 `result.fas` output into one FASTA file per gene, header renamed to `<sample_id>_<gene>`

---

## Directory structure (required)

The pipeline expects a fixed project layout, referenced everywhere via `--folder`:

```
<folder>/
├── Scripts/
│   ├── main.nf
│   ├── nextflow.config
│   ├── modules/                       # one .nf file per process
│   ├── bin/                           # custom scripts (see above)
│   ├── databases/
│   │   ├── bowtie/
│   │   │   ├── contam_bt2_index.*.bt2         # UniVec + bacteria/fungi/animal contaminant index
│   │   │   ├── <taxon>_rRNA.*.bt2             # rRNA index for the species of interest
│   │   │   └── <taxon>_mtDNA.*.bt2            # mitochondrial index for the species of interest
│   │   └── blast/
│   │       └── hg38Proteome.*                 # human proteome FASTA + blastp index (.phr/.pin/.psq)
│   ├── List_Ortholog_1to1_<ortholog_taxon>_Hsapiens_Ensembl.tsv   # one ENSG_GeneName per line
│   ├── List_Mito_Genes.txt                     # one mitochondrial gene/tRNA/rRNA/OH/OL name per line (as annotated by MITOS2)
│   ├── List_Mito_PCG.txt                       # subset of List_Mito_Genes.txt: the 13 protein-coding genes only
│   ├── List_sra.csv                            # for --download_sra true
│   └── List_samples.csv                        # for --download_sra false
└── IlluminaOutput/                    # <sample_id>_R1.fastq.gz / _R2.fastq.gz (if reads are local)
```

---

## Pipeline steps

| Step | Module | Tool(s) |
|---|---|---|
| SRA download | `WGETSRA` | ffq, sratoolkit |
| Error correction | `RCORRECT` | Rcorrector |
| Trimming | `TRIMGALORE` | Trim Galore |
| Contaminant decontamination | `DECONTAM_MAP` | bowtie2 vs UniVec + bacteria/fungi/animal DB |
| QC | `QCTRIM` | FastQC (run at multiple stages) |
| rRNA/mtDNA filtering | `BOWTIEFILTER` | bowtie2 vs species-specific rRNA/mtDNA |
| Overrepresented read filtering | `FILTEROVERREP` | FastQC overrepresented sequences + seqkit |
| *De novo* assembly | `TRINITY` | Trinity |
| Assembly stats | `STATSTRINITY` | TrinityStats.pl |
| rRNA contig removal | `REBOWTIE` | bowtie2 |
| Assembly quality / filtering | `TRANSRATE` | Transrate + custom score filter |
| Chimera detection | `CHIMERA` | blastx vs human proteome |
| Transcript clustering | `CORSET` | salmon + Corset |
| CDS prediction | `TRANSDECODER` | TransDecoder + blastp |
| Redundancy removal | `CDHIT` | cd-hit-est |
| Reciprocal best hit | `RCPBLAST` | blastx/tblastn vs human proteome |
| 1-to-1 ortholog filtering | `FILTERFASTA` | header editing, length/stop-codon filtering, longest-seq/dedup |
| BUSCO completeness | `BUSCO` | BUSCO (raw + filtered assembly) |
| Mitogenome assembly + annotation | `MITOGENOMES` | NOVOPlasty + MITOS2, per-gene FASTA split, whole-genome FASTA |
| **— multi-sample convergence, nuclear —** | | |
| Multifasta per gene | `MAKEMULTIFASTA` | concatenates 1-to-1 orthologs across samples |
| Alignment | `GUIDANCE` | MACSE (trim non-homologous, codon only) + Guidance3/PRANK — generic, reused for mito genes |
| Alignment cleaning | `FILTERALIGNMENTS` | MACSE, clipKit, CIAlign — generic, reused for mito genes; also produces a positions-1+2-only FASTA (`clipkit -m c3`) for codon alignments |
| Gene tree (pass 1) | `GENETREE` | IQ-TREE2 (≥10 taxa filter) |
| Tree/alignment QC | `QCTREE` | TreeShrink |
| Matrix concatenation | `CONCATENATION` | AMAS (nexus + RAxML partitions, by gene / by gene×codon) |
| Gene tree (pass 2) | `GENETREE2` | IQ-TREE2 on TreeShrink-filtered matrices |
| Alignment/tree statistics | `STATSALGTREE` | PhyKIT (saturation, treeness/RCV, LB score, DVMC, evolutionary rate, outlier taxa, RCV), PhiPack recombination test, internal stop codon check |
| **— multi-sample convergence, mitogenome —** | | |
| Multifasta per mito gene + whole genome | `MITO_MULTIFASTA` | concatenates per-gene and whole-mitogenome FASTAs across samples |
| Whole mitogenome alignment (optional) | `MITO_WHOLE_ALIGN` | MAFFT — not used downstream; QC/visualization only, disabled by default (see Notes) |
| Alignment (mito genes) | `GUIDANCE` | same module as nuclear genes; `codon` for the 13 protein-coding genes, `nuc` for tRNA/rRNA |
| Alignment cleaning (mito genes) | `FILTERALIGNMENTS` | same module as nuclear genes |
| Mitogenome matrix concatenation | `MITO_CONCAT` | AMAS — full matrix (PCG+tRNA+rRNA, excl. OH/OL) by gene and by gene×codon; PCG-only matrix by gene×codon; positions-1+2-only versions of both |

---

## Usage

```bash
nextflow run main.nf [options]
```

### Modes (`--mode`)

| Mode | Description |
|---|---|
| `full` (default) | Whole pipeline: per-sample assembly through nuclear + mitogenome phylogenomics, single batch |
| `mito_only` | Reads → trimming → NOVOPlasty/MITOS2 mitogenome assembly only (single sample batch, no convergence) |
| `per_sample` | Per-sample steps only (assembly through 1-to-1 orthologs + BUSCO + mitogenome assembly). Use this to process independent sequencing/SRA batches separately |
| `downstream` | Multi-sample phylogenomics steps only (nuclear orthologs + mitogenome), combining results from one or more previous `per_sample` runs (see `--previous_outdirs`) |

### Required parameters

| Parameter | Description |
|---|---|
| `--folder <dir>` | Project root folder (see directory structure above) |

### Input mode (choose one)

| Parameter | Description |
|---|---|
| `--download_sra true --sra_list <file>` | Download reads from SRA. CSV in `Scripts/`: `Genus_species_locality_code,SRRcode` |
| `--download_sra false --samp_list <file>` | Use local reads. CSV in `Scripts/`: `sample_id,...`. Files expected at `IlluminaOutput/<sample_id>_R1.fastq.gz` / `_R2.fastq.gz` |

### Databases / references

| Parameter | Default | Description |
|---|---|---|
| `--taxon` | `phyllostomidae` | Prefix for rRNA/mtDNA bowtie2 indexes (`databases/bowtie/<taxon>_rRNA`, `<taxon>_mtDNA`) |
| `--ortholog_taxon` | `Microbat` | Ensembl reference taxon name used in the 1-to-1 ortholog list |
| `--busco_lineage` | `mammalia_odb12` | BUSCO lineage dataset name (auto-downloaded) |

### Mitogenome assembly

| Parameter | Default | Description |
|---|---|---|
| `--novoplasty_seed <fasta>` | *required* | Mitochondrial seed sequence (e.g. COI) |
| `--novoplasty_genome_range` | `12000-22000` | Expected mitogenome size range |
| `--novoplasty_kmer` | `23` | NOVOPlasty k-mer size |
| `--novoplasty_read_length` | `151` | Read length |
| `--novoplasty_insert_size` | `300` | Insert size |
| `--mitos_refdir` | `/dlocal/.../MITOS2` | MITOS2 reference database directory |
| `--mitos_refseqver` | `refseq63m` | MITOS2 reference sequence version |
| `--mitos_code` | `2` | Genetic code |

### Downstream mode only

| Parameter | Description |
|---|---|
| `--previous_outdirs <list>` | Comma-separated list of previous `per_sample` output directories, e.g. `results_batch1,results_batch2` |

### Other

| Parameter | Default | Description |
|---|---|---|
| `--outdir` | `results` | Output directory |
| `--min_species_genetree` | `10` | Minimum number of species required to build a gene tree |
| `--min_seq_guidance` | `6` | Minimum number of sequences required to run an alignment (nuclear or mito gene) |
| `--help` / `--h` | | Show help and exit |

### Examples

```bash
# Whole pipeline, single batch, local reads
nextflow run main.nf --mode full --download_sra false --samp_list List_samples.csv \
    --novoplasty_seed Seed_COI.fasta --folder /path/to/PROJECT

# Per-sample only, first sequencing batch
nextflow run main.nf --mode per_sample --samp_list Batch1_samples.csv \
    --novoplasty_seed Seed_COI.fasta --outdir results_batch1 --folder /path/to/PROJECT

# Per-sample only, SRA batch
nextflow run main.nf --mode per_sample --download_sra true --sra_list Batch2_sra.csv \
    --novoplasty_seed Seed_COI.fasta --outdir results_batch2 --folder /path/to/PROJECT

# Combine batches for phylogenomics (nuclear + mitogenome)
nextflow run main.nf --mode downstream --previous_outdirs results_batch1,results_batch2 \
    --outdir results_final --folder /path/to/PROJECT

# Mitogenome only
nextflow run main.nf --mode mito_only --download_sra false --samp_list List_samples.csv \
    --novoplasty_seed Seed_COI.fasta --folder /path/to/PROJECT
```

---

## Output structure

```
<outdir>/
├── CombinedSRR/
│   └── <sample_id>/
│       ├── sra2/                          # raw downloaded reads (SRA mode)
│       ├── trimgalore/                    # corrected + trimmed reads
│       ├── fastqc/
│       │   ├── post_trimgalore/
│       │   ├── post_decontam/
│       │   └── post_bowtiefilter/
│       ├── decontam/                      # contaminant-filtered reads + mapping log
│       ├── bowtie/                        # rRNA/mtDNA-filtered reads, cleaned_edited.*.fq.gz
│       ├── trinity/                       # <sample_id>.Trinity.fasta
│       ├── rebowtie/                      # rRNA-filtered contigs
│       ├── transrate/                     # good/bad transcripts, scores
│       ├── chimera/                       # chimera-filtered transcripts
│       ├── salmon_corset/                 # clusters, representative transcripts
│       ├── transdecoder/                  # .pep.fa, .cds.fa
│       ├── cd-hit-est/                    # unique transcripts
│       ├── ReciprocalBLAST/               # candidate orthologs
│       ├── orthologs/                     # <sample_id>_<ENSG_gene>_CDS_OnetoOne.fasta (per gene)
│       ├── busco_Summaries/               # BUSCO short summary, full table, figure
│       └── novoplasty/
│           ├── log_<sample_id>.txt                    # NOVOPlasty status (circularized / merged / not circularized), assembly length, warnings
│           ├── Circularized_assembly_1_<sample_id>.fas
│           ├── <sample_id>_whole_mito.fasta           # header renamed to <sample_id>
│           ├── mitos_out/                             # MITOS2 annotation (result.fas, GFF, etc.)
│           └── genes/                                 # <sample_id>_<gene>.fasta, one file per annotated gene/tRNA/rRNA/OH/OL
│
├── multifasta/                            # <gene>_MultiFasta_forGuidance_GeneName.fasta (all samples, per nuclear gene)
├── Guidance-Prank/
│   ├── Alignment_Raw/
│   ├── Alignment_Cleaned/
│   ├── ImageAln/
│   └── MACSE/
├── Genetree/                              # first-pass nuclear gene trees (IQ-TREE2, GTR+I+G)
├── TreeShrink/                            # TreeShrink outputs, MatrixOK/*_MatrixOK.fasta
├── MatrixOK/
│   ├── Concatenated_Matrix.fasta
│   ├── Concatenated_Partition_byCodonPos.nex
│   ├── Partition_byGene.raxml
│   ├── Partition_byGeneCodon.raxml
│   └── Statistic*.tsv
├── Genetree2/                             # second-pass nuclear gene trees (IQ-TREE2, MFP, on TreeShrink matrices)
├── StatsAlgTree/
│   └── <gene>/
│       ├── <gene>_stats_summary.tsv       # saturation, treeness/RCV, LB score, DVMC, evo rate, RCV, outlier taxa count
│       ├── <gene>_long_branch_score.txt
│       ├── <gene>_outlier_taxa.txt
│       ├── <gene>_noPrematureSTOP.fasta
│       └── <gene>_recombination.txt       # PhiPack
├── mito_multifasta/
│   ├── <gene>_MitoMultiFasta.fasta                    # one per mito gene/tRNA/rRNA/OH/OL, all samples
│   ├── MitoGenome_Whole_MultiFasta.fasta              # whole mitogenome, all samples
│   ├── MitoGenome_Whole_aligned.fasta                 # only if MITO_WHOLE_ALIGN is enabled
│   ├── Statistic_MitoMultiFasta.tsv
│   └── Concatenated/
│       ├── MitoGenes_Concatenated.fasta               # PCG + tRNA + rRNA (excl. OH/OL)
│       ├── Partition_byGene.raxml
│       ├── MitoPCG_Concatenated.fasta                 # 13 protein-coding genes only
│       ├── Partition_byGeneCodon.raxml
│       ├── MitoPCG_pos12_Concatenated.fasta            # PCG, codon positions 1+2 only
│       ├── Partition_byGene_PCGpos12.raxml
│       ├── MitoGenes_pos12_Concatenated.fasta          # full matrix, PCG in positions 1+2 + tRNA/rRNA in full
│       ├── Partition_byGene_Fullpos12.raxml
│       └── Statistics_MitoAlignments.tsv
├── all_samples_contam_summary.tsv
├── all_samples_trinity_stats.tsv
├── all_samples_unique_transcripts.tsv
├── all_genes_stats_summary.tsv
└── pipeline_info/                         # Nextflow execution report, timeline, trace, DAG
```

---

## Notes

- Steps needing more than one conda environment switch inside the same process (e.g. `RCORRECT`, `CORSET`, `CHIMERA`, `RCPBLAST`) do so via `conda deactivate` / `conda activate` inside the `script:` block.
- `Trinity` and its bundled Perl utilities (`analyze_blastPlus_topHit_coverage.pl`, `TrinityStats.pl`) are kept in a dedicated `trinity_env`, isolated from the rest of the tools, due to its heavy and fragile dependency chain.
- Steps with fewer than `--min_seq_guidance` sequences (`GUIDANCE`) or fewer than `--min_species_genetree` species (`GENETREE`) are automatically skipped for that gene.
- `GUIDANCE` and `FILTERALIGNMENTS` are generic modules shared between the nuclear ortholog pipeline and the mitogenome pipeline, driven by a `seqType` (`codon` or `nuc`) carried through the channel tuples.
- `MITOGENOMES` handles three NOVOPlasty outcomes: a cleanly circularized single contig, a merged multi-contig assembly (first candidate option is used automatically — check `Merged_contigs_<sample_id>.txt` manually if in doubt), or a raw non-circularized contig (flagged in the log, may represent an incomplete mitogenome, typically missing the control region/D-loop). Assemblies shorter than 14 kb are flagged with a warning.
- The origin-of-replication regions annotated by MITOS2 (`OH`, `OL`) are excluded from `MITO_CONCAT` matrices, as they are not homologous coding/structural genes.
- `MITO_WHOLE_ALIGN` (MAFFT alignment of the entire assembled mitogenome, including the control region) is disabled by default: nothing downstream consumes it, and whole-genome alignment of a circular molecule is sensitive to inconsistent start-position rotation between samples. Enable it in `main.nf` only if you want this file for manual QC/visualization.

## Author

Nicolas Nesi — UMR 1311 Inserm DYNAMICURE, Université de Caen Normandie
