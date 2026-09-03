# PhyloTransPip

**Phylotranscriptomics Pipeline** — a Nextflow pipeline from RNA-seq reads to CDSs for phylogenomics and comparative genomics analyses.

From raw Illumina reads or SRA accessions, the pipeline performs read cleaning and decontamination, *de novo* transcriptome assembly with Trinity, chimera and redundancy filtering, ortholog identification against a human reference proteome, multi-sample gene alignment, gene tree inference, and phylogenomic quality control (saturation, treeness/RCV, recombination, etc.). It also builds a mitochondrial genome assembly per sample as an independent side branch.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 25.10
- A SLURM cluster (the pipeline is configured for SLURM by default; adjust `nextflow.config` for other executors)
- `mamba`/`conda`

### Conda/mamba environments

The pipeline expects the following environments to already exist (see `nextflow.config` for the exact activation commands used by each step):

| Environment | Used for |
|---|---|
| `ptp_env` | Core tools: fasterq-dump, clipkit, seqkit, csvtk, bowtie2, cd-hit, blast+, BioPerl, salmon, corset, CIAlign, clipkit, macse, phykit, PhiPack |
| `trinity_env` | Trinity assembly, TrinityStats.pl, analyze_blastPlus_topHit_coverage.pl, blast (recip. blast) |
| `busco_env` | BUSCO |
| `transdecoder_env` | TransDecoder |
| `python2_env` | Legacy Python 2 scripts (chimera detection) |
| `guidance3_env` | Guidance3 (alignment), MACSE, PRANK |
| `phylo_env` | IQ-TREE2, TreeShrink, NOVOPlasty, MITOS2 |

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
| Mitogenome assembly | `MITOGENOMES` | NOVOPlasty + MITOS2 |
| **— multi-sample convergence —** | | |
| Multifasta per gene | `MAKEMULTIFASTA` | concatenates 1-to-1 orthologs across samples |
| Alignment | `GUIDANCE` | MACSE (trim non-homologous) + Guidance3/PRANK |
| Alignment cleaning | `FILTERALIGNMENTS` | MACSE, clipKit, CIAlign |
| Gene tree (pass 1) | `GENETREE` | IQ-TREE2 (≥10 taxa filter) |
| Tree/alignment QC | `QCTREE` | TreeShrink |
| Matrix concatenation | `CONCATENATION` | AMAS (nexus + RAxML partitions, by gene / by gene×codon) |
| Gene tree (pass 2) | `GENETREE2` | IQ-TREE2 on TreeShrink-filtered matrices |
| Alignment/tree statistics | `STATSALGTREE` | PhyKIT (saturation, treeness/RCV, LB score, DVMC, evolutionary rate, outlier taxa, RCV), PhiPack recombination test, internal stop codon check |

---

## Usage

```bash
nextflow run main.nf [options]
```

### Modes (`--mode`)

| Mode | Description |
|---|---|
| `full` (default) | Whole pipeline: per-sample assembly through phylogenomics, single batch |
| `mito_only` | Reads → trimming → NOVOPlasty/MITOS2 mitogenome assembly only |
| `per_sample` | Per-sample steps only (assembly through 1-to-1 orthologs + BUSCO + mitogenome). Use this to process independent sequencing/SRA batches separately |
| `downstream` | Multi-sample phylogenomics steps only, combining orthologs from one or more previous `per_sample` runs (see `--previous_outdirs`) |

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

# Combine batches for phylogenomics
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
│       └── novoplasty/                    # mitogenome assembly + MITOS2 annotation
│
├── multifasta/                            # <gene>_MultiFasta_forGuidance_GeneName.fasta (all samples, per gene)
├── Guidance-Prank/
│   ├── Alignment_Raw/
│   ├── Alignment_Cleaned/
│   ├── ImageAln/
│   └── MACSE/
├── Genetree/                              # first-pass gene trees (IQ-TREE2, GTR+I+G)
├── TreeShrink/                            # TreeShrink outputs, MatrixOK/*_MatrixOK.fasta
├── MatrixOK/
│   ├── Concatenated_Matrix.fasta
│   ├── Concatenated_Partition_byCodonPos.nex
│   ├── Partition_byGene.raxml
│   ├── Partition_byGeneCodon.raxml
│   └── Statistic*.tsv
├── Genetree2/                             # second-pass gene trees (IQ-TREE2, MFP, on TreeShrink matrices)
├── StatsAlgTree/
│   └── <gene>/
│       ├── <gene>_stats_summary.tsv       # saturation, treeness/RCV, LB score, DVMC, evo rate, RCV
│       ├── <gene>_long_branch_score.txt
│       ├── <gene>_outlier_taxa.txt
│       ├── <gene>_noPrematureSTOP.fasta
│       └── <gene>_recombination.txt       # PhiPack
├── all_samples_contam_summary.tsv
├── all_samples_trinity_stats.tsv
├── all_samples_unique_transcripts.tsv
├── all_genes_stats_summary.tsv
└── pipeline_info/                         # Nextflow execution report, timeline, trace, DAG
```

---

## Notes

- Steps needing more than one conda environment switch inside the same process (e.g. `RCORRECT`, `CORSET`, `CHIMERA`) do so via `conda deactivate` / `conda activate` inside the `script:` block, or by activating two environments in sequence in `beforeScript` when both tools are needed simultaneously by the same command.
- `Trinity` and its bundled Perl utilities (`analyze_blastPlus_topHit_coverage.pl`, `TrinityStats.pl`) are kept in a dedicated `trinity_env`, isolated from the rest of the tools, due to its heavy and fragile dependency chain.
- Steps with fewer than 6 sequences (`GUIDANCE`) or fewer than `--min_species_genetree` species (`GENETREE`) are automatically skipped for that gene.

## Author

Nicolas Nesi — UMR 1311 Inserm DYNAMICURE, Université de Caen Normandie
