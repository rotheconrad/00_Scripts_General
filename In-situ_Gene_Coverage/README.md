# Workflow to calculate ANIr and sequence coverage (depth and breadth) of genome(s) / MAG(s) from metagenomes by gene, intergenic region, contig, and whole genome.

This workflow produces separate files in tab separated value (tsv) format for ANIr, sequence depth, and sequence breadth for the genes, intergenic regions, and contigs of a genome / MAG in fasta format. It also produces a file containing sequence depth at each position as well as a file with results calculated for the whole genome sequence. tsv files can be easily opened in Excel, imported into Python with Pandas, or read into R for further analysis. There is also an option to generate some summary plots.


#### Coverage calculated as Truncated Average Depth (TAD):
- Set TAD to 100 for no truncatation.
- TAD 80 removes the top 10% and bottom 10% of base pair depths and caluclates coverage from the middle 80% of values. Intended to reduce effects of conserved motif peaks and contig edge valleys.
- Coverage = base pairs recruited / length of genome, contig, or gene


#### Coverage calculated as Breadth:
- number of positions in reference sequence covered by at least one read alignment divided the length of the reference sequence.


#### Relative Abundance is calculated as:
- base pairs recruited / base pairs in metagenome * 100
- It is the percent of base pairs recruited out of the total base pairs sequenced in the metagenome.


#### ANIr is calculated as:
- average percent identity of sequence alignments for all reads (should be 1 blast match per read)


#### This workflow leads to the following result files:

- 3 column tsv output of Contig(or gene_name), coverage(or ANIr), sequence length.
- Writes 11 files total:
    - {out_file_prefix}_genome_by_bp.tsv
    - {out_file_prefix}_genome.tsv
    - {out_file_prefix}_contig_tad.tsv
    - {out_file_prefix}_contig_breadth.tsv
    - {out_file_prefix}_contig_anir.tsv
    - {out_file_prefix}_gene_tad.tsv
    - {out_file_prefix}_gene_breadth.tsv
    - {out_file_prefix}_gene_anir.tsv
    - {out_file_prefix}_intergene_tad.tsv
    - {out_file_prefix}_intergene_breadth.tsv
    - {out_file_prefix}_intergene_anir.tsv


*This workflow can also be used with Genomic FASTA and CDS from genomic FASTA files retrieved from the [NCBI assembly database](https://www.ncbi.nlm.nih.gov/assembly/). In this case, skip the renaming step for sequence names in the reference fasta file and skip Step 02. Use the -n flag for NCBI in Step 03.*

## Step 00: Required tools :: Blast+ (Blastp) and Kofamscan.


### Prodigal for protein coding gene prediction.
 
Information and installation instructions for Prodigal can be found [here](https://github.com/hyattpd/Prodigal). The publication is [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/).

Prodigal can also be installed with a [conda environment]():

```bash
conda create -n prodigal
conda activate prodigal
conda install -c bioconda prodigal
```


### Magic Blast for short read metagenome mapping to reference genome(s) / MAG(s).

Information and installation instructions for Magic Blast can be found [here](https://ncbi.github.io/magicblast/). The publication is [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2996-x). NOTE: If the latest version of Magic Blast gives errors on your system, navigate to the parent directory try a previous version.


## Step 01: Map metagenomic reads to reference genome(s) / MAG(s).

### Check metagenome read names and rename if needed. (fastq or fasta).

### Check sequence names in reference fasta files and rename if needed.

### For individual genome or MAG:

### Competitive read recruitment for multiple genomes or MAGs:


## Step 02: Predict protein coding genes with Prodigal.


## Step 03: Calculate ANIr and Coverage.


## Step 04: Generate summary plots.
