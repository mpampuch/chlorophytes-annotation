# Input Data Summary Table

Quick reference table for all input data requirements in the Chlorophytes Genome Annotation Pipeline.

## Primary Input Data Types

| Input Type | Format | Required | Description | Used By Programs |
|------------|--------|----------|-------------|------------------|
| **Reference Genome** | FASTA | ✅ | Assembled genome sequences | FreeBayes, Infernal, EDTA, STAR, Trinity, Braker2, MAKER2, Minimap2 |
| **Illumina Reads** | FASTQ (paired) | ✅ | Short reads for variant calling | FreeBayes → VCFtools |
| **Nanopore Reads** | FASTQ (5mCG) | ✅ | Long reads with methylation | Minimap2 → modkit → mbtools |
| **RNA-seq Reads** | FASTQ (paired) | ✅ | Transcriptome sequencing | STAR → Trinity → Braker2 |
| **Protein Database** | FASTA | ✅ | Stramenopiles proteins | DIAMOND, TransDecoder, Braker2, MAKER2 |
| **RFAM Database** | Covariance Model | ✅ | RNA family models | Infernal |
| **T. pseudonana CDS** | FASTA | ✅ | Related species CDS | EDTA |
| **Reference Genomes** | FASTA | ✅ | Multiple genomes for phylogeny | RAxML-NG |

## Program Input/Output Flow

```
INPUT DATA                    PROGRAMS                    OUTPUT DATA
┌─────────────────┐
│ Genome FASTA    │──────────► FreeBayes ──────────────► VCF variants
│ Illumina Reads  │           VCFtools ──────────────► Filtered VCF
└─────────────────┘

┌─────────────────┐
│ Genome FASTA    │──────────► Infernal ───────────────► RNA families
│ RFAM Database   │
└─────────────────┘

┌─────────────────┐
│ Genome FASTA    │──────────► EDTA ────────────────────► TE annotations
│ T.pseudo CDS    │                    │                  Masked genome
└─────────────────┘                    │
                                       │
┌─────────────────┐                    │
│ Genome FASTA    │──────────► STAR ───┼─────────────────► Gene index
│ RNA-seq Reads   │           Trinity ─┼─────────────────► Transcriptome
└─────────────────┘          CD-HIT ──┼─────────────────► Clustered seqs
                             TransDecoder ────────────────► Predicted ORFs
                                       │
┌─────────────────┐                    │
│ Protein DB      │──────────► DIAMOND ┼─────────────────► Homology hits
└─────────────────┘                    │
                                       │
                             Braker2 ──┼─────────────────► Gene predictions
                             MAKER2 ───┼─────────────────► Final annotation
                             BUSCO ────┼─────────────────► Quality metrics
                                       │
                             PANNZER2 ─┼─────────────────► GO annotations
                             KAAS ─────┼─────────────────► KEGG pathways
                                       │
┌─────────────────┐                    │
│ Nanopore Reads  │──────────► Minimap2┼─────────────────► Aligned reads
└─────────────────┘           modkit ──┼─────────────────► Methylation calls
                              mbtools ─┼─────────────────► Methylation stats
                                       │
┌─────────────────┐                    │
│ Reference       │──────────► RAxML-NG┼─────────────────► Phylogenetic tree
│ Genomes         │                    │
└─────────────────┘                    │
```

## Data Size Requirements

| Data Type | Typical Size | Storage Needs |
|-----------|--------------|---------------|
| Genome Assembly | 50-200 MB | 1-2 GB (with indices) |
| Illumina Reads | 5-50 GB | 50-100 GB (working space) |
| Nanopore Reads | 10-100 GB | 100-200 GB (working space) |
| RNA-seq Reads | 2-20 GB per sample | 20-100 GB total |
| Protein Database | 10-50 MB | 100 MB (with indices) |
| RFAM Database | 500 MB | 2 GB (with indices) |
| Work Directory | - | 500 GB - 2 TB |
| Final Results | - | 50-200 GB |

## Critical Data Quality Thresholds

| Metric | Minimum | Recommended |
|--------|---------|-------------|
| **Genome BUSCO Completeness** | >85% | >95% |
| **Genome N50** | >100 kb | >1 Mb |
| **Illumina Coverage** | 20x | 30-50x |
| **Nanopore Coverage** | 15x | 20-30x |
| **RNA-seq Reads** | 20M per sample | 50-100M per sample |
| **Illumina Q30** | >80% | >90% |

## File Naming Conventions

### Required Patterns

```bash
# Illumina paired-end reads
sample1_R1.fastq.gz  # Forward reads
sample1_R2.fastq.gz  # Reverse reads

# RNA-seq paired-end reads  
condition1_rep1_R1.fastq.gz
condition1_rep1_R2.fastq.gz
condition1_rep2_R1.fastq.gz
condition1_rep2_R2.fastq.gz

# Genome assembly
genome.fasta         # Main assembly
genome.fasta.fai     # FASTA index (auto-generated)

# Nanopore reads
nanopore_5mCG.fastq.gz  # Must include 5mCG modification calls
```

## Essential Database Downloads

```bash
# RFAM Database
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm

# T. pseudonana CDS (example - adjust URL)
wget https://genome.jgi.doe.gov/portal/Thaps3/download/Thaps3_GeneCatalog_CDS_20120524.fasta.gz

# Stramenopiles proteins (example - from UniProt)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_stramenopiles.dat.gz

# BUSCO lineage datasets
busco --download eukaryota_odb10
busco --download stramenopiles_odb10
```

## Quick Data Validation Commands

```bash
# Check FASTA files
grep -c "^>" genome.fasta          # Count sequences
head -n 20 genome.fasta           # Check format

# Check FASTQ files  
zcat reads_R1.fastq.gz | head -n 8  # Check format
fastqc *.fastq.gz                   # Quality control

# Check file sizes
ls -lh *.{fasta,fastq.gz}          # Check all input files

# Validate dependencies
which nextflow singularity          # Check software
```

This summary provides a quick reference for preparing and validating all input data required for the chlorophytes annotation pipeline!