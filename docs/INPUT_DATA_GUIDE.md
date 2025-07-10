# Input Data Requirements Guide

This document details the specific input data requirements for each bioinformatics program in the Chlorophytes Genome Annotation Pipeline.

## Required Input Data Types

### 1. Reference Genome Data

#### **Genome FASTA File** (`genome_fasta`)
- **Format**: FASTA (.fa, .fasta, .fna)
- **Content**: Assembled genome sequences
- **Requirements**: 
  - Single or multi-FASTA format
  - Chromosome/scaffold sequences
  - Headers should be simple (avoid special characters)
- **Used by**: FreeBayes, Infernal, EDTA, Asgart, STAR, Trinity, Braker2, MAKER2, Minimap2

**Example**:
```
>chromosome_1
ATCGATCGATCG...
>chromosome_2
GCTAGCTAGCTA...
```

### 2. Sequencing Read Data

#### **Illumina Paired-End Reads** (`illumina_reads`)
- **Format**: FASTQ (.fastq, .fq, .fastq.gz, .fq.gz)
- **Content**: Short-read sequencing data (typically 100-300bp)
- **Requirements**:
  - Paired-end reads: `*_R1.fastq.gz` and `*_R2.fastq.gz`
  - Quality scores in Phred+33 format
  - Adapter sequences removed (or will be trimmed)
- **Used by**: FreeBayes (for variant calling and ploidy evaluation)

**Example file pattern**:
```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

#### **Nanopore Long Reads** (`nanopore_reads`)
- **Format**: FASTQ (.fastq, .fq, .fastq.gz)
- **Content**: Long-read sequencing data with 5mCG methylation calls
- **Requirements**:
  - Single-end long reads (1kb-100kb+)
  - Must include base modification calls (5mCG)
  - Basecalled with Guppy or Dorado with methylation model
- **Used by**: Minimap2, modkit, mbtools (for methylation analysis)

**Example**:
```
@read_1 modification_calls=5mCG
ATCGATCGATCGATCG...
+
####################
```

#### **RNA-seq Paired-End Reads** (`rna_reads`)
- **Format**: FASTQ (.fastq, .fq, .fastq.gz, .fq.gz)
- **Content**: RNA sequencing data for transcriptome analysis
- **Requirements**:
  - Paired-end reads: `*_R1.fastq.gz` and `*_R2.fastq.gz`
  - Strand-specific or non-strand-specific
  - Multiple biological replicates recommended
- **Used by**: STAR, Trinity, Braker2 (for gene annotation)

### 3. Database and Reference Files

#### **Protein Database** (`protein_db`)
- **Format**: FASTA (.fa, .fasta, .faa)
- **Content**: Protein sequences from related organisms (Stramenopiles)
- **Requirements**:
  - High-quality, well-annotated protein sequences
  - Preferably from closely related species
  - Headers with meaningful identifiers
- **Used by**: DIAMOND BLAST, TransDecoder, Braker2, MAKER2

**Example**:
```
>sp|P12345|PROT1_THAPS Protein description [Thalassiosira pseudonana]
MKLLVVDEADRMLEVLANQAKDVKEALLKQAENLLKDRQS...
>sp|P67890|PROT2_THAPS Another protein [Thalassiosira pseudonana]
MGKQVWVDFSRENLHQFNDSRKELKHMWEEAFKALEQLS...
```

#### **RFAM Database** (`rfam_db`)
- **Format**: Covariance model (.cm)
- **Content**: RNA family database with structural models
- **Requirements**:
  - Current RFAM release (Dec 2021 or later)
  - Covariance model format for Infernal
  - Must be pressed/indexed for cmscan
- **Used by**: Infernal (for RNA family prediction)

**Download**: 
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
```

#### **T. pseudonana CDS Sequences** (`tpseudnana_cds`)
- **Format**: FASTA (.fa, .fasta, .fna)
- **Content**: Coding sequences from Thalassiosira pseudonana
- **Requirements**:
  - Complete CDS sequences (not genomic regions)
  - Used as closely related species for TE annotation
- **Used by**: EDTA (for transposable element detection)

#### **Reference Genomes** (`reference_genomes`)
- **Format**: FASTA (.fa, .fasta)
- **Content**: Multiple genome sequences for phylogenetic analysis
- **Requirements**:
  - Well-annotated genomes from related species
  - Protein sequences or whole genomes
  - Consistent sequence identifiers
- **Used by**: RAxML-NG (for phylogenomic tree construction)

## Program-Specific Input Requirements

### Ploidy Evaluation

#### **FreeBayes**
```yaml
Inputs:
  - Reference genome (FASTA)
  - Illumina paired-end reads (FASTQ)
Parameters:
  - min_mapping_quality: 30
  - min_base_quality: 30  
  - min_coverage: 4
Output:
  - Variant calls (VCF)
```

#### **VCFtools**
```yaml
Inputs:
  - VCF file from FreeBayes
Parameters:
  - min_depth: 10
  - min_quality: 10
Output:
  - Filtered VCF file
```

### RNA Family Prediction

#### **Infernal**
```yaml
Inputs:
  - Reference genome (FASTA)
  - RFAM database (covariance models)
Parameters:
  - evalue_threshold: 0.01
Output:
  - RNA family annotations (tblout)
```

### Transposable Elements

#### **EDTA**
```yaml
Inputs:
  - Reference genome (FASTA)
  - T. pseudonana CDS sequences (FASTA)
Parameters:
  - species: "others"
Output:
  - TE annotations (GFF3)
  - Masked genome (FASTA)
  - TE library (FASTA)
```

### Transcriptome Assembly

#### **STAR**
```yaml
Inputs:
  - Reference genome (FASTA)
  - RNA-seq paired-end reads (FASTQ)
Parameters:
  - sjdbOverhang: 99
  - max_intron_length: 100000
Output:
  - Genome index
  - Aligned reads (BAM)
```

#### **Trinity**
```yaml
Inputs:
  - Aligned RNA-seq reads (BAM)
  - Reference genome (FASTA)
Parameters:
  - max_memory: "50G"
  - CPU: 16
Output:
  - Transcriptome assembly (FASTA)
```

#### **CD-HIT-EST**
```yaml
Inputs:
  - Transcriptome sequences (FASTA)
Parameters:
  - identity: 0.95
  - global_alignment: 1
Output:
  - Clustered sequences (FASTA)
```

#### **TransDecoder**
```yaml
Inputs:
  - Clustered transcripts (FASTA)
  - DIAMOND BLAST results (optional)
Parameters:
  - min_protein_length: 20
Output:
  - Predicted proteins (FASTA)
  - ORF coordinates (GFF3)
```

### Gene Annotation

#### **Braker2**
```yaml
Inputs:
  - Masked genome (FASTA)
  - RNA-seq alignments (BAM)
  - Protein sequences (FASTA)
Parameters:
  - species: "generic"
Output:
  - Gene predictions (GFF3)
  - Protein sequences (FASTA)
```

#### **MAKER2**
```yaml
Inputs:
  - Masked genome (FASTA)
  - Braker2 predictions (GFF3)
  - Protein alignments
Parameters:
  - cpus: 8
Output:
  - Final gene annotation (GFF3)
  - Gene sequences (FASTA)
  - Protein sequences (FASTA)
```

### Quality Assessment

#### **BUSCO**
```yaml
Inputs:
  - Predicted proteins (FASTA)
Parameters:
  - lineage: "eukaryota_odb10"
  - lineage: "stramenopiles_odb10"
  - mode: "genome"
Output:
  - Completeness statistics
```

### Functional Annotation

#### **PANNZER2**
```yaml
Inputs:
  - Protein sequences (FASTA)
Parameters:
  - min_query_coverage: 0.4
  - min_subject_coverage: 0.4
  - min_alignment_length: 50
Output:
  - GO annotations (TSV)
  - Functional descriptions
```

#### **KAAS (KEGG)**
```yaml
Inputs:
  - Protein sequences (FASTA)
Parameters:
  - organisms: "cre,bmi,pti,fcy,tps,ath,olu,mpu"
Output:
  - KEGG pathway annotations
  - KO identifiers
```

### Methylation Analysis

#### **Minimap2**
```yaml
Inputs:
  - Nanopore reads with methylation calls (FASTQ)
  - Reference genome (FASTA)
Parameters:
  - preset: "map-ont"
Output:
  - Aligned reads (BAM)
```

#### **modkit**
```yaml
Inputs:
  - Aligned nanopore reads (BAM)
  - Gene coordinates (BED/GFF)
  - Repeat coordinates (BED/GFF)
Parameters:
  - percentile_threshold: 10
  - min_coverage: 1
Output:
  - Methylation calls (BED)
  - Summary statistics (TSV)
```

### Phylogenomics

#### **RAxML-NG**
```yaml
Inputs:
  - Protein sequences (FASTA)
  - Reference genome proteins (FASTA)
Parameters:
  - model: "LG+G4"
  - bootstrap: 100
Output:
  - Phylogenetic tree (Newick)
```

## Data Quality Requirements

### Genome Assembly
- **Completeness**: >90% BUSCO completeness recommended
- **Contiguity**: N50 >1Mb preferred
- **Size**: Expected genome size for the organism
- **Contamination**: <5% contamination

### Sequencing Data
- **Coverage**: 
  - Illumina: 30-50x for variant calling
  - Nanopore: 20-30x for methylation analysis  
  - RNA-seq: 50-100M reads per sample
- **Quality**: Q30 >90% for Illumina reads
- **Insert Size**: 300-500bp for Illumina paired-end

### Protein Database
- **Completeness**: High-quality, complete proteomes
- **Relevance**: Phylogenetically close species preferred
- **Size**: 10,000-50,000 protein sequences typical

## File Format Specifications

### FASTA Files
```
>sequence_identifier optional_description
ATCGATCGATCGATCG
>another_sequence
GCTAGCTAGCTAGCTA
```

### FASTQ Files
```
@read_identifier
ATCGATCGATCGATCG
+
####################
```

### VCF Files
```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	100	.	A	T	60	PASS	DP=30	GT:DP	1/1:30
```

### GFF3 Files
```
##gff-version 3
chr1	EDTA	transposable_element	1000	2000	.	+	.	ID=TE001;Name=LTR_retrotransposon
```

## Common Data Preparation Steps

### 1. Genome Preparation
```bash
# Remove redundant sequences
# Rename contigs with simple names
# Check for contamination
# Validate FASTA format
```

### 2. Read Quality Control
```bash
# Check read quality with FastQC
# Trim adapters if necessary
# Remove low-quality reads
# Check for contamination
```

### 3. Database Preparation
```bash
# Download latest versions
# Index databases as needed
# Validate file formats
# Check completeness
```

This comprehensive input guide ensures all required data is properly formatted and meets quality standards for successful pipeline execution!