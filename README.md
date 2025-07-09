# Chlorophytes Genome Annotation Pipeline

A comprehensive Nextflow pipeline for chlorophytes genome annotation implementing all the bioinformatics analyses described in the T. rotula genome annotation study.

## Overview

This pipeline implements a complete workflow for chlorophytes genome annotation including:

- **Ploidy Evaluation**: Variant calling with FreeBayes and filtering with VCFtools
- **RNA Families Prediction**: Sequence similarity searches with Infernal against RFAM database
- **Transposable Elements**: Annotation with EDTA software
- **Segmental Duplication**: Identification with Asgart software
- **Transcriptome Assembly**: Genome-guided assembly with STAR and Trinity
- **Gene Annotation**: Comprehensive annotation with Braker2, MAKER2, and supporting tools
- **Functional Annotation**: Annotation with PANNZER2 and KEGG (KAAS)
- **Methylation Analysis**: Nanopore methylation analysis with modkit and mbtools
- **Phylogenomics**: Tree inference with RAxML-NG

## Quick Start

1. **Install Nextflow** (>=22.10.1):
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. **Install dependencies** using conda/mamba:
   ```bash
   conda env create -f environment.yml
   conda activate chlorophytes-annotation
   ```

3. **Test the pipeline**:
   ```bash
   nextflow run . -profile test,conda
   ```

4. **Run with your data**:
   ```bash
   nextflow run . -params-file params.yaml -profile conda
   ```

## Pipeline Steps

### 1. Ploidy Evaluation
- **FreeBayes** (v1.3.4): Variant calling with configurable mapping quality, base quality, and coverage thresholds
- **VCFtools** (v0.1.17): VCF filtering by depth and quality score

### 2. RNA Families Prediction  
- **Infernal** (v1.1.4): Sequence similarity searches against RFAM database
- E-value filtering and removal of partial/fragmented matches

### 3. Transposable Elements
- **EDTA** (v2.0.0): Complete TE annotation pipeline
- Generates softmasked genome and non-redundant TE library

### 4. Segmental Duplication
- **Asgart** (v2.3): Identifies genome segmental duplications
- **BEDTools** (v2.30.0): Evaluates overlap with coding genes

### 5. Transcriptome Assembly
- **STAR** (v2.7.9a): RNA-seq mapping in double-pass mode
- **Trinity** (v2.11.0): Genome-guided transcriptome assembly with super-transcripts
- **CD-HIT-EST**: Clustering of transcriptome sequences
- **TransDecoder** (v5.5.0): ORF identification and prediction

### 6. Gene Annotation
- **Braker2** (v2.1.5): Gene prediction using RNA-seq and protein evidence
- **Exonerate** (v2.4.0): Protein-to-genome alignment
- **GeneMark-EP+**: Gene prediction with fungus branch point model
- **MAKER2**: Comprehensive gene annotation integration
- **AUGUSTUS** (v3.3): Gene prediction

### 7. Quality Assessment
- **BUSCO**: Annotation completeness assessment for Eukaryota and Stramenopiles lineages

### 8. Functional Annotation
- **PANNZER2**: Functional annotation with GO terms
- **KAAS**: KEGG pathway annotation

### 9. Methylation Analysis
- **Minimap2** (v2.28): Nanopore read mapping
- **modkit** (v0.3.0): Methylation calling and analysis
- **mbtools**: Methylation frequency calculation

### 10. Phylogenomics
- **RAxML-NG**: Maximum-likelihood phylogenetic tree inference

## Input Requirements

### Required Files
- `genome_fasta`: Reference genome in FASTA format
- `illumina_reads`: Paired-end Illumina reads for variant calling
- `nanopore_reads`: Nanopore reads with 5mCG modifications
- `rna_reads`: RNA-seq paired-end reads
- `protein_db`: Stramenopiles protein database for homology searches
- `rfam_db`: RFAM covariance model database
- `tpseudnana_cds`: T. pseudonana CDS sequences for EDTA
- `reference_genomes`: Reference genomes for phylogenomic analysis

### File Formats
- FASTA files for genomes and sequences
- FASTQ files for sequencing reads (gzipped supported)
- Covariance model (.cm) for RFAM database

## Configuration

### Resource Requirements
Default resource limits can be adjusted in `nextflow.config`:
```groovy
params {
    max_memory = '128.GB'
    max_cpus = 16
    max_time = '240.h'
}
```

### Tool Parameters
Key parameters can be modified in `params.yaml`:

```yaml
# Variant calling
min_mapping_quality: 30
min_base_quality: 30
min_coverage: 4

# RNA families
evalue_threshold: 0.01

# Transcriptome assembly
max_intron_length: 100000
trinity_max_memory: "50G"

# Methylation analysis
modkit_percentile: 10
```

## Output Structure

```
results/
├── ploidy_evaluation/          # FreeBayes and VCFtools outputs
├── rna_families/              # Infernal predictions
├── repeats/                   # EDTA transposable elements
├── segmental_duplications/    # Asgart duplications
├── transcriptome/            # STAR, Trinity, CD-HIT outputs
├── annotation/               # Gene annotation results
│   ├── braker2/
│   ├── maker2/
│   ├── busco/
│   └── transdecoder/
├── functional_annotation/    # PANNZER2 and KAAS outputs
├── methylation/             # Methylation analysis
├── phylogenomics/           # RAxML-NG tree
└── pipeline_info/           # Execution reports
```

## Execution Profiles

### Available Profiles
- `conda`: Use conda for dependency management
- `docker`: Use Docker containers
- `singularity`: Use Singularity containers
- `test`: Run with minimal test data
- `test_full`: Run with full test dataset

### Example Commands

```bash
# Run with conda
nextflow run . -profile conda -params-file params.yaml

# Run with Docker
nextflow run . -profile docker -params-file params.yaml

# Run with Singularity on HPC
nextflow run . -profile singularity -params-file params.yaml

# Test run
nextflow run . -profile test,conda
```

## Dependencies

### Core Tools
- FreeBayes (1.3.4)
- VCFtools (0.1.17)  
- Infernal (1.1.4)
- EDTA (2.0.0)
- Asgart (2.3)
- BEDTools (2.30.0)
- STAR (2.7.9a)
- Trinity (2.11.0)
- CD-HIT-EST
- TransDecoder (5.5.0)
- DIAMOND (0.9.29)
- Braker2 (2.1.5)
- Exonerate (2.4.0)
- GeneMark-EP+
- MAKER2
- AUGUSTUS (3.3)
- BUSCO (5.4.4)
- Minimap2 (2.28)
- modkit (0.3.0)
- RAxML-NG

### System Requirements
- Nextflow (>=22.10.1)
- Conda/Mamba or Docker/Singularity
- 64-bit Linux or macOS
- Minimum 16 GB RAM (128 GB recommended)
- Minimum 4 CPU cores (16+ recommended)

## Citation

If you use this pipeline, please cite:

1. The original T. rotula study (add reference)
2. nf-core: Philip Ewels et al. Nature Biotechnology 2020
3. Individual tool citations (automatically included in pipeline outputs)

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Support

For questions and support:
- Open an issue on GitHub
- Check the [nf-core documentation](https://nf-co.re/docs/)
- Join the nf-core Slack community

## License

This pipeline is distributed under the MIT License. See `LICENSE` for details.
