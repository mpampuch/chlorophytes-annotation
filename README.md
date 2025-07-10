# Nanopore vs Illumina Read Quality Assessment Pipeline

A Nextflow pipeline for assessing the quality of nanopore reads compared to illumina reads using k-mer analysis with the [YAK](https://github.com/lh3/yak) tool.

## Overview

This pipeline uses k-mer spectrum analysis to compare the quality of nanopore long reads against illumina short reads. The YAK tool enables reference-free quality assessment by analyzing k-mer frequencies and distributions between different sequencing technologies.

### Key Features

- **Reference-free quality assessment**: No reference genome required
- **K-mer spectrum comparison**: Compares k-mer distributions between illumina and nanopore reads
- **Quality value calculation**: Computes QV (Quality Values) for nanopore reads using illumina k-mers as reference
- **Comprehensive reporting**: Generates HTML and text summary reports
- **Scalable execution**: Supports local, HPC cluster, and containerized execution

### Pipeline Workflow

1. **K-mer counting**: Count k-mers from illumina (reference) and nanopore reads
2. **Spectrum comparison**: Compare k-mer frequency distributions
3. **Quality assessment**: Calculate quality values for nanopore reads
4. **Report generation**: Create comprehensive summary reports

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥ 22.10.0)
- [YAK](https://github.com/lh3/yak) k-mer analyzer
- Python 3.6+ (for report generation)

### Install YAK

```bash
# Clone and compile YAK
git clone https://github.com/lh3/yak.git
cd yak
make

# Add to PATH
export PATH=$PATH:/path/to/yak
```

### Download Pipeline

```bash
# Clone the pipeline repository
git clone <repository-url>
cd nanopore-illumina-quality-assessment

# Test the pipeline
nextflow run main.nf --help
```

## Usage

### Quick Start

```bash
nextflow run main.nf \
  --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \
  --nanopore_reads 'data/nanopore/*.fastq.gz' \
  --outdir results \
  --threads 8
```

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--illumina_reads` | Path to paired-end Illumina FASTQ files (use glob pattern) |
| `--nanopore_reads` | Path to nanopore FASTQ files (use glob pattern) |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--threads` | `4` | Number of threads |
| `--kmer_size` | `31` | K-mer size for analysis |
| `--bloom_bits` | `37` | Bloom filter bits for yak count |
| `--max_memory` | `32.GB` | Maximum memory usage |
| `--max_cpus` | `8` | Maximum CPU cores |
| `--max_time` | `48.h` | Maximum execution time |

### Example Commands

#### Basic analysis
```bash
nextflow run main.nf \
  --illumina_reads 'data/sr*_{1,2}.fastq.gz' \
  --nanopore_reads 'data/ont*.fastq.gz'
```

#### High-throughput analysis with custom parameters
```bash
nextflow run main.nf \
  --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \
  --nanopore_reads 'data/nanopore/*.fastq.gz' \
  --outdir quality_results \
  --threads 16 \
  --kmer_size 25 \
  --max_memory 64.GB
```

#### Using custom configuration
```bash
nextflow run main.nf \
  -c custom.config \
  --illumina_reads 'data/illumina/*.fastq.gz' \
  --nanopore_reads 'data/nanopore/*.fastq.gz'
```

## Execution Profiles

The pipeline supports multiple execution profiles:

### Local Execution
```bash
nextflow run main.nf -profile local [options]
```

### SLURM Cluster
```bash
nextflow run main.nf -profile slurm [options]
```

### Docker
```bash
nextflow run main.nf -profile docker [options]
```

### Singularity
```bash
nextflow run main.nf -profile singularity [options]
```

### Test Profile (minimal resources)
```bash
nextflow run main.nf -profile test [options]
```

## Input Data Format

### Illumina Reads
- Paired-end FASTQ files (R1 and R2)
- Compressed (.gz) or uncompressed formats supported
- File naming convention: `*_{1,2}.fastq.gz` or `*_R{1,2}.fastq.gz`

### Nanopore Reads  
- Single-end FASTQ files
- Compressed (.gz) or uncompressed formats supported
- Multiple files per sample supported

### Example Directory Structure
```
data/
├── illumina/
│   ├── sample1_1.fastq.gz
│   ├── sample1_2.fastq.gz
│   ├── sample2_1.fastq.gz
│   └── sample2_2.fastq.gz
└── nanopore/
    ├── sample1_nanopore.fastq.gz
    ├── sample2_nanopore.fastq.gz
    └── sample3_nanopore.fastq.gz
```

## Output Files

The pipeline generates the following output structure:

```
results/
├── kmer_counts/
│   ├── *.yak                    # K-mer count databases
│   └── *.hist                   # K-mer frequency histograms
├── comparisons/
│   ├── *_comparison.txt         # Detailed comparison reports
│   └── *_kqv.txt               # K-mer quality values
├── quality_assessment/
│   └── *_qv_assessment.txt     # Quality value assessments
└── quality_assessment_summary.html  # Main summary report
```

### Key Output Files

#### 1. K-mer Count Files (`.yak`)
Binary databases containing k-mer frequency tables for each sample.

#### 2. Frequency Histograms (`.hist`)
Text files showing the distribution of k-mer frequencies.

#### 3. Comparison Reports (`*_comparison.txt`)
Detailed comparisons between illumina and nanopore k-mer spectra.

#### 4. Quality Assessments (`*_qv_assessment.txt`)
Quality value calculations for nanopore reads using illumina k-mers as reference.

#### 5. Summary Report (`quality_assessment_summary.html`)
Comprehensive HTML report with analysis overview and interpretation guide.

## Interpreting Results

### Quality Values (QV)
- **QV30**: 99.9% accuracy (1 error per 1,000 bases)
- **QV40**: 99.99% accuracy (1 error per 10,000 bases)
- **QV50**: 99.999% accuracy (1 error per 100,000 bases)

Higher QV values indicate better read quality.

### K-mer Spectra Analysis
- **Peak alignment**: Well-aligned peaks between illumina and nanopore indicate good quality
- **Coverage depth**: Higher coverage in expected ranges suggests good sequencing depth
- **Error patterns**: Unusual distributions may indicate systematic errors or contamination

### Quality Assessment Guidelines
1. **QV > 30**: Excellent quality, suitable for most applications
2. **QV 20-30**: Good quality, acceptable for many analyses
3. **QV < 20**: Poor quality, may require additional filtering or re-sequencing

## Troubleshooting

### Common Issues

#### 1. YAK not found
```
ERROR: yak is not installed or not in PATH
```
**Solution**: Install YAK and ensure it's in your PATH

#### 2. Memory errors
```
Out of memory error during k-mer counting
```
**Solution**: Increase memory allocation or reduce k-mer size
```bash
--max_memory 64.GB --kmer_size 25
```

#### 3. Input file issues
```
No such file or directory
```
**Solution**: Check file paths and glob patterns
```bash
# Use absolute paths if needed
--illumina_reads '/full/path/to/illumina/*_{1,2}.fastq.gz'
```

#### 4. Long execution times
**Solution**: Increase CPU allocation and check resource usage
```bash
--threads 16 --max_cpus 16
```

### Performance Optimization

1. **Adjust k-mer size**: Smaller k-mers (21-25) run faster but may be less specific
2. **Increase threads**: Use more CPU cores for parallel processing
3. **Optimize memory**: Allocate sufficient memory to avoid swapping
4. **Use SSD storage**: Faster I/O improves performance significantly

## Advanced Usage

### Custom Configuration

Create a custom configuration file:

```nextflow
// custom.config
params {
    kmer_size = 25
    threads = 16
    max_memory = '64.GB'
}

process {
    withName: COUNT_ILLUMINA_KMERS {
        memory = '32.GB'
        time = '6.h'
    }
}
```

Run with custom config:
```bash
nextflow run main.nf -c custom.config [other options]
```

### Batch Processing

For multiple sample pairs, create a sample sheet and use parameter files:

```bash
# params.yaml
illumina_reads: 'batch1/illumina/*_{1,2}.fastq.gz'
nanopore_reads: 'batch1/nanopore/*.fastq.gz'
outdir: 'batch1_results'
threads: 8

# Run pipeline
nextflow run main.nf -params-file params.yaml
```

## Citation

If you use this pipeline in your research, please cite:

1. **YAK**: Li, H. (2020). Yak: yet another k-mer analyzer. Available at: https://github.com/lh3/yak
2. **Nextflow**: Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35, 316–319 (2017).

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License

This pipeline is distributed under the MIT License. See `LICENSE` file for details.

## Support

For questions or issues:
1. Check the troubleshooting section above
2. Review YAK documentation: https://github.com/lh3/yak
3. Submit an issue on the project repository

## Changelog

### Version 1.0.0
- Initial release
- Basic k-mer analysis workflow
- Quality assessment functionality
- HTML and text reporting