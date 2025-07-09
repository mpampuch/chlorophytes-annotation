#!/bin/bash

# Example script demonstrating how to run the Nanopore vs Illumina Quality Assessment Pipeline
# This script shows different usage patterns and configurations

set -euo pipefail

echo "=========================================="
echo "Nanopore vs Illumina Quality Assessment"
echo "Pipeline Execution Examples"
echo "=========================================="

# Check if YAK is installed
if ! command -v yak &> /dev/null; then
    echo "ERROR: YAK is not installed or not in PATH"
    echo "Please run: bash install_yak.sh"
    exit 1
fi

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "ERROR: Nextflow is not installed"
    echo "Please install Nextflow from: https://www.nextflow.io/"
    exit 1
fi

echo "Prerequisites check passed!"
echo ""

# Example 1: Basic usage with minimal parameters
echo "Example 1: Basic Analysis"
echo "========================="
echo "Command:"
echo "nextflow run main.nf \\"
echo "  --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \\"
echo "  --nanopore_reads 'data/nanopore/*.fastq.gz'"
echo ""

# Example 2: High-throughput analysis
echo "Example 2: High-Throughput Analysis"
echo "===================================="
echo "Command:"
echo "nextflow run main.nf \\"
echo "  --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \\"
echo "  --nanopore_reads 'data/nanopore/*.fastq.gz' \\"
echo "  --outdir quality_results \\"
echo "  --threads 16 \\"
echo "  --max_memory 64.GB \\"
echo "  --max_cpus 16"
echo ""

# Example 3: Low-memory configuration
echo "Example 3: Low-Memory Configuration"
echo "==================================="
echo "Command:"
echo "nextflow run main.nf \\"
echo "  --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \\"
echo "  --nanopore_reads 'data/nanopore/*.fastq.gz' \\"
echo "  --kmer_size 25 \\"
echo "  --max_memory 16.GB \\"
echo "  --bloom_bits 35"
echo ""

# Example 4: Using parameter file
echo "Example 4: Using Parameter File"
echo "==============================="
echo "Create params.yaml:"
echo "cat > params.yaml << EOF"
echo "illumina_reads: 'data/illumina/*_{1,2}.fastq.gz'"
echo "nanopore_reads: 'data/nanopore/*.fastq.gz'"
echo "outdir: 'results'"
echo "threads: 8"
echo "kmer_size: 31"
echo "EOF"
echo ""
echo "Run command:"
echo "nextflow run main.nf -params-file params.yaml"
echo ""

# Example 5: Cluster execution
echo "Example 5: SLURM Cluster Execution"
echo "==================================="
echo "Command:"
echo "nextflow run main.nf \\"
echo "  -profile slurm \\"
echo "  --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \\"
echo "  --nanopore_reads 'data/nanopore/*.fastq.gz' \\"
echo "  --threads 16"
echo ""

# Example 6: Test run
echo "Example 6: Test Run (Help)"
echo "=========================="
echo "Command:"
echo "nextflow run main.nf --help"
echo ""

# Interactive execution option
echo "=========================================="
echo "Interactive Execution"
echo "=========================================="
echo "Would you like to run a test to check if the pipeline works? (y/n)"
read -r response

if [[ "$response" =~ ^[Yy]$ ]]; then
    echo "Running pipeline help to test installation..."
    nextflow run main.nf --help
    echo ""
    echo "Pipeline test completed successfully!"
    echo ""
    echo "To run with your data, modify the file paths in the examples above"
    echo "and execute the appropriate command."
else
    echo "Skipping test run."
fi

echo ""
echo "=========================================="
echo "Quick Start Guide"
echo "=========================================="
echo "1. Prepare your data:"
echo "   - Illumina paired-end FASTQ files"
echo "   - Nanopore FASTQ files"
echo ""
echo "2. Run the pipeline:"
echo "   nextflow run main.nf \\"
echo "     --illumina_reads 'path/to/illumina/*_{1,2}.fastq.gz' \\"
echo "     --nanopore_reads 'path/to/nanopore/*.fastq.gz'"
echo ""
echo "3. Check results in the 'results' directory"
echo ""
echo "For more information, see the README.md file."
echo ""