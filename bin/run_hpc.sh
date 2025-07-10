#!/bin/bash

# Chlorophytes Genome Annotation Pipeline - HPC Execution Script
# This script runs Nextflow on a login node, which then submits jobs to SLURM
# DO NOT submit this script to SLURM - run it directly on the login node

set -euo pipefail

# ==============================================
# CONFIGURATION - MODIFY THESE PATHS
# ==============================================

# Pipeline parameters file (modify according to your data)
PARAMS_FILE="params_hpc.yaml"

# Nextflow work directory (use fast storage if available)
WORK_DIR="/scratch/$USER/chlorophytes-work"

# Singularity cache directory (use persistent storage)
SINGULARITY_CACHE="/scratch/$USER/singularity"

# Results output directory
OUTDIR="results"

# Maximum resources (adjust according to your HPC)
MAX_CPUS=256
MAX_MEMORY="1.TB"
MAX_TIME="72.h"

# Email for notifications (optional)
EMAIL="your.email@institution.edu"

# HPC configuration file (modify as needed - use slurm.config or your custom config)
HPC_CONFIG="conf/slurm.config"

# ==============================================
# SETUP
# ==============================================

echo "Starting Chlorophytes Annotation Pipeline at $(date)"
echo "Running on: $(hostname)"
echo "User: $(whoami)"

# Create necessary directories
mkdir -p logs
mkdir -p reports
mkdir -p "$WORK_DIR"
mkdir -p "$SINGULARITY_CACHE"

# Load required modules (uncomment and modify as needed for your HPC)
# module load nextflow/22.10.1
# module load singularity/3.8.0
# module load java/11

# Set environment variables
export NXF_WORK="$WORK_DIR"
export SINGULARITY_CACHEDIR="$SINGULARITY_CACHE"
export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHE"

# Prevent local conda/python environments from interfering
export PYTHONNOUSERSITE=1

# ==============================================
# PIPELINE EXECUTION
# ==============================================

echo "Launching Nextflow pipeline..."
echo "Nextflow will submit individual jobs to SLURM"
echo "Monitor with: squeue -u $(whoami)"

nextflow run . \
    -profile slurm \
    -c "$HPC_CONFIG" \
    -params-file "$PARAMS_FILE" \
    -work-dir "$WORK_DIR" \
    --outdir "$OUTDIR" \
    --max_cpus "$MAX_CPUS" \
    --max_memory "$MAX_MEMORY" \
    --max_time "$MAX_TIME" \
    --email "$EMAIL" \
    -with-report "reports/execution_report_$(date +%Y%m%d_%H%M%S).html" \
    -with-timeline "reports/execution_timeline_$(date +%Y%m%d_%H%M%S).html" \
    -with-trace "reports/execution_trace_$(date +%Y%m%d_%H%M%S).txt" \
    -with-dag "reports/pipeline_dag_$(date +%Y%m%d_%H%M%S).html" \
    -resume

# ==============================================
# CLEANUP AND REPORTING
# ==============================================

exit_code=$?

echo "Pipeline finished at $(date)"
echo "Exit code: $exit_code"

if [ $exit_code -eq 0 ]; then
    echo "✅ Pipeline completed successfully!"
    echo "Results are available in: $OUTDIR"
    echo "Reports are available in: reports/"
else
    echo "❌ Pipeline failed with exit code: $exit_code"
    echo "Check the logs for more details:"
    echo "  - Nextflow log: .nextflow.log"
    echo "  - Individual job logs in: $WORK_DIR/*/.command.log"
fi

# Optional: Clean up work directory if pipeline succeeded and cleanup is desired
# if [ $exit_code -eq 0 ] && [ "$CLEANUP_WORK" = "true" ]; then
#     echo "Cleaning up work directory..."
#     rm -rf "$WORK_DIR"
# fi

exit $exit_code