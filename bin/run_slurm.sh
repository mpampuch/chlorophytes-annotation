#!/bin/bash
#SBATCH --job-name=chlorophytes-annotation
#SBATCH --output=logs/chlorophytes_%j.out
#SBATCH --error=logs/chlorophytes_%j.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your.email@institution.edu

# Chlorophytes Genome Annotation Pipeline - SLURM Execution Script
# This script submits the Nextflow pipeline to SLURM

set -euo pipefail

# ==============================================
# CONFIGURATION - MODIFY THESE PATHS
# ==============================================

# Pipeline parameters file (modify according to your data)
PARAMS_FILE="params.yaml"

# Nextflow work directory (use fast storage if available)
WORK_DIR="/scratch/$USER/chlorophytes-work"

# Singularity cache directory (use persistent storage)
SINGULARITY_CACHE="/scratch/$USER/singularity"

# Results output directory
OUTDIR="results"

# Maximum resources (adjust according to your HPC)
MAX_CPUS=128
MAX_MEMORY="500.GB"
MAX_TIME="72.h"

# Email for notifications (optional)
EMAIL="your.email@institution.edu"

# ==============================================
# SETUP
# ==============================================

echo "Starting Chlorophytes Annotation Pipeline at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on: $(hostname)"

# Create necessary directories
mkdir -p logs
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

# ==============================================
# PIPELINE EXECUTION
# ==============================================

echo "Launching Nextflow pipeline..."

nextflow run . \
    -profile slurm \
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
    echo "  - SLURM output: logs/chlorophytes_${SLURM_JOB_ID}.out"
    echo "  - SLURM error:  logs/chlorophytes_${SLURM_JOB_ID}.err"
    echo "  - Nextflow log: .nextflow.log"
fi

# Optional: Clean up work directory if pipeline succeeded and cleanup is desired
# if [ $exit_code -eq 0 ] && [ "$CLEANUP_WORK" = "true" ]; then
#     echo "Cleaning up work directory..."
#     rm -rf "$WORK_DIR"
# fi

exit $exit_code