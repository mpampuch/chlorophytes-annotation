# SLURM Executor Example

This document shows a practical example of running the Chlorophytes Annotation Pipeline with Nextflow's SLURM executor.

## Understanding Nextflow SLURM Execution

**Key Concept**: Nextflow runs on the login node and submits individual process jobs to SLURM automatically.

```
Login Node          SLURM Cluster
┌─────────────┐    ┌─────────────────────────────────┐
│ Nextflow    │    │ ┌─────────┐ ┌─────────┐         │
│ Coordinator │───►│ │ Job 1   │ │ Job 2   │  ...    │
│             │    │ │FreeBayes│ │ EDTA    │         │
└─────────────┘    │ └─────────┘ └─────────┘         │
                   └─────────────────────────────────┘
```

## Step-by-Step Example

### 1. Prepare Your Environment

```bash
# Login to your HPC system
ssh user@hpc-login.university.edu

# Navigate to your working directory
cd /project/your_group/chlorophytes-annotation

# Load required modules
module load nextflow/22.10.1
module load singularity/3.8.0

# Check modules are loaded
which nextflow
which singularity
```

### 2. Configure for Your HPC System

```bash
# Copy and edit the HPC configuration template
cp conf/hpc_custom.config conf/my_hpc.config

# Edit the configuration file
nano conf/my_hpc.config
```

Example customization:
```groovy
params {
    hpc_account               = 'my_research_group'     // Your SLURM account
    hpc_partition_normal      = 'compute'               // Your normal partition
    hpc_partition_highmem     = 'bigmem'               // High memory partition
    hpc_modules               = 'nextflow,singularity' // Required modules
    scratch_dir               = '/scratch/$USER'       // Scratch directory
    singularity_cache_dir     = '/project/shared/containers'
}

singularity {
    runOptions = '--cleanenv --containall --bind /scratch --bind /project'
}
```

### 3. Prepare Your Data Parameters

```bash
# Copy and edit the HPC parameters template
cp params_hpc.yaml my_params.yaml

# Edit with your actual data paths
nano my_params.yaml
```

Example data configuration:
```yaml
# Input files - use absolute paths
genome_fasta: "/project/your_group/data/T_rotula_genome.fasta"
illumina_reads: "/project/your_group/data/illumina/*_R{1,2}.fastq.gz"
nanopore_reads: "/project/your_group/data/nanopore/ont_reads.fastq.gz"
rna_reads: "/project/your_group/data/rnaseq/*_R{1,2}.fastq.gz"

# Databases
protein_db: "/project/shared/databases/stramenopiles_proteins.fasta"
rfam_db: "/project/shared/databases/Rfam.cm"

# Output directory
outdir: "/project/your_group/results/chlorophytes_annotation"

# HPC-specific settings
max_memory: "512.GB"
max_cpus: 128
max_time: "48.h"
```

### 4. Test Configuration (Recommended)

```bash
# Test with a small subset first
nextflow run . \
    -profile slurm \
    -c conf/my_hpc.config \
    -params-file my_params.yaml \
    --max_cpus 4 \
    --max_memory 32.GB \
    --max_time 2.h \
    -with-trace test_trace.txt \
    -resume

# Check if jobs are being submitted to SLURM
squeue -u $USER
```

### 5. Run Full Pipeline

```bash
# Start a screen session for long-running jobs
screen -S chlorophytes

# Run the pipeline
nextflow run . \
    -profile slurm \
    -c conf/my_hpc.config \
    -params-file my_params.yaml \
    -work-dir /scratch/$USER/chlorophytes-work \
    -with-report reports/execution_report.html \
    -with-timeline reports/execution_timeline.html \
    -with-trace reports/execution_trace.txt \
    -with-dag reports/pipeline_dag.html \
    -resume

# Detach from screen: Ctrl+A, then D
# Reattach: screen -r chlorophytes
```

### 6. Monitor Progress

```bash
# Check SLURM queue
squeue -u $USER

# Check detailed job info
scontrol show job <JOBID>

# Monitor Nextflow progress
tail -f .nextflow.log

# Check resource usage
sacct -u $USER --format=JobID,JobName,Partition,State,Time,Start,End,NodeList,CPUTime,MaxRSS
```

## Process-Specific SLURM Configuration

Individual processes automatically get appropriate resources:

```groovy
// In your config file
process {
    withName:EDTA {
        queue = 'bigmem'
        cpus = 24
        memory = '200.GB'
        time = '48.h'
        clusterOptions = '--exclusive'
    }
    
    withName:FREEBAYES {
        queue = 'normal' 
        cpus = 6
        memory = '32.GB'
        time = '8.h'
    }
}
```

## Common SLURM Job States

When monitoring with `squeue -u $USER`:

- **PD** (Pending): Job is queued, waiting for resources
- **R** (Running): Job is currently executing
- **CG** (Completing): Job is finishing up
- **CD** (Completed): Job finished successfully

## Troubleshooting

### Job Failures
```bash
# Check failed job details
scontrol show job <JOBID>

# Check Nextflow work directory for specific job
ls work/XX/XXXXXXX/
cat work/XX/XXXXXXX/.command.log
cat work/XX/XXXXXXX/.command.err
```

### Resource Issues
```bash
# If jobs are pending due to resources, check:
sinfo -o "%P %C %m %G %l %N"

# Adjust your configuration accordingly
```

### Container Issues
```bash
# Pre-pull containers to avoid timeouts
singularity pull docker://biocontainers/freebayes:1.3.4--py39h5c33dc1_2
```

## Best Practices

1. **Start Small**: Test with minimal resources first
2. **Use Screen/Tmux**: For long-running pipelines  
3. **Monitor Resources**: Check actual vs requested resources
4. **Clean Work Directory**: Remove after successful runs
5. **Use Resume**: Take advantage of Nextflow's checkpointing
6. **Set Email Notifications**: Get notified of job status

## Example Complete Command

```bash
#!/bin/bash
# Example complete execution script

# Setup
module load nextflow/22.10.1 singularity/3.8.0
export SINGULARITY_CACHEDIR="/project/shared/containers"
export NXF_WORK="/scratch/$USER/chlorophytes-work"

# Run pipeline
nextflow run . \
    -profile slurm \
    -c conf/my_hpc.config \
    -params-file my_params.yaml \
    -work-dir $NXF_WORK \
    --outdir /project/your_group/results \
    --max_cpus 128 \
    --max_memory 512.GB \
    --max_time 48.h \
    --email your.email@university.edu \
    -with-report reports/report_$(date +%Y%m%d).html \
    -with-timeline reports/timeline_$(date +%Y%m%d).html \
    -with-trace reports/trace_$(date +%Y%m%d).txt \
    -resume
```

This approach correctly uses Nextflow's built-in SLURM executor for efficient HPC execution!