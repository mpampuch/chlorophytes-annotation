# HPC/SLURM Usage Guide

This guide provides detailed instructions for running the Chlorophytes Genome Annotation Pipeline on High Performance Computing (HPC) systems using SLURM job scheduler and Singularity containers.

**Important**: This pipeline uses Nextflow's built-in SLURM executor. You run Nextflow from a login node, and Nextflow automatically submits individual process jobs to SLURM. Do NOT submit the entire Nextflow run as a SLURM job.

## Quick Start for HPC

```bash
# 1. Clone and setup
git clone <repository>
cd chlorophytes-annotation

# 2. Configure for your HPC system
cp conf/hpc_custom.config conf/my_hpc.config
# Edit conf/my_hpc.config with your HPC settings

# 3. Prepare your data and parameters
cp params_hpc.yaml my_params.yaml
# Edit my_params.yaml with your data paths

# 4. Run on login node (Nextflow will submit jobs to SLURM)
./bin/run_hpc.sh

# 5. Monitor progress
squeue -u $USER
tail -f .nextflow.log
```

## HPC System Requirements

### Software Requirements
- **Nextflow** (≥22.10.1)
- **Singularity** (≥3.7.0)
- **SLURM** job scheduler
- **Java** (≥11)

### Compute Resources
- **Minimum**: 64 CPUs, 256 GB RAM
- **Recommended**: 128+ CPUs, 512+ GB RAM  
- **Storage**: 1-2 TB fast scratch space
- **Runtime**: 24-72 hours (depending on genome size)

## Configuration for HPC

### 1. Basic SLURM Configuration

Use the pre-configured SLURM profile:

```bash
nextflow run . -profile slurm -params-file params.yaml
```

### 2. Custom HPC Configuration

Copy and modify the HPC template:

```bash
cp conf/hpc_custom.config conf/my_hpc.config
```

Edit `conf/my_hpc.config` with your HPC specifications:

```groovy
params {
    hpc_account               = 'myproject'           // Your SLURM account
    hpc_partition_normal      = 'compute'             // Your compute partition
    hpc_partition_highmem     = 'highmem'            // Your high-memory partition
    hpc_modules               = 'nextflow/22.10.1,singularity/3.8.0'
    scratch_dir               = '/scratch/$USER'      // Your scratch directory
    singularity_cache_dir     = '/project/$USER/singularity'
}

singularity {
    runOptions = '--cleanenv --containall --bind /scratch --bind /project'
}
```

### 3. Resource Allocation by Process

The pipeline automatically allocates resources based on computational requirements:

| Process | CPUs | Memory | Time | Queue |
|---------|------|--------|------|-------|
| EDTA | 24 | 200 GB | 72h | highmem |
| TRINITY | 32 | 256 GB | 48h | highmem |
| MAKER2 | 24 | 200 GB | 72h | highmem |
| BRAKER2 | 16 | 128 GB | 48h | highmem |
| STAR | 12 | 48 GB | 12h | normal |
| FreeBayes | 6 | 36 GB | 8h | normal |

## How Nextflow SLURM Executor Works

Nextflow's SLURM executor works differently from traditional batch job submission:

1. **Main Process**: Nextflow runs on the login node as a lightweight coordinator
2. **Job Submission**: Each pipeline process becomes a separate SLURM job
3. **Dependency Management**: Nextflow handles job dependencies automatically
4. **Resource Management**: Each process can request different resources

### Method 1: Using the Provided Script (Recommended)

```bash
# Edit bin/run_hpc.sh with your settings
./bin/run_hpc.sh
```

### Method 2: Direct Command Line

```bash
# Load modules on login node
module load nextflow/22.10.1 singularity/3.8.0

# Run from login node (Nextflow submits jobs to SLURM)
nextflow run . \
    -profile slurm \
    -c conf/slurm.config \
    -params-file params_hpc.yaml \
    --max_cpus 256 \
    --max_memory 1.TB \
    --max_time 72.h
```

### Method 3: Custom Configuration

Create your own SLURM config:

```groovy
// my_slurm.config
executor {
    name = 'slurm'
    queueSize = 100
}

process {
    executor = 'slurm'
    queue = 'your_partition'
    clusterOptions = '--account=your_account'
    
    withLabel:big_job {
        queue = 'bigmem'
        cpus = 32
        memory = '256.GB'
        time = '48.h'
    }
}
```

### Method 4: Screen/Tmux Session (Long-running)

For very long pipelines, run in a persistent session:

```bash
# Start screen session
screen -S chlorophytes

# Load modules and run
module load nextflow/22.10.1 singularity/3.8.0
./bin/run_hpc.sh

# Detach: Ctrl+A, D
# Reattach: screen -r chlorophytes
```

## Singularity Configuration

### Container Cache Management

Set up persistent container cache:

```bash
# Create cache directory
mkdir -p /project/$USER/singularity

# Set environment variables
export SINGULARITY_CACHEDIR=/project/$USER/singularity
export NXF_SINGULARITY_CACHEDIR=/project/$USER/singularity
```

### Pre-pulling Containers (Recommended)

Pre-download containers to avoid network timeouts:

```bash
# Create container list
cat > containers.txt << 'EOF'
https://depot.galaxyproject.org/singularity/freebayes:1.3.4--py39h5c33dc1_2
https://depot.galaxyproject.org/singularity/vcftools:0.1.17--pl5321hdfd78af_4
https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0
https://depot.galaxyproject.org/singularity/trinity:2.11.0--h00214ad_1
https://depot.galaxyproject.org/singularity/busco:5.4.4--pyhdfd78af_0
EOF

# Pre-pull containers
while read container; do
    singularity pull $container
done < containers.txt
```

## Storage Considerations

### Work Directory
- Use fast scratch storage: `/scratch/$USER/work`
- Ensure sufficient space (500GB - 2TB)
- Clean up after successful runs

### Input Data
- Place data on accessible storage
- Use absolute paths in parameters
- Consider data transfer time

### Output Data
- Results directory: persistent storage
- Backup important results
- Use compression for large files

## Monitoring and Troubleshooting

### Monitoring Jobs

```bash
# Check SLURM queue
squeue -u $USER

# Monitor specific job
scontrol show job <JOBID>

# View resource usage
sacct -j <JOBID> --format=JobID,JobName,Partition,Account,AllocCPUS,State,ExitCode,Start,End,Elapsed,MaxRSS,MaxVMSize

# Check Nextflow progress
tail -f .nextflow.log
```

### Common Issues and Solutions

#### 1. Out of Memory Errors
```bash
# Increase memory limits in config
withName:PROBLEMATIC_PROCESS {
    memory = { check_max( 512.GB * task.attempt, 'memory' ) }
}
```

#### 2. Time Limit Exceeded
```bash
# Increase time limits
withName:SLOW_PROCESS {
    time = { check_max( 120.h * task.attempt, 'time' ) }
}
```

#### 3. Container Pull Failures
```bash
# Pre-pull containers or use local cache
singularity.cacheDir = '/project/$USER/singularity'
```

#### 4. SLURM Account Issues
```bash
# Add account to cluster options
clusterOptions = '--account=myproject --qos=normal'
```

### Log Files and Debugging

Important log locations:
- SLURM logs: `logs/chlorophytes_<JOBID>.{out,err}`
- Nextflow log: `.nextflow.log`
- Process logs: `work/*/` directories
- Pipeline reports: `results/pipeline_info/`

## Performance Optimization

### 1. Resource Scaling
Adjust resources based on your data size:

```groovy
// For large genomes (>1GB)
params {
    max_memory = '1.TB'
    max_cpus = 256
    max_time = '120.h'
}

// EDTA configuration for large genomes
withName:EDTA {
    memory = '512.GB'
    time = '120.h'
    cpus = 32
}
```

### 2. Parallel Execution
Optimize parallelization:

```bash
# Allow more concurrent jobs
executor.queueSize = 50
executor.submitRateLimit = '10 sec'
```

### 3. Storage Optimization
```bash
# Use local scratch for work directory
export NXF_WORK=/scratch/$USER/work

# Enable automatic cleanup
cleanup = true
```

## Example HPC Configurations

### SLAC/SLURM Example
```groovy
params {
    hpc_account = 'shared'
    hpc_partition_normal = 'shared'
    hpc_partition_highmem = 'bigmem'
    scratch_dir = '/scratch/$USER'
}
```

### TACC Stampede2 Example
```groovy
params {
    hpc_account = 'TG-MCB140063'
    hpc_partition_normal = 'normal'
    hpc_partition_highmem = 'large-shared'
    scratch_dir = '$SCRATCH'
    hpc_modules = 'nextflow,singularity'
}
```

### NERSC Cori Example
```groovy
params {
    hpc_account = 'm1234'
    hpc_partition_normal = 'regular'
    hpc_partition_highmem = 'bigmem'
    scratch_dir = '$SCRATCH'
    hpc_qos = 'normal'
}
```

## Best Practices

1. **Test First**: Run with small datasets using test profile
2. **Monitor Resources**: Check resource usage and adjust as needed
3. **Use Resume**: Take advantage of Nextflow's resume feature
4. **Backup Data**: Keep copies of important input data
5. **Document Setup**: Keep notes on your HPC configuration
6. **Clean Up**: Remove work directories after successful runs
7. **Share Configs**: Share working configurations with your team

## Support

For HPC-specific issues:
1. Check your HPC documentation
2. Contact your HPC support team
3. Open an issue on GitHub with your configuration
4. Join the nf-core Slack community

## Additional Resources

- [Nextflow on HPC](https://www.nextflow.io/docs/latest/executor.html#slurm)
- [Singularity User Guide](https://sylabs.io/guides/latest/user-guide/)
- [SLURM Documentation](https://slurm.schedmd.com/documentation.html)
- [nf-core HPC Configs](https://github.com/nf-core/configs)