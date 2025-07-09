#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    Nanopore vs Illumina Read Quality Assessment using YAK
========================================================================================
    Pipeline for assessing the quality of nanopore reads compared to illumina reads
    using k-mer analysis with the yak tool.
    
    Author: AI Assistant
    Version: 1.0.0
========================================================================================
*/

// Parameter definitions
params.illumina_reads = null
params.nanopore_reads = null
params.outdir = 'results'
params.threads = 4
params.kmer_size = 31
params.bloom_bits = 37
params.help = false

// Help message
def helpMessage() {
    log.info"""
    ===================================
    Nanopore vs Illumina Quality Assessment Pipeline
    ===================================
    
    Usage:
    nextflow run main.nf --illumina_reads 'path/to/illumina/*_{1,2}.fastq.gz' --nanopore_reads 'path/to/nanopore/*.fastq.gz'
    
    Required Parameters:
    --illumina_reads      Path to paired-end Illumina FASTQ files (use glob pattern)
    --nanopore_reads      Path to nanopore FASTQ files (use glob pattern)
    
    Optional Parameters:
    --outdir              Output directory (default: 'results')
    --threads             Number of threads (default: 4)
    --kmer_size           K-mer size for analysis (default: 31)
    --bloom_bits          Bloom filter bits for yak count (default: 37)
    --help                Show this help message
    
    Example:
    nextflow run main.nf \\
        --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \\
        --nanopore_reads 'data/nanopore/*.fastq.gz' \\
        --outdir results \\
        --threads 8
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.illumina_reads || !params.nanopore_reads) {
    log.error "Please provide both --illumina_reads and --nanopore_reads parameters"
    helpMessage()
    exit 1
}

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process CHECK_YAK {
    tag "check_yak"
    
    output:
    val true, emit: ready
    
    script:
    """
    if ! command -v yak &> /dev/null; then
        echo "ERROR: yak is not installed or not in PATH"
        echo "Please install yak from: https://github.com/lh3/yak"
        exit 1
    fi
    
    echo "yak version:"
    yak 2>&1 | head -3 || echo "yak found in PATH"
    """
}

process COUNT_ILLUMINA_KMERS {
    tag "illumina_kmers"
    publishDir "${params.outdir}/kmer_counts", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    val ready
    
    output:
    tuple val(sample_id), path("${sample_id}_illumina.yak"), emit: illumina_kmers
    path "${sample_id}_illumina.hist", emit: illumina_hist
    
    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    # Count k-mers from paired-end Illumina reads
    # Using process substitution to provide two identical streams for paired-end data
    yak count \\
        -k ${params.kmer_size} \\
        -b ${params.bloom_bits} \\
        -t ${params.threads} \\
        -o ${sample_id}_illumina.yak \\
        <(zcat ${read1} ${read2}) \\
        <(zcat ${read1} ${read2})
    
    # Generate histogram of k-mer frequencies
    yak inspect ${sample_id}_illumina.yak > ${sample_id}_illumina.hist
    """
}

process COUNT_NANOPORE_KMERS {
    tag "nanopore_kmers_${sample_id}"
    publishDir "${params.outdir}/kmer_counts", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    val ready
    
    output:
    tuple val(sample_id), path("${sample_id}_nanopore.yak"), emit: nanopore_kmers
    path "${sample_id}_nanopore.hist", emit: nanopore_hist
    
    script:
    """
    # Count k-mers from nanopore reads
    yak count \\
        -k ${params.kmer_size} \\
        -K1.5g \\
        -t ${params.threads} \\
        -o ${sample_id}_nanopore.yak \\
        <(zcat ${reads})
    
    # Generate histogram of k-mer frequencies
    yak inspect ${sample_id}_nanopore.yak > ${sample_id}_nanopore.hist
    """
}

process COMPARE_KMER_SPECTRA {
    tag "compare_${illumina_sample}_vs_${nanopore_sample}"
    publishDir "${params.outdir}/comparisons", mode: 'copy'
    
    input:
    tuple val(illumina_sample), path(illumina_yak), val(nanopore_sample), path(nanopore_yak)
    
    output:
    path "${illumina_sample}_vs_${nanopore_sample}_comparison.txt", emit: comparison
    path "${illumina_sample}_vs_${nanopore_sample}_kqv.txt", emit: kmer_qv
    
    script:
    """
    # Compare k-mer spectra between illumina and nanopore
    echo "K-mer spectrum comparison: ${illumina_sample} (Illumina) vs ${nanopore_sample} (Nanopore)" > ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    echo "Generated on: \$(date)" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    echo "Parameters: k=${params.kmer_size}, threads=${params.threads}" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    echo "" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    
    # Generate k-mer quality values
    yak inspect ${illumina_yak} ${nanopore_yak} > ${illumina_sample}_vs_${nanopore_sample}_kqv.txt
    
    # Add summary statistics to comparison file
    echo "=== K-mer Quality Assessment ===" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    echo "Illumina reference: ${illumina_sample}" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    echo "Nanopore sample: ${nanopore_sample}" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    echo "" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    
    # Extract key statistics
    echo "K-mer statistics (first 20 lines):" >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    head -20 ${illumina_sample}_vs_${nanopore_sample}_kqv.txt >> ${illumina_sample}_vs_${nanopore_sample}_comparison.txt
    """
}

process COMPUTE_READ_QUALITY {
    tag "quality_${illumina_sample}_vs_${nanopore_sample}"
    publishDir "${params.outdir}/quality_assessment", mode: 'copy'
    
    input:
    tuple val(illumina_sample), path(illumina_yak), val(nanopore_sample), path(nanopore_reads)
    
    output:
    path "${nanopore_sample}_qv_assessment.txt", emit: quality_report
    
    script:
    """
    # Compute quality values for nanopore reads using illumina k-mers as reference
    yak qv \\
        -t ${params.threads} \\
        -p \\
        -K3.2g \\
        -l100k \\
        ${illumina_yak} \\
        ${nanopore_reads} > ${nanopore_sample}_qv_assessment.txt
    
    # Add header information
    sed -i '1i\\# Quality assessment of ${nanopore_sample} nanopore reads using ${illumina_sample} illumina reference' ${nanopore_sample}_qv_assessment.txt
    sed -i '2i\\# Generated on: '"\$(date)" ${nanopore_sample}_qv_assessment.txt
    sed -i '3i\\# Parameters: k=${params.kmer_size}, threads=${params.threads}' ${nanopore_sample}_qv_assessment.txt
    sed -i '4i\\#' ${nanopore_sample}_qv_assessment.txt
    """
}

process GENERATE_SUMMARY_REPORT {
    tag "summary_report"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path illumina_hists
    path nanopore_hists
    path comparisons
    path quality_reports
    
    output:
    path "quality_assessment_summary.html", emit: summary_html
    path "quality_assessment_summary.txt", emit: summary_txt
    
    script:
    """
    #!/usr/bin/env python3
    
    import os
    import glob
    from datetime import datetime
    
    # Generate text summary
    with open('quality_assessment_summary.txt', 'w') as f:
        f.write("Nanopore vs Illumina Read Quality Assessment Summary\\n")
        f.write("=" * 55 + "\\n")
        f.write(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
        f.write(f"Parameters:\\n")
        f.write(f"  - K-mer size: ${params.kmer_size}\\n")
        f.write(f"  - Threads: ${params.threads}\\n")
        f.write(f"  - Bloom filter bits: ${params.bloom_bits}\\n\\n")
        
        # List processed files
        illumina_files = glob.glob('*illumina.hist')
        nanopore_files = glob.glob('*nanopore.hist')
        comparison_files = glob.glob('*comparison.txt')
        quality_files = glob.glob('*qv_assessment.txt')
        
        f.write(f"Processed Files:\\n")
        f.write(f"  - Illumina samples: {len(illumina_files)}\\n")
        f.write(f"  - Nanopore samples: {len(nanopore_files)}\\n")
        f.write(f"  - Comparisons generated: {len(comparison_files)}\\n")
        f.write(f"  - Quality assessments: {len(quality_files)}\\n\\n")
        
        f.write("Files located in:\\n")
        f.write("  - K-mer counts: ${params.outdir}/kmer_counts/\\n")
        f.write("  - Comparisons: ${params.outdir}/comparisons/\\n")
        f.write("  - Quality assessments: ${params.outdir}/quality_assessment/\\n")
    
    # Generate HTML summary
    with open('quality_assessment_summary.html', 'w') as f:
        f.write('''<!DOCTYPE html>
<html>
<head>
    <title>Nanopore vs Illumina Quality Assessment</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; }
        .highlight { background-color: #f39c12; color: white; padding: 5px; }
        .info-box { background-color: #ecf0f1; padding: 15px; border-left: 4px solid #3498db; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>Nanopore vs Illumina Read Quality Assessment</h1>
    <div class="info-box">
        <p><strong>Analysis completed:</strong> ''' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '''</p>
        <p><strong>Pipeline version:</strong> 1.0.0</p>
        <p><strong>K-mer size:</strong> ${params.kmer_size}</p>
        <p><strong>Threads used:</strong> ${params.threads}</p>
    </div>
    
    <h2>Analysis Overview</h2>
    <p>This pipeline compared nanopore and illumina sequencing data using k-mer analysis with the YAK tool.</p>
    
    <h2>Output Files</h2>
    <ul>
        <li><strong>K-mer counts:</strong> .yak files containing k-mer frequency tables</li>
        <li><strong>Histograms:</strong> .hist files showing k-mer frequency distributions</li>
        <li><strong>Comparisons:</strong> Comparative analysis between illumina and nanopore k-mer spectra</li>
        <li><strong>Quality assessments:</strong> QV (Quality Value) calculations for nanopore reads</li>
    </ul>
    
    <h2>Interpretation Guide</h2>
    <div class="info-box">
        <p><strong>K-mer Quality Values (QV):</strong> Higher QV values indicate better read quality. 
        QV30 means 99.9% accuracy, QV40 means 99.99% accuracy.</p>
        <p><strong>K-mer spectra comparison:</strong> Shows the overlap and differences in k-mer distributions 
        between illumina (reference) and nanopore reads.</p>
    </div>
    
    <h2>Next Steps</h2>
    <ol>
        <li>Review the quality assessment files in <code>${params.outdir}/quality_assessment/</code></li>
        <li>Examine k-mer frequency histograms for data quality indicators</li>
        <li>Compare QV values across different nanopore samples</li>
        <li>Use results to inform downstream analysis decisions</li>
    </ol>
</body>
</html>''')
    """
}

/*
========================================================================================
    WORKFLOWS
========================================================================================
*/

workflow {
    // Check if yak is available
    CHECK_YAK()
    
    // Parse input files
    Channel
        .fromFilePairs(params.illumina_reads, checkIfExists: true)
        .set { illumina_ch }
    
    Channel
        .fromPath(params.nanopore_reads, checkIfExists: true)
        .map { file -> 
            def sample_id = file.getBaseName().replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
            tuple(sample_id, file)
        }
        .set { nanopore_ch }
    
    // Count k-mers for illumina reads
    COUNT_ILLUMINA_KMERS(illumina_ch, CHECK_YAK.out.ready)
    
    // Count k-mers for nanopore reads  
    COUNT_NANOPORE_KMERS(nanopore_ch, CHECK_YAK.out.ready)
    
    // Create all pairwise comparisons between illumina and nanopore samples
    illumina_nanopore_pairs = COUNT_ILLUMINA_KMERS.out.illumina_kmers
        .combine(COUNT_NANOPORE_KMERS.out.nanopore_kmers)
        .map { illumina_sample, illumina_yak, nanopore_sample, nanopore_yak ->
            tuple(illumina_sample, illumina_yak, nanopore_sample, nanopore_yak)
        }
    
    // Compare k-mer spectra
    COMPARE_KMER_SPECTRA(illumina_nanopore_pairs)
    
    // Compute read quality for nanopore reads using illumina as reference
    illumina_nanopore_quality = COUNT_ILLUMINA_KMERS.out.illumina_kmers
        .combine(nanopore_ch)
        .map { illumina_sample, illumina_yak, nanopore_sample, nanopore_reads ->
            tuple(illumina_sample, illumina_yak, nanopore_sample, nanopore_reads)
        }
    
    COMPUTE_READ_QUALITY(illumina_nanopore_quality)
    
    // Generate summary report
    GENERATE_SUMMARY_REPORT(
        COUNT_ILLUMINA_KMERS.out.illumina_hist.collect(),
        COUNT_NANOPORE_KMERS.out.nanopore_hist.collect(),
        COMPARE_KMER_SPECTRA.out.comparison.collect(),
        COMPUTE_READ_QUALITY.out.quality_report.collect()
    )
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    if (workflow.success) {
        log.info"""
        ===================================
        Pipeline completed successfully!
        ===================================
        
        Results are available in: ${params.outdir}
        
        Key outputs:
        - K-mer counts: ${params.outdir}/kmer_counts/
        - Comparisons: ${params.outdir}/comparisons/
        - Quality assessments: ${params.outdir}/quality_assessment/
        - Summary report: ${params.outdir}/quality_assessment_summary.html
        
        """.stripIndent()
    } else {
        log.info"""
        ===================================
        Pipeline failed!
        ===================================
        
        Please check the error messages above and ensure:
        1. YAK is properly installed
        2. Input files exist and are accessible
        3. Sufficient compute resources are available
        
        """.stripIndent()
    }
}