process MODKIT {
    tag "$bam"
    label 'process_medium'

    conda "bioconda::modkit=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/modkit:0.3.0--h9ee0642_0' :
        'biocontainers/modkit:0.3.0--h9ee0642_0' }"

    input:
    path bam
    path genes
    path repeats

    output:
    path "*.modified.bam"       , emit: modified_bam
    path "*.methylation.bed"    , emit: methylation
    path "*.summary.tsv"        , emit: summary
    path "*.entropy.bed"        , emit: entropy
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = bam.baseName
    """
    # Calculate methylation probability distribution
    modkit sample-probs \\
        --no-sampling \\
        --only-mapped \\
        $bam > ${prefix}.prob_dist.txt

    # Calculate percentile threshold
    THRESHOLD=\$(awk -v p=${params.modkit_percentile} 'NR==int(NR*p/100)+1{print \$1; exit}' ${prefix}.prob_dist.txt)

    # Call modifications with threshold
    modkit call-mods \\
        --filter-threshold \$THRESHOLD \\
        --threads $task.cpus \\
        $bam \\
        ${prefix}.modified.bam

    # Generate methylation summary
    modkit summary \\
        --no-sampling \\
        --only-mapped \\
        --tsv \\
        --include-bed $genes \\
        --filter-threshold \$THRESHOLD \\
        ${prefix}.modified.bam > ${prefix}.summary.tsv

    # Generate methylation calls
    modkit pileup \\
        --threads $task.cpus \\
        --filter-threshold \$THRESHOLD \\
        ${prefix}.modified.bam \\
        ${prefix}.methylation.bed

    # Calculate methylation entropy
    modkit entropy \\
        --regions $genes \\
        --filter-threshold \$THRESHOLD \\
        --threads $task.cpus \\
        ${prefix}.modified.bam \\
        ${prefix}.entropy.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version | head -n1 | sed 's/modkit //')
    END_VERSIONS
    """
}