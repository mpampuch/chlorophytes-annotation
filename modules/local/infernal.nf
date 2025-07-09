process INFERNAL {
    tag "$genome"
    label 'process_high'

    conda "bioconda::infernal=1.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/infernal:1.1.4--h779adbc_0' :
        'biocontainers/infernal:1.1.4--h779adbc_0' }"

    input:
    path genome
    path rfam_db

    output:
    path "*.tblout"     , emit: families
    path "*.txt"        , emit: filtered_output
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = genome.baseName
    """
    # Run cmscan against RFAM database
    cmscan \\
        --cpu $task.cpus \\
        --tblout ${prefix}.tblout \\
        $args \\
        $rfam_db \\
        $genome

    # Filter results by E-value threshold
    awk '\$3 <= ${params.evalue_threshold}' ${prefix}.tblout > ${prefix}.filtered.tblout

    # Remove partial/fragmented matches
    awk '\$8 == \$9' ${prefix}.filtered.tblout > ${prefix}.complete_matches.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: \$(cmscan -h | grep "# INFERNAL" | sed 's/# INFERNAL //')
    END_VERSIONS
    """
}