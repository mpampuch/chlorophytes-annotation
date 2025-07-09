process DIAMOND_BLAST {
    tag "$query"
    label 'process_high'

    conda "bioconda::diamond=0.9.29"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:0.9.29--ha0959df_0' :
        'biocontainers/diamond:0.9.29--ha0959df_0' }"

    input:
    path query
    path database

    output:
    path "*.diamond.txt"    , emit: blast_results
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = query.baseName
    """
    # Create diamond database if needed
    if [ ! -f "${database}.dmnd" ]; then
        diamond makedb --in $database -d ${database}
    fi

    # Run DIAMOND BLAST
    diamond blastp \\
        --query $query \\
        --db ${database}.dmnd \\
        --out ${prefix}.diamond.txt \\
        --outfmt 6 \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version | head -n1 | sed 's/diamond version //')
    END_VERSIONS
    """
}