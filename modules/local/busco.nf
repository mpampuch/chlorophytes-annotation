process BUSCO {
    tag "$proteins"
    label 'process_medium'

    conda "bioconda::busco=5.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.4.4--pyhdfd78af_0' :
        'biocontainers/busco:5.4.4--pyhdfd78af_0' }"

    input:
    path proteins

    output:
    path "busco_*"              , emit: results
    path "short_summary*.txt"   , emit: summary
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = proteins.baseName
    """
    # Run BUSCO for Eukaryota
    busco \\
        -i $proteins \\
        -l ${params.busco_lineage_eukaryota} \\
        -o busco_eukaryota \\
        -m protein \\
        --cpu $task.cpus \\
        $args

    # Run BUSCO for Stramenopiles  
    busco \\
        -i $proteins \\
        -l ${params.busco_lineage_stramenopiles} \\
        -o busco_stramenopiles \\
        -m protein \\
        --cpu $task.cpus \\
        $args

    # Copy summary files
    cp busco_eukaryota/short_summary.*.txt short_summary_eukaryota.txt
    cp busco_stramenopiles/short_summary.*.txt short_summary_stramenopiles.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>&1 | sed 's/BUSCO //')
    END_VERSIONS
    """
}