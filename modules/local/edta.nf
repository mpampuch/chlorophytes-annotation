process EDTA {
    tag "$genome"
    label 'process_high_memory'

    conda "bioconda::edta=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/edta:2.0.0--hdfd78af_0' :
        'biocontainers/edta:2.0.0--hdfd78af_0' }"

    input:
    path genome
    path cds_sequences

    output:
    path "*.EDTA.TElib.fa"     , emit: te_library
    path "*.EDTA.TEanno.gff3"  , emit: te_annotation
    path "*.EDTA.raw"          , emit: raw_output
    path "*.EDTA.final"        , emit: final_output
    path "*.EDTA.mod"          , emit: masked_genome
    path "*.EDTA.out"          , emit: repeats
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = genome.baseName
    """
    EDTA.pl \\
        --genome $genome \\
        --cds $cds_sequences \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EDTA: \$(EDTA.pl --version 2>&1 | head -n1 | sed 's/EDTA version: //')
    END_VERSIONS
    """
}