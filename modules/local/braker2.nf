process BRAKER2 {
    tag "$genome"
    label 'process_high_memory'

    conda "bioconda::braker2=2.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/braker2:2.1.5--py38pl5321hdfd78af_3' :
        'biocontainers/braker2:2.1.5--py38pl5321hdfd78af_3' }"

    input:
    path genome
    path bam_file
    path proteins

    output:
    path "braker.gff3"      , emit: gff
    path "braker.gtf"       , emit: gtf
    path "braker.aa"        , emit: proteins
    path "braker.codingseq" , emit: cds
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create output directory
    mkdir -p braker_output

    # Run BRAKER2
    braker.pl \\
        --genome=$genome \\
        --bam=$bam_file \\
        --prot_seq=$proteins \\
        --workingdir=braker_output \\
        --threads=$task.cpus \\
        $args

    # Move output files
    cp braker_output/braker.gff3 .
    cp braker_output/braker.gtf .
    cp braker_output/braker.aa .
    cp braker_output/braker.codingseq .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        braker2: \$(braker.pl --version 2>&1 | grep -o 'BRAKER version [0-9.]*' | sed 's/BRAKER version //')
    END_VERSIONS
    """
}