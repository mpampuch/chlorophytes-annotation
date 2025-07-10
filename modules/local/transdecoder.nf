process TRANSDECODER_LONGORFS {
    tag "$transcripts"
    label 'process_medium'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5321hdfd78af_4' :
        'biocontainers/transdecoder:5.5.0--pl5321hdfd78af_4' }"

    input:
    path transcripts

    output:
    path "*.transdecoder_dir/longest_orfs.pep"  , emit: longest_orfs
    path "*.transdecoder_dir/longest_orfs.cds"  , emit: longest_cds
    path "*.transdecoder_dir/longest_orfs.gff3" , emit: longest_gff
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    TransDecoder.LongOrfs \\
        -t $transcripts \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(TransDecoder.LongOrfs --version 2>&1 | head -n1 | sed 's/TransDecoder.LongOrfs //')
    END_VERSIONS
    """
}

process TRANSDECODER_PREDICT {
    tag "$transcripts"
    label 'process_medium'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5321hdfd78af_4' :
        'biocontainers/transdecoder:5.5.0--pl5321hdfd78af_4' }"

    input:
    path transcripts
    path blast_results

    output:
    path "*.transdecoder.pep"   , emit: proteins
    path "*.transdecoder.cds"   , emit: cds
    path "*.transdecoder.gff3"  , emit: gff
    path "*.transdecoder.bed"   , emit: bed
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def blast_arg = blast_results ? "--retain_blastp_hits $blast_results" : ""
    """
    # First run LongOrfs if not already done
    TransDecoder.LongOrfs -t $transcripts

    # Predict coding regions
    TransDecoder.Predict \\
        -t $transcripts \\
        $blast_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(TransDecoder.Predict --version 2>&1 | head -n1 | sed 's/TransDecoder.Predict //')
    END_VERSIONS
    """
}