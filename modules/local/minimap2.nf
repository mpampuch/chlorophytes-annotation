process MINIMAP2 {
    tag "$nanopore_reads"
    label 'process_high'

    conda "bioconda::minimap2=2.28 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0' :
        'biocontainers/minimap2:2.28--he4a0461_0' }"

    input:
    path nanopore_reads
    path genome

    output:
    path "*.bam"        , emit: bam
    path "*.bai"        , emit: bai
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = nanopore_reads.baseName
    """
    # Index genome if needed
    if [ ! -f "${genome}.mmi" ]; then
        minimap2 -d ${genome}.mmi $genome
    fi

    # Map nanopore reads
    minimap2 \\
        $args \\
        -t $task.cpus \\
        ${genome}.mmi \\
        $nanopore_reads | \\
    samtools view -bS - | \\
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -

    # Index BAM file
    samtools index ${prefix}.sorted.bam

    mv ${prefix}.sorted.bam ${prefix}.bam
    mv ${prefix}.sorted.bam.bai ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
}