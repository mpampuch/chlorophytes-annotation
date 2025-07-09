process FREEBAYES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::freebayes=1.3.4 bioconda::bwa=0.7.17 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ad317f19f5881324e963f6a2d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0' :
        'biocontainers/mulled-v2-ad317f19f5881324e963f6a2d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0' }"

    input:
    tuple val(meta), path(reads)
    path genome

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Index genome if not already indexed
    if [ ! -f "${genome}.fai" ]; then
        samtools faidx $genome
    fi

    # Align reads to genome
    bwa mem -t $task.cpus $genome ${reads[0]} ${reads[1]} | \\
    samtools view -bS - | \\
    samtools sort -@ $task.cpus -o ${prefix}.bam -

    # Index BAM file
    samtools index ${prefix}.bam

    # Run FreeBayes
    freebayes \\
        $args \\
        -f $genome \\
        ${prefix}.bam > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(freebayes --version | head -n 1 | sed 's/version:  //')
        bwa: \$(bwa 2>&1 | grep -E '^Version' | sed 's/Version: //')
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
}