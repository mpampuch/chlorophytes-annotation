process STAR_GENOMEGENERATE {
    tag "$genome"
    label 'process_high'

    conda "bioconda::star=2.7.9a"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0' :
        'biocontainers/star:2.7.9a--h9ee0642_0' }"

    input:
    path genome

    output:
    path "star"         , emit: index
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir star

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $genome \\
        --runThreadN $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star=2.7.9a bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0' :
        'biocontainers/star:2.7.9a--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.bai")        , emit: bai
    tuple val(meta), path("*Log.final.out"), emit: log_final
    tuple val(meta), path("*Log.out")     , emit: log_out
    tuple val(meta), path("*SJ.out.tab")  , emit: junction
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix ${prefix}. \\
        $args

    # Convert SAM to BAM and sort
    samtools view -bS ${prefix}.Aligned.out.sam | \\
    samtools sort -@ $task.cpus -o ${prefix}.sorted.bam -

    # Index BAM
    samtools index ${prefix}.sorted.bam

    mv ${prefix}.sorted.bam ${prefix}.bam
    mv ${prefix}.sorted.bam.bai ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
}