process TRINITY {
    tag "$bam"
    label 'process_high_memory'

    conda "bioconda::trinity=2.11.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trinity:2.11.0--h00214ad_1' :
        'biocontainers/trinity:2.11.0--h00214ad_1' }"

    input:
    path bam
    path genome

    output:
    path "Trinity.fasta"        , emit: transcripts
    path "Trinity.fasta.gene_trans_map" , emit: gene_map
    path "trinity_stats.txt"    , emit: stats
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Convert BAM to FASTQ
    samtools sort -n $bam -o sorted.bam
    bedtools bamtofastq -i sorted.bam -fq reads_1.fq -fq2 reads_2.fq

    # Run Trinity
    Trinity \\
        --genome_guided_bam $bam \\
        --genome_guided_max_intron 100000 \\
        --left reads_1.fq \\
        --right reads_2.fq \\
        --CPU $task.cpus \\
        --max_memory ${task.memory.toGiga()}G \\
        $args

    # Generate statistics
    \$TRINITY_HOME/util/TrinityStats.pl Trinity.fasta > trinity_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trinity: \$(Trinity --version 2>&1 | grep -o 'Trinity version: Trinity-v[0-9.]*' | sed 's/Trinity version: Trinity-v//')
    END_VERSIONS
    """
}