process VCFTOOLS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vcftools=0.1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcftools:0.1.17--pl5321hdfd78af_4' :
        'quay.io/biocontainers/vcftools:0.1.17--pl5321hdfd78af_4' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filtered.vcf"), emit: filtered_vcf
    path "*.log"                           , emit: log
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcftools \\
        --vcf $vcf \\
        $args \\
        --recode \\
        --recode-INFO-all \\
        --out ${prefix}.filtered

    mv ${prefix}.filtered.recode.vcf ${prefix}.filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(vcftools --version 2>&1 | head -n1 | sed 's/^.*VCFtools (//' | sed 's/).*//')
    END_VERSIONS
    """
}