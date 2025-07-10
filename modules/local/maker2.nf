process MAKER2 {
    tag "$genome"
    label 'process_high_memory'

    conda "bioconda::maker=3.01.03"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/maker:3.01.03--pl5321h9ee0642_0' :
        'biocontainers/maker:3.01.03--pl5321h9ee0642_0' }"

    input:
    path genome
    path braker_gff
    path protein_alignments

    output:
    path "*.gff"        , emit: gff
    path "*.fasta"      , emit: proteins
    path "*genes.bed"   , emit: genes
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = genome.baseName
    """
    # Create MAKER control files
    maker -CTL

    # Configure maker_opts.ctl
    sed -i "s|genome=|genome=$genome|" maker_opts.ctl
    sed -i "s|pred_gff=|pred_gff=$braker_gff|" maker_opts.ctl
    sed -i "s|protein=|protein=$protein_alignments|" maker_opts.ctl
    sed -i "s|cpus=1|cpus=$task.cpus|" maker_opts.ctl

    # Run MAKER
    maker $args

    # Process outputs
    gff3_merge -d ${prefix}.maker.output/${prefix}_master_datastore_index.log
    fasta_merge -d ${prefix}.maker.output/${prefix}_master_datastore_index.log

    # Generate gene BED file
    maker2bed < ${prefix}.all.gff > ${prefix}.genes.bed

    # Rename outputs
    mv ${prefix}.all.gff ${prefix}.maker.gff
    mv ${prefix}.all.maker.proteins.fasta ${prefix}.proteins.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maker: \$(maker --version 2>&1 | head -n1 | sed 's/MAKER version //')
    END_VERSIONS
    """
}