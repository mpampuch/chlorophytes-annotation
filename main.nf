#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FREEBAYES }                    from './modules/local/freebayes'
include { VCFTOOLS }                     from './modules/local/vcftools'
include { INFERNAL }                     from './modules/local/infernal'
include { EDTA }                         from './modules/local/edta'
include { ASGART }                       from './modules/local/asgart'
include { BEDTOOLS_INTERSECT }           from './modules/local/bedtools'
include { STAR_GENOMEGENERATE }          from './modules/local/star'
include { STAR_ALIGN }                   from './modules/local/star'
include { TRINITY }                      from './modules/local/trinity'
include { CDHIT_EST }                    from './modules/local/cdhit'
include { TRANSDECODER_LONGORFS }        from './modules/local/transdecoder'
include { DIAMOND_BLAST }                from './modules/local/diamond'
include { TRANSDECODER_PREDICT }         from './modules/local/transdecoder'
include { BRAKER2 }                      from './modules/local/braker2'
include { EXONERATE }                    from './modules/local/exonerate'
include { GENEMARK }                     from './modules/local/genemark'
include { MAKER2 }                       from './modules/local/maker2'
include { BUSCO }                        from './modules/local/busco'
include { PANNZER2 }                     from './modules/local/pannzer2'
include { KAAS }                         from './modules/local/kaas'
include { MINIMAP2 }                     from './modules/local/minimap2'
include { MODKIT }                       from './modules/local/modkit'
include { MBTOOLS }                      from './modules/local/mbtools'
include { RAXMLNG }                      from './modules/local/raxmlng'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHLOROPHYTES_ANNOTATION {
    
    take:
    genome_fasta           // Path to reference genome FASTA
    illumina_reads         // Channel of Illumina paired-end reads
    nanopore_reads         // Path to Nanopore reads with 5mCG modifications
    rna_reads              // Channel of RNA-seq reads
    protein_db             // Path to Stramenopiles protein database
    rfam_db                // Path to RFAM database

    main:
    
    // ====================
    // PLOIDY EVALUATION
    // ====================
    
    FREEBAYES(
        illumina_reads,
        genome_fasta
    )
    
    VCFTOOLS(
        FREEBAYES.out.vcf
    )
    
    // ====================
    // RNA FAMILIES PREDICTION
    // ====================
    
    INFERNAL(
        genome_fasta,
        rfam_db
    )
    
    // ====================
    // TRANSPOSABLE ELEMENTS
    // ====================
    
    EDTA(
        genome_fasta,
        params.tpseudnana_cds
    )
    
    // ====================
    // SEGMENTAL DUPLICATION
    // ====================
    
    ASGART(
        genome_fasta
    )
    
    BEDTOOLS_INTERSECT(
        ASGART.out.duplications,
        MAKER2.out.genes
    )
    
    // ====================
    // TRANSCRIPTOME ASSEMBLY
    // ====================
    
    STAR_GENOMEGENERATE(
        genome_fasta
    )
    
    STAR_ALIGN(
        rna_reads,
        STAR_GENOMEGENERATE.out.index
    )
    
    TRINITY(
        STAR_ALIGN.out.bam,
        genome_fasta
    )
    
    CDHIT_EST(
        TRINITY.out.transcripts
    )
    
    TRANSDECODER_LONGORFS(
        CDHIT_EST.out.clustered
    )
    
    DIAMOND_BLAST(
        TRANSDECODER_LONGORFS.out.longest_orfs,
        protein_db
    )
    
    TRANSDECODER_PREDICT(
        CDHIT_EST.out.clustered,
        DIAMOND_BLAST.out.blast_results
    )
    
    // ====================
    // GENOME ANNOTATION
    // ====================
    
    BRAKER2(
        EDTA.out.masked_genome,
        STAR_ALIGN.out.bam,
        TRANSDECODER_PREDICT.out.proteins
    )
    
    EXONERATE(
        TRANSDECODER_PREDICT.out.proteins,
        genome_fasta
    )
    
    GENEMARK(
        EDTA.out.masked_genome
    )
    
    MAKER2(
        EDTA.out.masked_genome,
        BRAKER2.out.gff,
        EXONERATE.out.alignments
    )
    
    BUSCO(
        MAKER2.out.proteins
    )
    
    // ====================
    // FUNCTIONAL ANNOTATION
    // ====================
    
    PANNZER2(
        MAKER2.out.proteins
    )
    
    KAAS(
        MAKER2.out.proteins
    )
    
    // ====================
    // NANOPORE MAPPING & METHYLATION
    // ====================
    
    MINIMAP2(
        nanopore_reads,
        genome_fasta
    )
    
    MODKIT(
        MINIMAP2.out.bam,
        MAKER2.out.genes,
        EDTA.out.repeats
    )
    
    MBTOOLS(
        MODKIT.out.modified_bam,
        MAKER2.out.genes
    )
    
    // ====================
    // PHYLOGENOMICS
    // ====================
    
    RAXMLNG(
        MAKER2.out.proteins,
        params.reference_genomes
    )
    
    emit:
    vcf_filtered           = VCFTOOLS.out.filtered_vcf
    rna_families           = INFERNAL.out.families
    repeats                = EDTA.out.repeats
    masked_genome          = EDTA.out.masked_genome
    segmental_duplications = ASGART.out.duplications
    transcriptome          = TRINITY.out.transcripts
    gene_annotation        = MAKER2.out.gff
    proteins               = MAKER2.out.proteins
    busco_results          = BUSCO.out.results
    functional_annotation  = PANNZER2.out.annotation
    kegg_annotation        = KAAS.out.kegg
    methylation_calls      = MODKIT.out.methylation
    phylogenetic_tree      = RAXMLNG.out.tree
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    // Input channels
    genome_fasta = file(params.genome_fasta)
    
    illumina_reads = Channel
        .fromFilePairs(params.illumina_reads, checkIfExists: true)
        .map { meta, reads -> [meta, reads] }
    
    nanopore_reads = file(params.nanopore_reads)
    
    rna_reads = Channel
        .fromFilePairs(params.rna_reads, checkIfExists: true)
        .map { meta, reads -> [meta, reads] }
    
    protein_db = file(params.protein_db)
    rfam_db = file(params.rfam_db)
    
    // Run main workflow
    CHLOROPHYTES_ANNOTATION(
        genome_fasta,
        illumina_reads,
        nanopore_reads,
        rna_reads,
        protein_db,
        rfam_db
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}