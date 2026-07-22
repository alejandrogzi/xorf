#!/usr/bin/env nextflow

/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the Apache License, Version 2.0.
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xorf

    end-to-end robust and comprehensive ORF prediction pipeline
    Authors: Alejandro Gonzales-Irribarren, Michael Hiller

    GitHub:  https://github.com/hillerlab/xorf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    log.info """
    xorf v${workflow.manifest.version}
    End-to-end robust and comprehensive ORF prediction pipeline

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

    Usage (full run):
        nextflow run main.nf \\
            --regions        /path/to/regions.{bed, gtf, gff} \\
            --sequence       /path/to/genome.{fa, 2bit, fa.gz} \\
            --outdir         results/ \\
            -profile         apptainer,slurm

    Pass all parameters from a JSON file (replaces old --params_from_file):
        nextflow run main.nf -params-file my_params.json

    Required parameters (full run + fill/clean):
        --regions            PATH    Path to genomic regions (BED, GTF, or GFF)
        --sequence           PATH    Path to genome sequence (FASTA or 2bit)

    Optional parameters (common):
        --outdir              PATH    Output directory [default: ./results]
        --custom_database     PATH    Path to custom protein database (.dmnd) [default: null]
        --chunk_size          INT     Chunk size for parallel processing [default: 20]
        --predict_keep_raw    BOOL      Keep raw predictions [default: false]
        --selenocysteine_sites PATH      Selenocysteine masking [default: null]
        --predict_min_score_max_predictions FLOAT   Minimum score for ORF predictions [default: 0.50]
        --predict_max_predictions INT  Maximum number of ORF predictions [default: 3]
        --skip_netstart       BOOL      Skip netstart [default: false]

    Profiles:
        local       Run on local machine (default)
        slurm       Submit jobs to SLURM cluster
        conda       Use conda environments
        apptainer   Use Apptainer containers
        singularity Use Singularity containers
        docker      Use Docker containers
        test        Run with bundled test data

    Use --help to show this message.
    """.stripIndent()
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EMAIL }        from './modules/email/main.nf'
include { XORF as MAIN }         from './subworkflows/xorf/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from './modules/wget/main.nf'
include { UNTAR } from './modules/untar/main.nf'
include { WGET as WGET_PROTEIN_DATABASE } from './modules/wget/main.nf'
include { GUNZIP as GUNZIP_DATABASE } from './modules/gunzip/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateRun() {
    def errors = []
    if (!params.regions)   errors << "  --regions is required"
    if (!params.sequence)  errors << "  --sequence is required"

    if (errors) {
        log.error "Parameter validation failed:\n${errors.join('\n')}"
        System.exit(1)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
      email
      email_on_fail
      plaintext_email
      outdir
      use_mailx
      files
      counts
      ch_versions

    main:

      if (params.sent_email) {
          EMAIL (
              email,
              email_on_fail,
              plaintext_email,
              outdir,
              use_mailx,
              files,
              counts,
              ch_versions
          )
      }

      workflow.onError {
          log.error "ERROR: Pipeline failed!"
          log.error "ERROR: Please check the following error message:\n${workflow.errorMessage}"
          log.error "ERROR: Refer to github issues: https://github.com/alejandrogzi/HLorf/issues"
      }

      workflow.onComplete {
          log.info "\nPipeline completed successfully!"
      }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    XORF()
}

workflow XORF {
    
    log.info """
    > xorf v${workflow.manifest.version}
    > End-to-end robust and comprehensive ORF prediction pipeline

    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

      Regions:   ${params.regions}
      Sequence:  ${params.sequence}
      Database:  ${params.database} (custom: ${params.custom_database})
      Outdir:    ${params.outdir}
      Profile:   ${workflow.profile}
    """.stripIndent()

    validateRun()

    ch_samba_weights = Channel.empty()
    if (params.samba_local_weights) {
        ch_samba_weights = Channel.value(
          file(params.samba_local_weights, checkIfExists: true)
        ).map { path -> [ [id : path.baseName ], path ] }
    } else {
        WGET_SAMBA_WEIGHTS(
          Channel.value(
            params.samba_weights
          ).map { url -> [ [id : url.tokenize('/')[-1]], url ] }
        )
        ch_samba_weights = WGET_SAMBA_WEIGHTS.out.outfile
    }

    ch_database = Channel.empty()
    if (params.custom_database) {
      if (params.custom_database.endsWith('.gz')) {
          GUNZIP_DATABASE(
              Channel.value(
                  [ [id: params.custom_database.tokenize('/')[-1]], params.custom_database ]
              )
          )
          GUNZIP_DATABASE.out.gunzip
            .map { meta, it -> it }
            .set { ch_database }
      } else {
          ch_database = Channel.fromPath(params.custom_database)
      }
    } else {
      WGET_PROTEIN_DATABASE(
        Channel.value(params.database)
        .map { it -> [ [ id: 'uniprot_sprot.tar.gz' ], it ] }
      )
      UNTAR(WGET_PROTEIN_DATABASE.out.outfile)
      ch_database = UNTAR.out.contents.map { meta, it -> it }
    }

    MAIN (
       Channel.fromPath(params.regions).map { it -> [ [id: it.baseName, chr: randomHash()], it ] },
       Channel.fromPath(params.sequence),
       ch_database,
       params.outdir,
       params.chunk_size,
       ch_samba_weights,
       params.predict_keep_raw,
       params.selenocysteine_sites,
       params.skip_netstart,
       params.rename_deactivate
    )

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.use_mailx,
        MAIN.out.files,
        MAIN.out.counts,
        MAIN.out.versions
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def randomHash() {
    def chars = ('0'..'9') + ('a'..'z') + ('A'..'Z')
    def random = new Random()
    return (1..6).collect { chars[random.nextInt(chars.size())] }.join()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (workflow.success) {
        def results_dir = new File(params.outdir as String, '04_results')
        def final_beds = results_dir.exists() ? (results_dir.listFiles()?.findAll { it.name.endsWith('.bed') } ?: []) : []
        log.info "Pipeline completed successfully!"
        if (final_beds) {
            log.info "Final predictions: ${final_beds.collect { it.toString() }.join(', ')}"
        } else {
            log.warn "Pipeline reported success but final bed file was not produced - check that all steps ran"
        }
        log.info "Run time   : ${workflow.duration}"
    } else {
        log.error "Pipeline FAILED — ${workflow.errorMessage}"
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
