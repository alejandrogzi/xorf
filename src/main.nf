#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE — Orchestrates XORF analysis and sends email notifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EMAIL }        from './modules/email/main.nf'
include { XORF }         from './subworkflows/xorf/main.nf'
include { WGET as WGET_SAMBA_WEIGHTS } from './modules/wget/main.nf'

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

    XORF (
       Channel.fromPath(params.regions).map { it -> [ [id: it.baseName, chr:'xorf'], it ] },
       Channel.fromPath(params.sequence),
       Channel.fromPath(params.database),
       params.outdir,
       params.chunk_size,
       ch_samba_weights,
       params.predict_keep_raw,
       params.selenocysteine_sites
    )

    PIPELINE_COMPLETION (
        params.email_to,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.use_mailx,
        XORF.out.files,
        XORF.out.counts,
        XORF.out.versions
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
