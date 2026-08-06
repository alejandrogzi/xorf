// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EMAIL — Sends pipeline completion notifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process EMAIL {
    tag "email"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    val email
    val email_on_fail
    val plaintext_email
    val outdir
    val use_mailx
    tuple path(bed), path(tsv)
    path counts
    path yml

    when:
    task.ext.when == null || task.ext.when

    script:
    // INFO: if use_mailx, all smpt options are ignored
    if (use_mailx) {
        """
        email.py \\
            --bed ${bed} \\
            --tsv ${tsv} \\
            --counts ${counts} \\
            --yml ${yml} \\
            --email ${email} \\
            --email-on-fail ${email_on_fail} \\
            --outdir ${outdir} \\
            --use-mailx
        """
    } else {
        """
        email.py \\
            --bed ${bed} \\
            --tsv ${tsv} \\
            --counts ${counts} \\
            --yml ${yml} \\
            --email ${email} \\
            --email-on-fail ${email_on_fail} \\
            --outdir ${outdir} \\
            --smtp-server ${params.smtp_server} \\
            --smtp-port ${params.smtp_port} \\
            --smtp-user ${params.smtp_user} \\
            --smtp-password ${params.smtp_password} \\
            --smtp-security ${params.smtp_security}
        """
    }
}
