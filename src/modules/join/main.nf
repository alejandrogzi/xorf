// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    JOIN — Merges NetStart and TransAID predictions into a unified net file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process JOIN {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/alejandrogzi/orf-chunk:latest'

    input:
    tuple val(meta), path(transaid), path(bed)
    tuple val(meta1), path(netstart)

    output:
    tuple val(meta), path('net/merged.net'), emit: net
    tuple val(meta), env(PREDICTION_COUNT), emit: count
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def nt = ${netstart} ? '--netstart $netstart' : ''
    """
    orf net \\
    $args \\
    --bed $bed \\
    --transaid $transaid \\
    $nt \\
    --outdir .

    PREDICTION_COUNT=\$(cat net/merged.net | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-net: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir net
    touch net/merged.net

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-chunk: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
    END_VERSIONS
    """
}
