// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRANSAID — Predicts translation initiation sites using deep learning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRANSAID {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/alejandrogzi/orf-net:latest'

    input:
    tuple val(meta), path(sequence)
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*csv")            , optional: true, emit: transaid
    tuple val(meta), env(PREDICTION_COUNT)   , optional: true, emit: count
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    transaid \\
    --input $sequence \\
    --gpu -1 \\
    --output ${meta.id}_${meta.name}_transaid \\
    $args

    mv *csv ${meta.id}_${meta.name}_transaid.csv
    PREDICTION_COUNT=\$(wc -l < ${meta.id}_${meta.name}_transaid.csv)
    
    rm *.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transaid: \$(transaid --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         transaid: \$(transaid --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """
}
