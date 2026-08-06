// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the GNU General Public License, Version 3.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONCAT — Concatenates BED and TSV files from parallel processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CONCAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'quay.io/biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(beds, stageAs: 'bed/*'), path(tsvs, stageAs: 'tsv/*')

    output:
    tuple val(meta), path("*bed"), path("*tsv"), optional: true, emit: files
    tuple val(meta), path("*.bed"), optional:true, emit: bed
    tuple val(meta), path("*.tsv"), optional:true, emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: [meta.id, meta.name].findAll { it }.join('.')
    """
    if [ -d bed ]; then
      cat bed/*.bed > ${prefix}.bed
      awk 'NR==1 || FNR>1' tsv/*.tsv > ${prefix}.tsv

      rm -rf bed
      rm -rf tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.bed
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/cat (GNU coreutils) //' )
    END_VERSIONS
    """
}
