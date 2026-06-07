// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RNASAMBA — Classifies ORFs as coding or non-coding using machine learning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RNASAMBA {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/alejandrogzi/orf-samba:latest'

    input:
    tuple val(meta), path(bed), path(sequence)
    tuple val(_), path(weights)

    output:
    tuple val(meta), path("${meta.id}_${meta.name}/*tsv")      , optional: true, emit: samba
    tuple val(meta), path("${meta.id}_${meta.name}/*strip.fa") , optional: true, emit: fasta
    tuple val(meta), path(bed)                    , optional: true, emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.name}"
    """
    orf samba \\
    --fasta $sequence \\
    --outdir ${prefix} \\
    --upstream-flank $upstream \\
    --downstream-flank $downstream \\
    --weights $weights \\
    $args

    mv ${prefix}/samba/*tsv ${prefix}/${prefix}.${meta.name}.samba.tsv && rm -rf ${prefix}/samba
    mv ${meta.name}.tmp.strip.fa ${prefix}/${prefix}.${meta.name}.strip.fa 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}
    touch ${prefix}/*strip.fa
    touch ${prefix}/samba
    touch ${prefix}/samba/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """
}
