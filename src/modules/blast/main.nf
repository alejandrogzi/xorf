// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BLAST — Validates ORFs against a reference database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BLAST {
    tag "$meta.id:$meta.name"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/alejandrogzi/orf-blast:latest'

    input:
    tuple val(meta), path(bed), path(sequence), path(predictions), path(net)
    each path(database)

    output:
    tuple val(meta), path(bed), path("${meta.id}/*result"),             optional: true, emit: blast
    tuple val(meta), env(INITIAL_REGION_COUNT), env(TRANSLATION_COUNT), optional: true, emit: counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def orf_min_len = task.ext.orf_min_len ?: 50
    def orf_min_percent = task.ext.orf_min_percent ?: 0.25
    def upstream = task.ext.upstream ?: 1000
    def downstream = task.ext.downstream ?: 1000
    """
    set +e
    orf blast \\
    --fasta $sequence \\
    --bed $bed \\
    --tai $predictions \\
    --net $net \\
    --outdir ${meta.id} \\
    --orf-min-len $orf_min_len \\
    --orf-min-percent $orf_min_percent \\
    --database $database \\
    --upstream-flank $upstream \\
    --downstream-flank $downstream \\
    $args > blast.stdout 2> blast.stderr
    blast_status=\$?
    set -e

    cat blast.stdout
    cat blast.stderr >&2

    shopt -s nullglob
    result_files=( ${meta.id}/orf/*result )
    blast_empty=false

    if (( \${#result_files[@]} > 0 )); then
      mv "\${result_files[@]}" ${meta.id}/
      rm -rf ${meta.id}/orf
    elif [[ -f ${meta.id}/orf/orf.dedup.pep && ! -s ${meta.id}/orf/orf.dedup.pep ]]; then
      if [[ \$blast_status -eq 0 ]] || grep -q 'Input file seems to be empty' blast.stderr; then
        blast_empty=true
        echo "INFO: BLAST produced no candidates after deduplication for ${meta.id}:${meta.name}; skipping output emission." >&2
      fi
    fi

    if [[ \$blast_status -ne 0 && \$blast_empty != true ]]; then
      exit \$blast_status
    fi

    if [[ \$blast_status -eq 0 && \${#result_files[@]} -eq 0 && \$blast_empty != true ]]; then
      echo "ERROR: BLAST finished without result files for ${meta.id}:${meta.name}." >&2
      exit 1
    fi

    if [[ \${#result_files[@]} -gt 0 ]]; then
      INITIAL_REGION_COUNT=\$(wc -l < ${bed})
      TRANSLATION_COUNT=\$(wc -l < ${predictions})
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-blast: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        diamond: \$(diamond version 2>&1 | sed 's/^.*diamond version //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}
    touch ${meta.id}/orf
    touch ${meta.id}/orf/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-blast: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        diamond: \$(diamond version 2>&1 | sed 's/^.*diamond version //' )
    END_VERSIONS
    """
}
