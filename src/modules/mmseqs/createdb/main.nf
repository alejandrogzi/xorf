/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the GNU General Public License, Version 3.0.
*/

process MMSEQS_CREATEDB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0'
        : 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("mmseqsdb"), emit: db
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir mmseqsdb
    mmseqs createdb \\
        ${fasta} \\
        mmseqsdb/protein \\
        --threads ${task.cpus} \\
        ${args}

    mmseqs createindex \\
        mmseqsdb/protein \\
        tmp \\
        --search-type 1 \\
        --threads ${task.cpus}

    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version 2>&1 | sed 's/^.*MMseqs Version: //; s/ .*\$//' | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    mkdir mmseqsdb
    touch mmseqsdb/protein

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs version 2>&1 | sed 's/^.*MMseqs Version: //; s/ .*\$//' | head -n 1)
    END_VERSIONS
    """
}
