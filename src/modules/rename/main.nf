/*
Copyright (c) 2026 The Hiller Lab at the Senckenberg Gessellschaft für Naturforschung
Distributed under the terms of the GNU General Public License, Version 3.0.
*/

process RENAME_PREDICTIONS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta1), path(tsv)

    output:
    tuple val(meta), path("*.renamed.bed"), path("*.renamed.tsv"), emit: files
    tuple val(meta), path("*.renamed.bed")  , emit: bed
    tuple val(meta), path("*.renamed.tsv")  , emit: tsv
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    def unmask_tag = meta.masked ? "" : "--unmask"
    """
    rename_predictions.py \
      --bed $bed \
      --tsv $tsv \
      --output-bed ${prefix}.renamed.bed \
      --output-tsv ${prefix}.renamed.tsv \
      $unmask_tag \
      $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rename_predictions: \$(rename_predictions.py --version | sed 's/rename_predictions.py //g')
    END_VERSIONS
    """

    stub:
    """
    touch dummy.renamed.bed
    touch dummy.renamed.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rename_predictions: \$(rename_predictions.py --version | sed 's/rename_predictions.py //g')
    END_VERSIONS
    """
}
