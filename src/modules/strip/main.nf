process STRIP_OCCURRENCES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-22105f082fdf56207f1dcc5b6da71a14394e28d7:387d955c0a2cdb831ec519d636e4ffd7062d6ae1-0':
        'biocontainers/mulled-v2-22105f082fdf56207f1dcc5b6da71a14394e28d7:387d955c0a2cdb831ec519d636e4ffd7062d6ae1-0' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta1), path(descriptor)
    val hint

    output:
    tuple val(meta), path("*.hq.bed"),   optional: true, emit: hq
    tuple val(meta2), path("*.discard.bed"),    optional: true, emit: discard
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    meta2 = meta.clone()
    meta2.id  = "${meta.id}.discard"
    """
    grep '${hint}' ${descriptor} | cut -f1 | \
    awk 'NR==FNR {ids[\$1]; next} \$4 in ids' - ${bed} \
    > ${prefix}.discard.bed   || [[ \$? == 1 ]]

    grep -v '${hint}' ${descriptor} | cut -f1 | \
    awk 'NR==FNR {ids[\$1]; next} \$4 in ids' - ${bed} \
    > ${prefix}.hq.bed || [[ \$? == 1 ]]

    if [[ ! -s ${prefix}.hq.bed ]]; then
        rm ${prefix}.hq.bed
    fi

    if [[ ! -s ${prefix}.discard.bed ]]; then
        rm ${prefix}.discard.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | sed 's/grep (GNU grep) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hq.bed
    touch ${prefix}.discard.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | sed 's/grep (GNU grep) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
