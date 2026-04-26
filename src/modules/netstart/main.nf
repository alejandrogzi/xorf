process NETSTART {
    tag "$meta.id:$meta.name"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/alejandrogzi/orf-net:latest'

    input:
    tuple val(meta), path(sequence)
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.id}*csv"), optional: true, emit: netstart
    tuple val(meta), env(PREDICTION_COUNT),  optional: true, emit: count
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export HF_HOME="\${PWD}/.hf_cache"
    export TRANSFORMERS_CACHE="\${PWD}/.hf_cache"
    export HF_HUB_CACHE="\${PWD}/.hf_cache"
    mkdir -p "\$HF_HOME"

    netstart2 \\
    -in $sequence \\
    -compute_device cpu \\
    -o chordata \\
    -out ${meta.id}_netstart
    $args

    if [ -f ${meta.id}_netstart.csv ]; then
        PREDICTION_COUNT=\$(wc -l < ${meta.id}_netstart.csv)
    else
        PREDICTION_COUNT=0
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netstart2: \$(netstart2 --version 2>&1 | sed 's/.*Version: //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orf-samba: \$(orf --version 2>&1 | sed 's/^.*orf //; s/ .*\$//')
        rnasamba: \$(rnasamba --version 2>&1 | tail -n 1 | sed 's/^rnasamba //')
    END_VERSIONS
    """
}
