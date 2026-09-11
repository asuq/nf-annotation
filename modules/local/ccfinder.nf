/*
 * Run CRISPRCasFinder with the locked standalone invocation semantics.
 */
process CCFINDER {
    tag "${meta.accession}"
    label 'process_small'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/ccfinder" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['result.json', 'ccfinder.log']
                ? filename
                : null
        },
    )

    input:
    tuple val(meta), path(genome), val(gcode)

    output:
    tuple val(meta), path('ccfinder'), path('result.json'), path('ccfinder.log'), emit: results
    tuple val(meta), path('result.json'), emit: result_json
    path 'versions.yml', emit: versions

    script:
    def ccfinderRoot = '/usr/local/CRISPRCasFinder'
    def extraArgs = (params.ccfinder_extra_args ?: '').toString()
    """
    export TMPDIR="\$PWD/ccfinder_tmp"
    mkdir -p "\${TMPDIR}"

    ccfinder_root='${ccfinderRoot}'
    if [[ ! -f "\${ccfinder_root}/CRISPRCasFinder.pl" ]]; then
        echo "CRISPRCasFinder.pl not found under \${ccfinder_root}." >&2
        exit 1
    fi

    task_root="\$PWD"
    run_root="\${task_root}/ccfinder_run"
    tool_output_root="\${task_root}/ccfinder_raw"
    genome_path="\$(cd "\$(dirname "${genome}")" && pwd)/\$(basename "${genome}")"
    genome_name="\$(basename "${genome}")"

    mkdir -p "\${run_root}"

    awk -v output_root="\${task_root}" '
        /^>/ {
            if (output_path != "") {
                close(output_path)
            }
            contig_id = \$1
            sub(/^>/, "", contig_id)
            sub(/\\.[0-9]+\$/, "", contig_id)
            output_path = output_root "/" contig_id ".fna"
            print \$0 > output_path
            next
        }
        {
            sequence_line = \$0
            while (length(sequence_line) > 60) {
                print substr(sequence_line, 1, 60) >> output_path
                sequence_line = substr(sequence_line, 61)
            }
            if (length(sequence_line) > 0) {
                print sequence_line >> output_path
            }
        }
    ' "\${genome_path}"

    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    wrapper_log="ccfinder_wrapper.log"
    : > "\${wrapper_log}"
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> "\${wrapper_log}"
        rm -rf "\${run_root}" "\${tool_output_root}"
        mkdir -p "\${run_root}"
        pushd "\${run_root}" >/dev/null
        cp "\${genome_path}" "\${genome_name}"
        : > index.html
        set +e
        perl "\${ccfinder_root}/CRISPRCasFinder.pl" -in "\${genome_name}" \
            -outdir "\${tool_output_root}" \
            -soFile "\${ccfinder_root}/sel392v2.so" \
            -DBcrispr "\${ccfinder_root}/supplementary_files/CRISPR_crisprdb.csv" \
            -repeats "\${ccfinder_root}/supplementary_files/Repeat_List.csv" \
            -DIRrepeat "\${ccfinder_root}/supplementary_files/repeatDirection.tsv" \
            -cpuMacSyFinder ${task.cpus} -cpuProkka ${task.cpus} \
            -log -html -levelMin 3 \
            -cas -ccvRep -getSummaryCasfinder -gcode "${gcode}" \
            ${extraArgs} \
            >> "\${task_root}/\${wrapper_log}" 2>&1
        exit_code=\$?
        set -e
        popd >/dev/null

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_ccfinder=%s\n' "\${attempt}" >> "\${wrapper_log}"
        (( attempt += 1 ))
    done

    mkdir -p ccfinder
    result_json_path=\$(find "\${tool_output_root}" "\${run_root}" -type f -name 'result.json' | head -n 1 || true)
    internal_log_path=\$(find "\${tool_output_root}" "\${run_root}" -type f \\( -name 'ccfinder.log' -o -name 'logFile_*' \\) | head -n 1 || true)
    if [[ "\${exit_code}" -eq 0 && -n "\${result_json_path}" ]]; then
        perl -0pi -e 's/:\\s*(?=,|\\}|\\])/: null/g' "\${result_json_path}"
    fi
    if [[ "\${exit_code}" -eq 0 && -n "\${result_json_path}" ]]; then
        cp "\${result_json_path}" result.json
    else
        : > result.json
    fi

    if [[ -n "\${internal_log_path}" ]]; then
        cp "\${internal_log_path}" ccfinder.log
    else
        : > ccfinder.log
    fi
    cat "\${wrapper_log}" >> ccfinder.log
    printf 'exit_code=%s\n' "\$exit_code" >> ccfinder.log

    rm -rf ccfinder
    mkdir -p ccfinder
    if [[ -s result.json ]]; then
        cp result.json ccfinder/
    fi

    ccfinder_version="\$(perl "\${ccfinder_root}/CRISPRCasFinder.pl" -v 2>&1 | sed -n 's/.*version \\([^,[:space:]]*\\).*/\\1/p' | head -n 1 || true)"
    ccfinder_version="\${ccfinder_version:-NA}"
    printf '"%s":\n  crisprcasfinder: "%s"\n' \
      "${task.process}" \
      "\${ccfinder_version}" \
      > versions.yml
    """

    stub:
    '''
    mkdir -p ccfinder
    cat <<'EOF' > result.json
    {"Sequences": []}
    EOF
    cp result.json ccfinder/result.json
    : > ccfinder.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      crisprcasfinder: "stub"
    EOF
    '''
}
