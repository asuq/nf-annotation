/*
 * Run all-vs-all FastANI from the prepared path list and normalise the matrix
 * output to a stable `fastani.matrix` filename.
 */
process FASTANI {
    tag "fastani"
    label 'process_high'
    publishDir(
        { "${params.outdir}/cohort/fastani" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    path fastani_inputs
    path fastani_paths

    output:
    path 'fastani.matrix', emit: matrix
    path 'fastani.tsv', emit: raw_output
    path 'fastani.log', emit: log
    path 'versions.yml', emit: versions

    script:
    """
    if [[ ! -s "${fastani_paths}" ]]; then
        : > fastani.matrix
        : > fastani.tsv
        echo 'No ANI-eligible genomes; FastANI was not run.' > fastani.log
    else
        set +e
        fastANI \
            --rl "${fastani_paths}" \
            --ql "${fastani_paths}" \
            --matrix \
            -t ${task.cpus} \
            -o fastani.tsv \
            > fastani.log 2>&1
        exit_code=\$?
        set -e

        matrix_file=\$(find . -maxdepth 1 -type f \\( -name 'fastani.tsv.matrix' -o -name '*.matrix' \\) | head -n 1 || true)
        if [[ -n "\${matrix_file}" ]]; then
            cp "\${matrix_file}" fastani.matrix
        else
            : > fastani.matrix
        fi

        if [[ ! -f fastani.tsv ]]; then
            : > fastani.tsv
        fi

        if grep -q 'Could not open ' fastani.log; then
            echo "FastANI could not open one or more staged inputs." >&2
            sed -n '1,120p' fastani.log >&2
            exit 1
        fi

        if [[ -s "${fastani_paths}" && ! -s fastani.matrix ]]; then
            echo "FastANI did not produce a matrix for a non-empty input list." >&2
            sed -n '1,120p' fastani.log >&2
            exit 1
        fi
        if [[ "\${exit_code}" -ne 0 ]]; then
            echo "FastANI failed with exit code \${exit_code}; see fastani.log." >&2
            exit "\${exit_code}"
        fi
    fi

    cat <<EOF > versions.yml
    "${task.process}":
      fastani: "\$(command -v fastANI >/dev/null 2>&1 && fastANI --version 2>&1 | head -n 1 || echo NA)"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > fastani.matrix
    1
    fastani_inputs/sample_a.fasta
    EOF
    : > fastani.tsv
    : > fastani.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      fastani: "stub"
    EOF
    '''
}
