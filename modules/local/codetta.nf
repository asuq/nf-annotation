/*
 * Run Codetta for every staged genome and retain all raw outputs per sample.
 */
process CODETTA {
    tag "${meta.accession}"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(genome), path(codetta_db)

    output:
    tuple val(meta), path('codetta'), emit: results
    tuple val(meta), path('codetta'), emit: summary_input
    path 'versions.yml', emit: versions

    script:
    def prefix = (meta.internal_id ?: meta.accession).toString()
    def inferExtraArgs = (params.codetta_extra_args ?: '').toString()
    """
    codetta_root="\${CODETTA_HOME:-/opt/codetta}"
    if [[ ! -d "\${codetta_root}/resources" ]]; then
        echo "Codetta resources directory not found under \${codetta_root}." >&2
        exit 1
    fi

    align_prefix="codetta/${prefix}"
    inference_output="codetta/codetta_inference.txt"
    results_summary="codetta/codetta_results_summary.csv"
    resource_directory="\$PWD/codetta_resources"

    mkdir -p codetta
    rm -rf "\${resource_directory}"
    mkdir -p "\${resource_directory}"
    cp -R "\${codetta_root}/resources/". "\${resource_directory}/"
    cp -R "${codetta_db}/". "\${resource_directory}/"

    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    : > codetta/codetta.log
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> codetta/codetta.log
        rm -f codetta/*.alignment_output.txt codetta/*.alignment_output.txt.gz codetta/*.genetic_code.out
        rm -rf codetta/*.temp_files
        rm -f "\${inference_output}" "\${results_summary}"

        set +e
        codetta_align.py "${genome}" \
            --align_output "\${align_prefix}" \
            -p Pfam-A_enone.hmm \
            --resource_directory "\${resource_directory}" \
            >> codetta/codetta.log 2>&1
        exit_code=\$?
        if [[ "\${exit_code}" -eq 0 ]]; then
            codetta_summary.py \
                "\${align_prefix}" \
                -p Pfam-A_enone.hmm \
                --resource_directory "\${resource_directory}" \
                >> codetta/codetta.log 2>&1
            exit_code=\$?
        fi
        if [[ "\${exit_code}" -eq 0 ]]; then
            codetta_infer.py \
                "\${align_prefix}" \
                --inference_output "\${inference_output}" \
                --results_summary "\${results_summary}" \
                -p Pfam-A_enone.hmm \
                --resource_directory "\${resource_directory}" \
                ${inferExtraArgs} \
                >> codetta/codetta.log 2>&1
            exit_code=\$?
        fi
        set -e

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_codetta=%s\n' "\${attempt}" >> codetta/codetta.log
        (( attempt += 1 ))
    done

    if [[ ! -f "\${inference_output}" ]]; then
        : > "\${inference_output}"
    fi
    if [[ ! -f "\${results_summary}" ]]; then
        : > "\${results_summary}"
    fi
    printf 'exit_code=%s\n' "\${exit_code}" >> codetta/codetta.log
    if [[ "\${exit_code}" -eq 0 && -s "\${inference_output}" && -s "\${results_summary}" ]]; then
        rm -f codetta/*.alignment_output.txt codetta/*.alignment_output.txt.gz codetta/*.genetic_code.out
        rm -rf codetta/*.temp_files
        rm -rf "\${resource_directory}"
    fi

    codetta_version='NA'
    codetta_source_commit='NA'
    if [[ -f "\${codetta_root}/.codetta_version" ]]; then
        codetta_version="\$(tr -d '\n' < "\${codetta_root}/.codetta_version")"
    fi
    if [[ -f "\${codetta_root}/.codetta_source_commit" ]]; then
        codetta_source_commit="\$(tr -d '\n' < "\${codetta_root}/.codetta_source_commit")"
    fi
    printf '"%s":\n  codetta: "%s"\n  codetta_source_commit: "%s"\n' \
      "${task.process}" \
      "\${codetta_version}" \
      "\${codetta_source_commit}" \
      > versions.yml
    """

    stub:
    '''
    mkdir -p codetta
    cat <<'EOF' > codetta/codetta_inference.txt
    # Analysis arguments
    alignment_prefix   codetta/sample_a
    profile_database   Pfam-A_enone.hmm
    results_summary    codetta/codetta_results_summary.csv
    #
    # Final genetic code inference
    FFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    EOF
    cat <<'EOF' > codetta/codetta_results_summary.csv
    prefix,profile_db,evalue_threshold,prob_threshold,max_fraction,excluded_domains,inferred_gencode,inferred_gencode_best_models
    codetta/sample_a,Pfam-A_enone.hmm,1e-10,0.9999,0.01,,FFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG,FFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    EOF
    cat <<'EOF' > codetta/codetta.log
    Genetic code: FFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    exit_code=0
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      codetta: "v2.0"
      codetta_source_commit: "863359ed326276602d44e48227b6003ac6ffd266"
    EOF
    '''
}
