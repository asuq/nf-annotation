/*
 * Run one CheckM2 translation-table prediction per sample. Sample-level tool
 * failures are converted into empty reports so downstream QC summarisation can
 * record them without killing the whole pipeline.
 */
process CHECKM2 {
    tag "${meta.accession} / ttable ${translation_table}"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/checkm2_gcode${translation_table}" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['quality_report.tsv', 'checkm2.log']
                ? filename
                : null
        },
    )

    input:
    tuple val(meta), path(genome), path(checkm2_db)
    val translation_table

    output:
    tuple val(meta), val(translation_table), path("checkm2_gcode${translation_table}"), path('quality_report.tsv'), path('checkm2.log'), emit: results
    tuple val(meta), val(translation_table), path('quality_report.tsv'), emit: quality_report
    path 'versions.yml', emit: versions

    script:
    """
    database_path="${checkm2_db}"
    if [[ ! -d "\${database_path}" ]]; then
        echo "params.checkm2_db must point to a CheckM2 database directory." >&2
        exit 1
    fi

    shopt -s nullglob
    dmnd_candidates=("\${database_path}"/*.dmnd)
    shopt -u nullglob
    if [[ \${#dmnd_candidates[@]} -ne 1 ]]; then
        echo "CHECKM2 expected exactly one .dmnd file in \${database_path}." >&2
        exit 1
    fi
    database_path="\${dmnd_candidates[0]}"

    output_dir="checkm2_gcode${translation_table}"

    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    : > checkm2.log
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> checkm2.log
        rm -rf "\${output_dir}"
        set +e
        checkm2 predict \
            --extension fasta \
            --input "${genome}" \
            --output-directory "\${output_dir}" \
            --database_path "\${database_path}" \
            --threads ${task.cpus} \
            --ttable ${translation_table} \
            --force \
            >> checkm2.log 2>&1
        exit_code=\$?
        set -e

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_checkm2=%s\n' "\${attempt}" >> checkm2.log
        (( attempt += 1 ))
    done

    mkdir -p "\${output_dir}"
    if [[ "\${exit_code}" -eq 0 && -f "\${output_dir}/quality_report.tsv" ]]; then
        cp "\${output_dir}/quality_report.tsv" quality_report.tsv
    else
        : > quality_report.tsv
    fi

    rm -rf "\${output_dir}"
    mkdir -p "\${output_dir}"
    if [[ -s quality_report.tsv ]]; then
        cp quality_report.tsv "\${output_dir}/"
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> checkm2.log

    checkm2_version="\$(command -v checkm2 >/dev/null 2>&1 && checkm2 --version 2>&1 | awk 'NF { print; exit }' || echo NA)"
    printf '"%s":\n  checkm2: "%s"\n' \
      "${task.process}" \
      "\${checkm2_version}" \
      > versions.yml
    """

    stub:
    def stubMeanLength = translation_table.toString() == '4' ? 300 : 150
    def stubCodingSequences = translation_table.toString() == '4' ? 3000 : 5300
    """
    mkdir -p "checkm2_gcode${translation_table}"
    cat <<'EOF' > quality_report.tsv
    Name	Completeness	Contamination	Coding_Density	Average_Gene_Length	Total_Coding_Sequences	Genome_Size	GC_Content	Contig_N50	Translation_Table_Used
    ${meta.internal_id}	95	2	0.9	${stubMeanLength}	${stubCodingSequences}	3000000	0.3	500000	${translation_table}
    EOF
    cp quality_report.tsv "checkm2_gcode${translation_table}/quality_report.tsv"
    : > checkm2.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      checkm2: "stub"
    EOF
    """
}
