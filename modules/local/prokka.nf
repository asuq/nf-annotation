/*
 * Run Prokka for gcode-qualified samples and emit stable file handles for the
 * downstream protein-bundle validation and annotation.
 */
process PROKKA {
    tag "${meta.accession}"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/prokka" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['prokka.gff', 'prokka.faa', 'prokka.gbk', 'prokka.log']
                ? filename
                : null
        },
    )

    input:
    tuple val(meta), path(genome), val(gcode)

    output:
    tuple val(meta), path('prokka'), path('prokka.gff'), path('prokka.faa'), path('prokka.gbk'), path('prokka.log'), emit: results
    tuple val(meta), path(genome), val(gcode), path('prokka.faa'), path('prokka.gff'), path('prokka.gbk'), path('prokka.log'), emit: bundle_inputs
    path 'versions.yml', emit: versions

    script:
    def internalId = (meta.internal_id ?: meta.accession).toString()
    def rawLocustag = internalId.toUpperCase().replaceAll(/[^A-Z0-9]/, '')
    def locustag = rawLocustag ? rawLocustag.take(20) : 'PROKKA'
    if (!(locustag ==~ /^[A-Z].*/)) {
        locustag = "L${locustag}".take(20)
    }
    """
    export TMPDIR="\$PWD/prokka_tmp"
    mkdir -p "\${TMPDIR}"

    max_attempts="${params.soft_fail_attempts}"
    if [[ "\${max_attempts}" -lt 1 ]]; then
        max_attempts=1
    fi

    attempt=1
    exit_code=1
    : > prokka.log
    while (( attempt <= max_attempts )); do
        printf 'attempt=%s/%s\n' "\${attempt}" "\${max_attempts}" >> prokka.log
        rm -rf prokka
        set +e
        prokka "${genome}" \
            --outdir prokka \
            --prefix "${internalId}" \
            --locustag "${locustag}" \
            --compliant \
            --gcode "${gcode}" \
            --cpus ${task.cpus} \
            --rfam \
            >> prokka.log 2>&1
        exit_code=\$?
        set -e

        if [[ "\${exit_code}" -eq 0 ]]; then
            break
        fi
        if (( attempt == max_attempts )); then
            break
        fi
        printf 'retrying_prokka=%s\n' "\${attempt}" >> prokka.log
        (( attempt += 1 ))
    done

    mkdir -p prokka

    prokka_gff=\$(find prokka -maxdepth 1 -type f -name '*.gff' | head -n 1 || true)
    if [[ "\${exit_code}" -eq 0 && -n "\${prokka_gff}" ]]; then
        cp "\${prokka_gff}" prokka.gff
    else
        : > prokka.gff
    fi

    prokka_faa=\$(find prokka -maxdepth 1 -type f -name '*.faa' | head -n 1 || true)
    if [[ "\${exit_code}" -eq 0 && -n "\${prokka_faa}" ]]; then
        cp "\${prokka_faa}" prokka.faa
    else
        : > prokka.faa
    fi

    prokka_gbk=\$(find prokka -maxdepth 1 -type f -name '*.gbk' | head -n 1 || true)
    if [[ "\${exit_code}" -eq 0 && -n "\${prokka_gbk}" ]]; then
        cp "\${prokka_gbk}" prokka.gbk
    else
        : > prokka.gbk
    fi

    rm -rf prokka
    mkdir -p prokka
    if [[ -s prokka.gff ]]; then
        cp prokka.gff prokka/
    fi
    if [[ -s prokka.faa ]]; then
        cp prokka.faa prokka/
    fi
    if [[ -s prokka.gbk ]]; then
        cp prokka.gbk prokka/
    fi

    printf 'exit_code=%s\n' "\$exit_code" >> prokka.log

    prokka_version="\$(command -v prokka >/dev/null 2>&1 && prokka --version 2>&1 | awk 'NF { value=\$0 } END { if (value) print value }' || true)"
    prokka_version="\${prokka_version:-NA}"
    printf '"%s":\n  prokka: "%s"\n' \
      "${task.process}" \
      "\${prokka_version}" \
      > versions.yml
    """

    stub:
    """
    mkdir -p prokka
    python3 - '${genome}' '${gcode}' <<'PY'
    import sys
    from pathlib import Path
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    records = list(SeqIO.parse(sys.argv[1], 'fasta'))
    code = int(sys.argv[2])
    assert str(records[0].seq).startswith('ATGTGAGCTTAA'), 'Unexpected staged stub sequence'
    end, translation = (12, 'MWA') if code == 4 else (6, 'M')
    for record in records:
        record.annotations['molecule_type'] = 'DNA'
    records[0].features = [SeqFeature(SimpleLocation(0, end, strand=1), type='CDS',
        qualifiers={'locus_tag': ['gene_1'], 'translation': [translation],
                    'transl_table': [str(code)], 'codon_start': ['1']})]
    SeqIO.write(records, 'prokka.gbk', 'genbank')
    Path('prokka.faa').write_text('>gene_1 stub protein\\n' + translation + '\\n')
    with Path('prokka.gff').open('w') as handle:
        handle.write('##gff-version 3\\n')
        handle.write(f'{records[0].id}\\tProkka\\tCDS\\t1\\t{end}\\t.\\t+\\t0\\tID=gene_1;locus_tag=gene_1\\n')
        handle.write('##FASTA\\n')
        SeqIO.write(records, handle, 'fasta')
    PY
    cp prokka.gff prokka/prokka.gff
    cp prokka.faa prokka/prokka.faa
    cp prokka.gbk prokka/prokka.gbk
    printf 'exit_code=0\\n' > prokka.log
    cat <<'EOF' > versions.yml
    "${task.process}":
      prokka: "stub"
    EOF
    """
}
