#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly PIPELINE_WRAPPER="bin/run_pipeline_test.sh"
readonly VALIDATOR="bin/validate_hpc_matrix.py"
readonly ACCEPTANCE_HARNESS="bin/run_acceptance_tests.py"
readonly MEDIUM_SOURCE_CATALOG="${REPO_ROOT}/assets/tables/medium/source_catalog.tsv"
readonly MEDIUM_COHORT_PLAN="${REPO_ROOT}/assets/tables/medium/cohort_plan.tsv"
readonly MEDIUM_SAMPLE_COUNT="20"
readonly DB_READY_MARKER=".nf_myco_ready.json"
readonly VALID_MODES="prepare medium-prepare db-full db-reuse db-matrix p1 p2 all"
readonly MINIMUM_HOST_PYTHON="3.12"

show_usage() {
    cat <<'EOF'
Usage:
  bin/run_oist_hpc_matrix.sh [--dry-run] [--resume] --hpc-root PATH <prepare|medium-prepare|db-full|db-reuse|db-matrix|p1|p2|all> [options]

Run the OIST HPC validation campaign through the medium real-data case.

Required:
  --hpc-root PATH             Campaign root on HPC-native storage.

Options:
  --dry-run                   Print commands and filesystem actions without executing them.
  --resume                    Pass -resume to Nextflow and wrapper runs where supported.
  --gcode-rule RULE           Optional gcode rule: mean_gene_length_ratio (default), strict_delta, delta_then_11.
  --medium-sample-csv PATH    Medium-cohort sample sheet override.
  --medium-metadata PATH      Medium-cohort metadata TSV override.
  --singularity-cache PATH    Override the Singularity cache root.
  -h, --help                  Show this help message.

Modes:
  prepare     Prepare the tracked 9-sample cohort on the login node.
  medium-prepare
              Prepare the fixed medium Mycoplasmatota/Bacillota cohort.
  db-full     Run the full runtime-database prep gate on SLURM.
  db-reuse    Re-run the full runtime-database prep gate and expect reuse.
  db-matrix   Run disposable database existence-state cases.
  p1          Run the tracked 9-sample real pipeline gate.
  p2          Run the medium Mycoplasmatota/Bacillota real-data gate.
  all         Run prepare, db-full, db-reuse, db-matrix, p1, and p2 in order.

Notes:
  - Resource settings come from the coded OIST profile and process defaults.
  - Host python3 must be >= 3.12 for the harness and matrix validators.
  - Runtime database gates cover taxdump, CheckM2, Codetta, BUSCO, and eggNOG.
  - PADLOC is not part of runtime database prep because its fixed image bundles the DB.
EOF
}

print_command() {
    local escaped
    local parts=()
    local arg
    for arg in "$@"; do
        printf -v escaped '%q' "${arg}"
        parts+=("${escaped}")
    done
    printf '%s\n' "${parts[*]}"
}

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

preflight_host_python() {
    local version_text

    if ! command -v python3 >/dev/null 2>&1; then
        fail "python3 was not found on PATH. Load Python ${MINIMUM_HOST_PYTHON} or adjust PATH before running the HPC wrapper."
    fi

    version_text="$(python3 --version 2>&1 | tr -d '\n')"
    if ! python3 -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 12) else 1)' >/dev/null 2>&1; then
        fail \
            "Host python3 must be >= ${MINIMUM_HOST_PYTHON} for the acceptance harness and matrix validators. Current interpreter: ${version_text}. Load Python ${MINIMUM_HOST_PYTHON} or adjust PATH before rerunning."
    fi
}

announce_test() {
    local label="$1"
    local description="$2"

    printf 'INFO: [%s] Running %s\n' "${label}" "${description}"
}

is_valid_mode() {
    case "$1" in
        prepare | medium-prepare | db-full | db-reuse | db-matrix | p1 | p2 | all)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

run_or_print() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        print_command "$@"
        return 0
    fi
    "$@"
}

reset_case_root() {
    local resolved_case_root

    if [[ -z "${CASE_ROOT}" ]]; then
        fail "CASE_ROOT must not be empty"
    fi

    resolved_case_root="$(python3 -c 'import pathlib, sys; print(pathlib.Path(sys.argv[1]).resolve())' "${CASE_ROOT}")"
    if [[ "${resolved_case_root}" == "/" ]]; then
        fail "Refusing to remove root directory for db-matrix reset"
    fi
    if [[ "${resolved_case_root}" != */db_cases ]]; then
        fail "Refusing to reset non-disposable case root: ${resolved_case_root}"
    fi

    if [[ "${DRY_RUN}" == "true" ]]; then
        print_command rm -rf "${resolved_case_root}"
        print_command mkdir -p "${resolved_case_root}"
        return 0
    fi

    rm -rf "${resolved_case_root}"
    mkdir -p "${resolved_case_root}"
}

run_expect_failure() {
    local expected_text="$1"
    local log_path="$2"
    shift 2

    if [[ "${DRY_RUN}" == "true" ]]; then
        print_command "$@"
        return 0
    fi

    mkdir -p "$(dirname "${log_path}")"
    if "$@" >"${log_path}" 2>&1; then
        fail "Command unexpectedly succeeded: $(print_command "$@")"
    fi
    if ! grep -Fq "${expected_text}" "${log_path}"; then
        fail "Failure log ${log_path} does not contain expected text: ${expected_text}"
    fi
}

ensure_medium_inputs() {
    if [[ -n "${MEDIUM_SAMPLE_CSV}" || -n "${MEDIUM_METADATA}" ]]; then
        if [[ -z "${MEDIUM_SAMPLE_CSV}" || -z "${MEDIUM_METADATA}" ]]; then
            fail "Both --medium-sample-csv and --medium-metadata are required together"
        fi
        return 0
    fi

    MEDIUM_SAMPLE_CSV="${MEDIUM_ROOT}/generated/sample_sheet.csv"
    MEDIUM_METADATA="${MEDIUM_ROOT}/generated/metadata.tsv"
    run_medium_prepare
}

run_medium_prepare() {
    announce_test "medium-prepare" "fixed medium Mycoplasmatota/Bacillota cohort prepare"
    run_or_print \
        python3 "${ACCEPTANCE_HARNESS}" prepare \
        --work-root "${MEDIUM_ROOT}" \
        --source-catalog "${MEDIUM_SOURCE_CATALOG}" \
        --cohort-plan "${MEDIUM_COHORT_PLAN}"
}

require_path() {
    local path="$1"
    local description="$2"
    if [[ ! -e "${path}" ]]; then
        fail "Missing ${description}: ${path}"
    fi
}

require_golden_db_tree() {
    require_path "${TAXDUMP_DIR}/${DB_READY_MARKER}" "taxdump ready marker"
    require_path "${CHECKM2_DIR}/${DB_READY_MARKER}" "CheckM2 ready marker"
    require_path "${CODETTA_DIR}/${DB_READY_MARKER}" "Codetta ready marker"
    require_path "${BUSCO_DIR}/${DB_READY_MARKER}" "BUSCO ready marker"
    require_path "${EGGNOG_DIR}/${DB_READY_MARKER}" "eggNOG ready marker"
}

require_tracked_cohort() {
    require_path "${ACCEPT_ROOT}/generated/sample_sheet.csv" "tracked cohort sample sheet"
    require_path "${ACCEPT_ROOT}/generated/metadata.tsv" "tracked cohort metadata"
}

run_prepare() {
    announce_test "prepare" "tracked 9-sample cohort prepare"
    run_or_print "${PIPELINE_WRAPPER}" prepare --work-root "${ACCEPT_ROOT}"
}

set_dbfull_expected_status() {
    local component path
    local -a ready_components=()
    local -a unready_components=()
    local -a db_roots=(
        "taxdump:${TAXDUMP_DIR}"
        "checkm2:${CHECKM2_DIR}"
        "codetta:${CODETTA_DIR}"
        "busco_root:${BUSCO_DIR}"
        "eggnog:${EGGNOG_DIR}"
    )

    for entry in "${db_roots[@]}"; do
        component="${entry%%:*}"
        path="${entry#*:}"
        if [[ -f "${path}/${DB_READY_MARKER}" ]]; then
            ready_components+=("${component}")
        else
            unready_components+=("${component}")
        fi
    done

    if [[ "${#ready_components[@]}" -eq 0 ]]; then
        DBFULL_EXPECTED_STATUS="prepared"
        return 0
    fi

    if [[ "${#ready_components[@]}" -eq "${#db_roots[@]}" ]]; then
        DBFULL_EXPECTED_STATUS="present"
        return 0
    fi

    fail \
        "Inconsistent runtime DB root under ${DB_ROOT}: ready [${ready_components[*]}], " \
        "not ready [${unready_components[*]}]. Clean the root or rebuild it explicitly."
}

run_dbprep_wrapper() {
    local label="$1"
    local description="$2"
    local expected_status="$3"
    local resume_args=()

    announce_test "${label}" "${description}"

    if [[ "${RESUME}" == "true" ]]; then
        resume_args+=(--resume)
    fi

    run_or_print \
        "${PIPELINE_WRAPPER}" dbprep-slurm \
        "${resume_args[@]}" \
        --dbprep-profile oist \
        --work-root "${ACCEPT_ROOT}" \
        --taxdump "${TAXDUMP_DIR}" \
        --checkm2-db "${CHECKM2_DIR}" \
        --codetta-db "${CODETTA_DIR}" \
        --busco-db "${BUSCO_DIR}" \
        --eggnog-db "${EGGNOG_DIR}" \
        --singularity-cache-dir "${SINGULARITY_CACHE}"

    run_or_print \
        python3 "${VALIDATOR}" dbprep \
        --results-dir "${ACCEPT_ROOT}/runs/dbprep-slurm/results" \
        --taxdump "${TAXDUMP_DIR}" \
        --checkm2-db "${CHECKM2_DIR}" \
        --codetta-db "${CODETTA_DIR}" \
        --busco-db "${BUSCO_DIR}" \
        --eggnog-db "${EGGNOG_DIR}" \
        --expected-component taxdump \
        --expected-component checkm2 \
        --expected-component codetta \
        --expected-component busco_root \
        --expected-component eggnog \
        --expected-status "${expected_status}" \
        --expected-arg=--taxdump \
        --expected-arg=--checkm2_db \
        --expected-arg=--codetta_db \
        --expected-arg=--busco_db \
        --expected-arg=--eggnog_db
}

run_partial_prepare_case() {
    local work_dir="$1"
    local outdir="$2"
    shift 2
    local resume_args=()
    if [[ "${RESUME}" == "true" ]]; then
        resume_args+=(-resume)
    fi
    run_or_print \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${work_dir}" \
        "${resume_args[@]}" \
        "$@" \
        --outdir "${outdir}" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
}

seed_valid_without_marker_cases() {
    local root="$1"

    run_or_print mkdir -p "${root}/checkm2" "${root}/busco/bacillota_odb12" \
        "${root}/busco/mycoplasmatota_odb12" "${root}/codetta" "${root}/eggnog"
    run_or_print cp "${CHECKM2_DIR}/"*.dmnd "${root}/checkm2/"
    run_or_print cp "${BUSCO_DIR}/bacillota_odb12/dataset.cfg" "${root}/busco/bacillota_odb12/"
    run_or_print cp "${BUSCO_DIR}/mycoplasmatota_odb12/dataset.cfg" "${root}/busco/mycoplasmatota_odb12/"
    run_or_print cp "${CODETTA_DIR}/Pfam-A_enone.hmm" "${root}/codetta/"
    run_or_print cp "${CODETTA_DIR}/Pfam-A_enone.hmm.h3f" "${root}/codetta/"
    run_or_print cp "${CODETTA_DIR}/Pfam-A_enone.hmm.h3i" "${root}/codetta/"
    run_or_print cp "${CODETTA_DIR}/Pfam-A_enone.hmm.h3m" "${root}/codetta/"
    run_or_print cp "${CODETTA_DIR}/Pfam-A_enone.hmm.h3p" "${root}/codetta/"
    run_or_print cp "${EGGNOG_DIR}/eggnog.db" "${root}/eggnog/"
    run_or_print cp "${EGGNOG_DIR}/eggnog_proteins.dmnd" "${root}/eggnog/"
}

run_db_matrix() {
    local case_root="${CASE_ROOT}"
    local fail_log

    announce_test "db-matrix" "runtime database existence-state matrix"

    if [[ "${DRY_RUN}" != "true" ]]; then
        require_golden_db_tree
    fi

    reset_case_root

    local valid_root="${case_root}/db3_valid_without_marker"
    announce_test "db-matrix" "db3_valid_without_marker/checkm2"
    seed_valid_without_marker_cases "${valid_root}"
    run_partial_prepare_case \
        "${valid_root}/work_checkm2" \
        "${valid_root}/results_checkm2" \
        --checkm2_db "${valid_root}/checkm2" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${valid_root}/results_checkm2" \
        --checkm2-db "${valid_root}/checkm2" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db
    announce_test "db-matrix" "db3_valid_without_marker/codetta"
    fail_log="${valid_root}/codetta_missing_marker.log"
    run_expect_failure \
        "Destination is valid but incomplete for codetta" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${valid_root}/work_codetta" \
        --codetta_db "${valid_root}/codetta" \
        --download_missing_databases true \
        --outdir "${valid_root}/results_codetta" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db3_valid_without_marker/busco"
    run_partial_prepare_case \
        "${valid_root}/work_busco" \
        "${valid_root}/results_busco" \
        --busco_db "${valid_root}/busco" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${valid_root}/results_busco" \
        --busco-db "${valid_root}/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db
    announce_test "db-matrix" "db3_valid_without_marker/eggnog"
    run_partial_prepare_case \
        "${valid_root}/work_eggnog" \
        "${valid_root}/results_eggnog" \
        --eggnog_db "${valid_root}/eggnog" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${valid_root}/results_eggnog" \
        --eggnog-db "${valid_root}/eggnog" \
        --expected-component eggnog \
        --expected-status prepared \
        --expected-arg=--eggnog_db

    announce_test "db-matrix" "db4_taxdump_missing_with_download"
    run_partial_prepare_case \
        "${case_root}/db4_taxdump/work" \
        "${case_root}/db4_taxdump/results" \
        --taxdump "${case_root}/db4_taxdump/db/ncbi_taxdump_20240914" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_taxdump/results" \
        --taxdump "${case_root}/db4_taxdump/db/ncbi_taxdump_20240914" \
        --expected-component taxdump \
        --expected-status prepared \
        --expected-arg=--taxdump
    announce_test "db-matrix" "db4_checkm2_missing_with_download"
    run_partial_prepare_case \
        "${case_root}/db4_checkm2/work" \
        "${case_root}/db4_checkm2/results" \
        --checkm2_db "${case_root}/db4_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_checkm2/results" \
        --checkm2-db "${case_root}/db4_checkm2/db/checkm2/CheckM2_database" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db
    announce_test "db-matrix" "db4_codetta_missing_with_download"
    run_partial_prepare_case \
        "${case_root}/db4_codetta/work" \
        "${case_root}/db4_codetta/results" \
        --codetta_db "${case_root}/db4_codetta/db/codetta/Pfam-A_enone" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_codetta/results" \
        --codetta-db "${case_root}/db4_codetta/db/codetta/Pfam-A_enone" \
        --expected-component codetta \
        --expected-status prepared \
        --expected-arg=--codetta_db
    announce_test "db-matrix" "db4_busco_missing_with_download"
    run_partial_prepare_case \
        "${case_root}/db4_busco/work" \
        "${case_root}/db4_busco/results" \
        --busco_db "${case_root}/db4_busco/db/busco" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_busco/results" \
        --busco-db "${case_root}/db4_busco/db/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db
    announce_test "db-matrix" "db4_eggnog_missing_with_download"
    run_partial_prepare_case \
        "${case_root}/db4_eggnog/work" \
        "${case_root}/db4_eggnog/results" \
        --eggnog_db "${case_root}/db4_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db4_eggnog/results" \
        --eggnog-db "${case_root}/db4_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --expected-component eggnog \
        --expected-status prepared \
        --expected-arg=--eggnog_db

    announce_test "db-matrix" "db5_taxdump_missing_no_download"
    fail_log="${case_root}/db5_taxdump_missing_no_download.log"
    run_expect_failure \
        "No local source was supplied for taxdump, and remote download is disabled." \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_taxdump/work" \
        --taxdump "${case_root}/db5_taxdump/db/ncbi_taxdump_20240914" \
        --outdir "${case_root}/db5_taxdump/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db5_checkm2_missing_no_download"
    fail_log="${case_root}/db5_checkm2_missing_no_download.log"
    run_expect_failure \
        "CheckM2 destination is missing and download is disabled" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_checkm2/work" \
        --checkm2_db "${case_root}/db5_checkm2/db/checkm2/CheckM2_database" \
        --outdir "${case_root}/db5_checkm2/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db5_codetta_missing_no_download"
    fail_log="${case_root}/db5_codetta_missing_no_download.log"
    run_expect_failure \
        "No local source was supplied for codetta, and remote download is disabled." \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_codetta/work" \
        --codetta_db "${case_root}/db5_codetta/db/codetta/Pfam-A_enone" \
        --outdir "${case_root}/db5_codetta/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db5_busco_missing_no_download"
    fail_log="${case_root}/db5_busco_missing_no_download.log"
    run_expect_failure \
        "BUSCO destination is missing and download is disabled" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_busco/work" \
        --busco_db "${case_root}/db5_busco/db/busco" \
        --outdir "${case_root}/db5_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db5_eggnog_missing_no_download"
    fail_log="${case_root}/db5_eggnog_missing_no_download.log"
    run_expect_failure \
        "eggNOG destination is missing and download is disabled" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db5_eggnog/work" \
        --eggnog_db "${case_root}/db5_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --outdir "${case_root}/db5_eggnog/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"

    run_or_print mkdir -p "${case_root}/db6_taxdump/db" "${case_root}/db6_checkm2/db/checkm2" \
        "${case_root}/db6_codetta/db/codetta" "${case_root}/db6_busco/db" "${case_root}/db6_eggnog/db/Eggnog_db"
    run_or_print touch "${case_root}/db6_taxdump/db/ncbi_taxdump_20240914"
    run_or_print touch "${case_root}/db6_checkm2/db/checkm2/CheckM2_database"
    run_or_print touch "${case_root}/db6_codetta/db/codetta/Pfam-A_enone"
    run_or_print touch "${case_root}/db6_busco/db/busco"
    run_or_print touch "${case_root}/db6_eggnog/db/Eggnog_db/Eggnog_Diamond_db"
    announce_test "db-matrix" "db6_taxdump_file_not_directory"
    fail_log="${case_root}/db6_taxdump_file.log"
    run_expect_failure \
        "Destination must be a directory for taxdump" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_taxdump/work" \
        --taxdump "${case_root}/db6_taxdump/db/ncbi_taxdump_20240914" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_taxdump/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db6_checkm2_file_not_directory"
    fail_log="${case_root}/db6_checkm2_file.log"
    run_expect_failure \
        "CheckM2 destination must be a directory" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_checkm2/work" \
        --checkm2_db "${case_root}/db6_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_checkm2/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db6_codetta_file_not_directory"
    fail_log="${case_root}/db6_codetta_file.log"
    run_expect_failure \
        "Destination must be a directory for codetta" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_codetta/work" \
        --codetta_db "${case_root}/db6_codetta/db/codetta/Pfam-A_enone" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_codetta/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db6_busco_file_not_directory"
    fail_log="${case_root}/db6_busco_file.log"
    run_expect_failure \
        "BUSCO destination must be a directory" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_busco/work" \
        --busco_db "${case_root}/db6_busco/db/busco" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db6_eggnog_file_not_directory"
    fail_log="${case_root}/db6_eggnog_file.log"
    run_expect_failure \
        "eggNOG destination must be a directory" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db6_eggnog/work" \
        --eggnog_db "${case_root}/db6_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true \
        --outdir "${case_root}/db6_eggnog/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"

    run_or_print mkdir -p "${case_root}/db7_checkm2/db/checkm2/CheckM2_database"
    run_or_print mkdir -p "${case_root}/db7_codetta/db/codetta/Pfam-A_enone"
    run_or_print mkdir -p "${case_root}/db7_busco/db/busco/bacillota_odb12"
    run_or_print mkdir -p "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db"
    run_or_print touch "${case_root}/db7_checkm2/db/checkm2/CheckM2_database/broken.txt"
    run_or_print touch "${case_root}/db7_codetta/db/codetta/Pfam-A_enone/broken.txt"
    run_or_print touch "${case_root}/db7_busco/db/busco/bacillota_odb12/dataset.cfg"
    run_or_print touch "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db/broken.txt"
    announce_test "db-matrix" "db7_checkm2_invalid_no_force"
    fail_log="${case_root}/db7_checkm2_invalid.log"
    run_expect_failure \
        "CheckM2 destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db7_checkm2/work" \
        --checkm2_db "${case_root}/db7_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true \
        --outdir "${case_root}/db7_checkm2/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db7_codetta_invalid_no_force"
    fail_log="${case_root}/db7_codetta_invalid.log"
    run_expect_failure \
        "Destination is invalid or partially prepared for codetta" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db7_codetta/work" \
        --codetta_db "${case_root}/db7_codetta/db/codetta/Pfam-A_enone" \
        --download_missing_databases true \
        --outdir "${case_root}/db7_codetta/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db7_busco_invalid_no_force"
    fail_log="${case_root}/db7_busco_invalid.log"
    run_expect_failure \
        "BUSCO destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db7_busco/work" \
        --busco_db "${case_root}/db7_busco/db/busco" \
        --download_missing_databases true \
        --outdir "${case_root}/db7_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db7_eggnog_invalid_no_force"
    fail_log="${case_root}/db7_eggnog_invalid.log"
    run_expect_failure \
        "eggNOG destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db7_eggnog/work" \
        --eggnog_db "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true \
        --outdir "${case_root}/db7_eggnog/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"

    announce_test "db-matrix" "db8_checkm2_invalid_with_force"
    run_partial_prepare_case \
        "${case_root}/db7_checkm2/work_force" \
        "${case_root}/db7_checkm2/results_force" \
        --checkm2_db "${case_root}/db7_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db7_checkm2/results_force" \
        --checkm2-db "${case_root}/db7_checkm2/db/checkm2/CheckM2_database" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db
    announce_test "db-matrix" "db8_codetta_invalid_with_force"
    run_partial_prepare_case \
        "${case_root}/db7_codetta/work_force" \
        "${case_root}/db7_codetta/results_force" \
        --codetta_db "${case_root}/db7_codetta/db/codetta/Pfam-A_enone" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db7_codetta/results_force" \
        --codetta-db "${case_root}/db7_codetta/db/codetta/Pfam-A_enone" \
        --expected-component codetta \
        --expected-status prepared \
        --expected-arg=--codetta_db
    announce_test "db-matrix" "db8_busco_invalid_with_force"
    run_partial_prepare_case \
        "${case_root}/db7_busco/work_force" \
        "${case_root}/db7_busco/results_force" \
        --busco_db "${case_root}/db7_busco/db/busco" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db7_busco/results_force" \
        --busco-db "${case_root}/db7_busco/db/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db
    announce_test "db-matrix" "db8_eggnog_invalid_with_force"
    run_partial_prepare_case \
        "${case_root}/db7_eggnog/work_force" \
        "${case_root}/db7_eggnog/results_force" \
        --eggnog_db "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db7_eggnog/results_force" \
        --eggnog-db "${case_root}/db7_eggnog/db/Eggnog_db/Eggnog_Diamond_db" \
        --expected-component eggnog \
        --expected-status prepared \
        --expected-arg=--eggnog_db

    run_or_print mkdir -p "${case_root}/db9_busco/db/busco/bacillota_odb12"
    run_or_print touch "${case_root}/db9_busco/db/busco/bacillota_odb12/dataset.cfg"
    announce_test "db-matrix" "db9_busco_missing_lineage_no_force"
    fail_log="${case_root}/db9_busco_missing_lineage.log"
    run_expect_failure \
        "BUSCO destination exists but is not ready" \
        "${fail_log}" \
        nextflow run prepare_databases.nf -profile oist \
        -work-dir "${case_root}/db9_busco/work" \
        --busco_db "${case_root}/db9_busco/db/busco" \
        --download_missing_databases true \
        --outdir "${case_root}/db9_busco/results" \
        --singularity_cache_dir "${SINGULARITY_CACHE}"
    announce_test "db-matrix" "db9_busco_missing_lineage_with_force"
    run_partial_prepare_case \
        "${case_root}/db9_busco/work_force" \
        "${case_root}/db9_busco/results_force" \
        --busco_db "${case_root}/db9_busco/db/busco" \
        --download_missing_databases true \
        --force_runtime_database_rebuild true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db9_busco/results_force" \
        --busco-db "${case_root}/db9_busco/db/busco" \
        --expected-component busco_root \
        --expected-status prepared \
        --expected-arg=--busco_db

    run_or_print mkdir -p "${case_root}/db10_checkm2/db/checkm2/CheckM2_database/CheckM2_database"
    run_or_print touch "${case_root}/db10_checkm2/db/checkm2/CheckM2_database/CheckM2_database/CheckM2_database.dmnd"
    announce_test "db-matrix" "db10_checkm2_nested_layout"
    run_partial_prepare_case \
        "${case_root}/db10_checkm2/work" \
        "${case_root}/db10_checkm2/results" \
        --checkm2_db "${case_root}/db10_checkm2/db/checkm2/CheckM2_database" \
        --download_missing_databases true
    run_or_print python3 "${VALIDATOR}" dbprep \
        --results-dir "${case_root}/db10_checkm2/results" \
        --checkm2-db "${case_root}/db10_checkm2/db/checkm2/CheckM2_database" \
        --expected-component checkm2 \
        --expected-status prepared \
        --expected-arg=--checkm2_db

    fail_log="${case_root}/db11_readonly_checkout.log"
    announce_test "db-matrix" "db11_readonly_software_checkout"
    if [[ "${DRY_RUN}" == "true" ]]; then
        print_command test ! -e "${REPO_ROOT}/download_eggnog_data.py.patched"
    else
        if [[ -e "${REPO_ROOT}/download_eggnog_data.py.patched" ]]; then
            fail "Unexpected eggNOG patched downloader in project checkout"
        fi
    fi
}

run_real_case() {
    local work_dir="$1"
    local outdir="$2"
    local sample_csv="$3"
    local metadata_tsv="$4"
    local resume_args=()
    local gcode_args=()

    if [[ "${RESUME}" == "true" ]]; then
        resume_args+=(-resume)
    fi
    if [[ -n "${GCODE_RULE}" ]]; then
        gcode_args+=(--gcode_rule "${GCODE_RULE}")
    fi

    run_or_print \
        nextflow run . -profile oist \
        -work-dir "${work_dir}" \
        "${resume_args[@]}" \
        --sample_csv "${sample_csv}" \
        --metadata "${metadata_tsv}" \
        --taxdump "${TAXDUMP_DIR}" \
        --checkm2_db "${CHECKM2_DIR}" \
        --codetta_db "${CODETTA_DIR}" \
        --busco_db "${BUSCO_DIR}" \
        --eggnog_db "${EGGNOG_DIR}" \
        "${gcode_args[@]}" \
        --singularity_cache_dir "${SINGULARITY_CACHE}" \
        --outdir "${outdir}"
}

run_p1() {
    announce_test "p1" "tracked 9-sample edge-case pipeline test"
    if [[ "${DRY_RUN}" != "true" ]]; then
        require_tracked_cohort
        require_golden_db_tree
    fi
    run_real_case \
        "${RESULT_ROOT}/p1/work" \
        "${RESULT_ROOT}/p1/out" \
        "${ACCEPT_ROOT}/generated/sample_sheet.csv" \
        "${ACCEPT_ROOT}/generated/metadata.tsv"
    run_or_print \
        python3 "${VALIDATOR}" tracked-run \
        --outdir "${RESULT_ROOT}/p1/out" \
        --metadata-tsv "${ACCEPT_ROOT}/generated/metadata.tsv" \
        --db-root "${DB_ROOT}"
}

run_p2() {
    local fixed_medium_sample_csv="${MEDIUM_ROOT}/generated/sample_sheet.csv"
    local fixed_medium_metadata="${MEDIUM_ROOT}/generated/metadata.tsv"

    announce_test "p2" "medium Mycoplasmatota/Bacillota pipeline test"
    ensure_medium_inputs
    if [[ "${DRY_RUN}" != "true" ]]; then
        require_golden_db_tree
    fi
    run_real_case \
        "${RESULT_ROOT}/p2/work" \
        "${RESULT_ROOT}/p2/out" \
        "${MEDIUM_SAMPLE_CSV}" \
        "${MEDIUM_METADATA}"

    if [[ "${MEDIUM_SAMPLE_CSV}" == "${fixed_medium_sample_csv}" && "${MEDIUM_METADATA}" == "${fixed_medium_metadata}" ]]; then
        run_or_print \
            python3 "${VALIDATOR}" medium-run \
            --outdir "${RESULT_ROOT}/p2/out" \
            --db-root "${DB_ROOT}" \
            --sample-count "${MEDIUM_SAMPLE_COUNT}" \
            --allowed-phylum Mycoplasmatota \
            --allowed-phylum Bacillota \
            --metadata-tsv "${MEDIUM_METADATA}" \
            --cohort-plan "${MEDIUM_COHORT_PLAN}" \
            --source-catalog "${MEDIUM_SOURCE_CATALOG}"
        return 0
    fi

    run_or_print \
        python3 "${VALIDATOR}" medium-run \
        --outdir "${RESULT_ROOT}/p2/out" \
        --db-root "${DB_ROOT}" \
        --sample-count "${MEDIUM_SAMPLE_COUNT}" \
        --allowed-phylum Mycoplasmatota \
        --allowed-phylum Bacillota
}

main() {
    local mode=""

    DRY_RUN=false
    RESUME=false
    HPC_ROOT=""
    MEDIUM_SAMPLE_CSV=""
    MEDIUM_METADATA=""
    GCODE_RULE=""
    SINGULARITY_CACHE=""

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --resume)
                RESUME=true
                shift
                ;;
            --gcode-rule)
                GCODE_RULE="$2"
                shift 2
                ;;
            --hpc-root)
                HPC_ROOT="$2"
                shift 2
                ;;
            --medium-sample-csv)
                MEDIUM_SAMPLE_CSV="$2"
                shift 2
                ;;
            --medium-metadata)
                MEDIUM_METADATA="$2"
                shift 2
                ;;
            --singularity-cache)
                SINGULARITY_CACHE="$2"
                shift 2
                ;;
            -h|--help)
                show_usage
                return 0
                ;;
            --)
                shift
                break
                ;;
            -*)
                fail "Unknown option: $1"
                ;;
            *)
                if [[ -n "${mode}" ]]; then
                    fail "Multiple modes provided: ${mode} and $1"
                fi
                mode="$1"
                shift
                ;;
        esac
    done

    if [[ -z "${mode}" ]]; then
        show_usage >&2
        return 1
    fi
    if ! is_valid_mode "${mode}"; then
        fail "Unsupported mode ${mode}. Expected one of: ${VALID_MODES}"
    fi
    if [[ -z "${HPC_ROOT}" ]]; then
        fail "--hpc-root is required"
    fi
    if [[ -n "${GCODE_RULE}" ]]; then
        case "${GCODE_RULE}" in
            mean_gene_length_ratio|strict_delta|delta_then_11)
                ;;
            *)
                fail "--gcode-rule must be one of: mean_gene_length_ratio, strict_delta, delta_then_11"
                ;;
        esac
    fi

    preflight_host_python

    MODE="${mode}"
    ACCEPT_ROOT="${HPC_ROOT}/acceptance"
    MEDIUM_ROOT="${HPC_ROOT}/medium"
    DB_ROOT="${HPC_ROOT}/db"
    RESULT_ROOT="${HPC_ROOT}/results"
    CASE_ROOT="${HPC_ROOT}/db_cases"
    if [[ -z "${SINGULARITY_CACHE}" ]]; then
        SINGULARITY_CACHE="${HPC_ROOT}/singularity_cache"
    fi
    TAXDUMP_DIR="${DB_ROOT}/ncbi_taxdump_20240914"
    CHECKM2_DIR="${DB_ROOT}/checkm2/CheckM2_database"
    CODETTA_DIR="${DB_ROOT}/codetta/Pfam-A_enone"
    BUSCO_DIR="${DB_ROOT}/busco"
    EGGNOG_DIR="${DB_ROOT}/Eggnog_db/Eggnog_Diamond_db"

    if [[ "${DRY_RUN}" != "true" ]]; then
        mkdir -p "${ACCEPT_ROOT}" "${MEDIUM_ROOT}" "${DB_ROOT}" "${RESULT_ROOT}" "${CASE_ROOT}" "${SINGULARITY_CACHE}"
    fi

    case "${MODE}" in
        prepare)
            run_prepare
            ;;
        medium-prepare)
            run_medium_prepare
            ;;
        db-full)
            set_dbfull_expected_status
            run_dbprep_wrapper "db-full" "full runtime DB prep gate" "${DBFULL_EXPECTED_STATUS}"
            ;;
        db-reuse)
            run_dbprep_wrapper "db-reuse" "runtime DB prep reuse gate" present
            ;;
        db-matrix)
            run_db_matrix
            ;;
        p1)
            run_p1
            ;;
        p2)
            run_p2
            ;;
        all)
            run_prepare
            run_medium_prepare
            set_dbfull_expected_status
            run_dbprep_wrapper "db-full" "full runtime DB prep gate" "${DBFULL_EXPECTED_STATUS}"
            run_dbprep_wrapper "db-reuse" "runtime DB prep reuse gate" present
            run_db_matrix
            run_p1
            run_p2
            ;;
    esac
}

main "$@"
