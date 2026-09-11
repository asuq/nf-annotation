#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"

fail() {
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

assert_file_contains() {
    local file=$1
    local pattern=$2

    grep -F "$pattern" "$file" >/dev/null 2>&1 || fail "expected $file to contain $pattern"
}

assert_path_exists() {
    [ -e "$1" ] || fail "expected path to exist: $1"
}

assert_path_missing() {
    [ ! -e "$1" ] || fail "expected path to be absent: $1"
}

assert_file_contains "$REPO_ROOT/.gitmodules" "path = external/nf-helper"
assert_file_contains "$REPO_ROOT/.gitmodules" "url = https://github.com/asuq/nf-helper.git"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/oist.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/gwdg.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/marmic.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/conf/sites/viper-cpu.config"
assert_path_exists "$REPO_ROOT/external/nf-helper/helpers/cleanup_processed_sample_workdirs.sh"
assert_path_exists "$REPO_ROOT/external/nf-helper/helpers/gwdg_promote_2h_qos.sh"
assert_path_missing "$REPO_ROOT/helpers/cleanup_processed_sra_workdirs.sh"
assert_file_contains "$REPO_ROOT/external/nf-helper/conf/sites/viper-cpu.config" "viper_slurm_queue_size = 250"
assert_file_contains "$REPO_ROOT/external/nf-helper/conf/sites/viper-cpu.config" "queueSize = params.viper_slurm_queue_size as int"

assert_file_contains "$REPO_ROOT/conf/oist.config" "external/nf-helper/conf/sites/oist.config"
assert_file_contains "$REPO_ROOT/conf/gwdg.config" "external/nf-helper/conf/sites/gwdg.config"
assert_file_contains "$REPO_ROOT/conf/marmic.config" "external/nf-helper/conf/sites/marmic.config"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "external/nf-helper/conf/sites/viper-cpu.config"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "'viper-cpu' {"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "DOWNLOAD_BUSCO_DATASET"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "PREP_TAXDUMP_DATABASE"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "PREP_CODETTA_DATABASE"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "DOWNLOAD_CHECKM2_DATABASE"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "DOWNLOAD_BUSCO_DATABASES"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "PREP_ANNOTATION_RESOURCE"
assert_file_contains "$REPO_ROOT/conf/viper-cpu.config" "executor = 'local'"
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig 'conf/marmic.config'"
assert_file_contains "$REPO_ROOT/nextflow.config" "includeConfig 'conf/viper-cpu.config'"

bash -n "$REPO_ROOT/helpers/cleanup_processed_sample_workdirs.sh"
bash -n "$REPO_ROOT/helpers/gwdg_promote_2h_qos.sh"
"$REPO_ROOT/helpers/cleanup_processed_sample_workdirs.sh" --help >/dev/null
"$REPO_ROOT/helpers/gwdg_promote_2h_qos.sh" --help >/dev/null

printf 'nf-helper integration tests passed\n'
