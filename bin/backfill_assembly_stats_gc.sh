#!/usr/bin/env bash
set -euo pipefail

usage() {
	cat <<'EOF'
Usage: backfill_assembly_stats_gc.sh --results-dir PATH [--input PATH] [--output PATH] [--seqtk-binary CMD] [--jobs INT]

Backfill gc_content into an existing assembly_stats.tsv for one previous run.
The script preserves the existing table rows and columns, appending gc_content
derived from the published staged genomes under samples/<accession>/staged/.
EOF
}

cleanup_temp_files() {
	if [[ -n "${OUTPUT_TMP}" ]]; then
		rm -f "${OUTPUT_TMP}"
	fi
	if ((${#CLEANUP_DIRS[@]} > 0)); then
		rm -rf "${CLEANUP_DIRS[@]}"
	fi
}

log_error() {
	printf 'ERROR: %s\n' "$1" >&2
}

log_info() {
	printf 'INFO: %s\n' "$1" >&2
}

validate_jobs() {
	local jobs="$1"

	if ! [[ "${jobs}" =~ ^[1-9][0-9]*$ ]]; then
		log_error "--jobs must be a positive integer."
		return 1
	fi
}

resolve_script_path() {
	local invocation="$1"

	if [[ "${invocation}" == */* ]]; then
		printf '%s\n' "${invocation}"
		return 0
	fi

	command -v "${invocation}"
}

find_column_index() {
	local -n header_ref=$1
	local target="$2"
	local index=0
	local normalised_target

	normalised_target="$(printf '%s' "$target" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')"
	for column in "${header_ref[@]}"; do
		if [[ "$(printf '%s' "$column" | tr '[:upper:]' '[:lower:]' | tr -cd '[:alnum:]')" == "${normalised_target}" ]]; then
			printf '%s\n' "${index}"
			return 0
		fi
		index=$((index + 1))
	done

	return 1
}

resolve_staged_fasta() {
	local results_dir="$1"
	local accession="$2"
	local staged_dir fasta_candidates
	local -a fasta_paths

	staged_dir="${results_dir}/samples/${accession}/staged"
	if [[ ! -d "${staged_dir}" ]]; then
		log_error "Missing staged directory for ${accession}: ${staged_dir}"
		return 1
	fi

	shopt -s nullglob
	fasta_paths=("${staged_dir}"/*.fasta)
	shopt -u nullglob

	if ((${#fasta_paths[@]} == 0)); then
		log_error "Expected exactly one staged FASTA for ${accession}, found none in ${staged_dir}."
		return 1
	fi
	if ((${#fasta_paths[@]} > 1)); then
		fasta_candidates="$(printf '%s, ' "${fasta_paths[@]}")"
		fasta_candidates="${fasta_candidates%, }"
		log_error "Expected exactly one staged FASTA for ${accession}, found multiple: ${fasta_candidates}"
		return 1
	fi

	printf '%s\n' "${fasta_paths[0]}"
}

run_seqtk_comp() {
	local fasta_path="$1"

	SEQTK_BINARY="${seqtk_binary}" FASTA_PATH="${fasta_path}" \
		bash -c 'eval "set -- $SEQTK_BINARY"; "$@" comp "$FASTA_PATH"'
}

compute_gc_content() {
	local fasta_path="$1"
	local gc_content

	gc_content="$(
		run_seqtk_comp "${fasta_path}" | awk '
            BEGIN {
                adenine = 0
                cytosine = 0
                guanine = 0
                thymine = 0
            }
            {
                if ($2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ && $6 ~ /^[0-9]+$/) {
                    adenine += $3
                    cytosine += $4
                    guanine += $5
                    thymine += $6
                }
            }
            END {
                canonical = adenine + cytosine + guanine + thymine
                if (canonical == 0) {
                    exit 1
                }
                value = ((cytosine + guanine) / canonical) * 100
                formatted = sprintf("%.6f", value)
                sub(/0+$/, "", formatted)
                sub(/\.$/, "", formatted)
                print formatted
            }
        '
	)"

	if [[ -z "${gc_content}" ]]; then
		return 1
	fi
	printf '%s\n' "${gc_content}"
}

run_worker() {
	local row_index="$1"
	local accession="$2"
	local fasta_path="$3"
	local worker_dir="$4"
	local seqtk_command="$5"
	local result_path error_path gc_content

	result_path="${worker_dir}/results/${row_index}.tsv"
	error_path="${worker_dir}/errors/${row_index}.log"
	seqtk_binary="${seqtk_command}"

	log_info "Starting accession ${accession}: row=${row_index} staged_fasta=${fasta_path}"
	if ! gc_content="$(compute_gc_content "${fasta_path}")"; then
		log_info "Backfill failed for accession ${accession}: row=${row_index}"
		printf '%s\n' "Failed to calculate gc_content for ${accession} from ${fasta_path}." >"${error_path}"
		return 255
	fi

	printf '%s\n' "${gc_content}" >"${result_path}"
	log_info "Backfilled accession ${accession}: row=${row_index} gc_content=${gc_content}"
}

collect_worker_error() {
	local worker_dir="$1"
	local task_count="$2"
	local row_index error_path

	for ((row_index = 1; row_index <= task_count; row_index++)); do
		error_path="${worker_dir}/errors/${row_index}.log"
		if [[ -s "${error_path}" ]]; then
			cat "${error_path}"
			return 0
		fi
	done

	return 1
}

merge_worker_results() {
	local worker_dir="$1"
	local task_count="$2"
	local input_path="$3"
	local output_path="$4"
	local output_tmp row_index result_path input_line

	output_tmp="${worker_dir}/assembly_stats.tsv"
	{
		printf '%s\tgc_content\n' "$(
			IFS=$'\t'
			printf '%s' "${HEADER_COLUMNS[*]}"
		)"
		row_index=0
		while IFS= read -r input_line; do
			input_line="${input_line%$'\r'}"
			if [[ -z "${input_line}" ]]; then
				continue
			fi
			row_index=$((row_index + 1))
			result_path="${worker_dir}/results/${row_index}.tsv"
			if [[ ! -f "${result_path}" ]]; then
				log_error "GC backfill worker did not produce row ${row_index}."
				return 1
			fi
			printf '%s\t%s\n' "${input_line}" "$(cat "${result_path}")"
		done < <(tail -n +2 "${input_path}")
		if ((row_index != task_count)); then
			log_error "GC backfill merge saw ${row_index} data rows, expected ${task_count}."
			return 1
		fi
	} >"${output_tmp}"

	mv "${output_tmp}" "${output_path}"
}

main() {
	local results_dir="" input="" output="" output_dir jobs="1"
	local accession_index required_column
	local duplicate_accessions="" accession staged_fasta input_line
	local worker_dir task_file worker_script task_count
	local -a header_array fields

	seqtk_binary="seqtk"

	if [[ "${1:-}" == "--worker" ]]; then
		shift
		if (($# != 5)); then
			log_error "Worker mode expects row index, accession, FASTA path, worker dir, and seqtk command."
			return 1
		fi
		run_worker "$1" "$2" "$3" "$4" "$5"
		return $?
	fi

	while (($# > 0)); do
		case "$1" in
		--results-dir)
			results_dir="$2"
			shift 2
			;;
		--input)
			input="$2"
			shift 2
			;;
		--output)
			output="$2"
			shift 2
			;;
		--seqtk-binary)
			seqtk_binary="$2"
			shift 2
			;;
		--jobs)
			jobs="$2"
			shift 2
			;;
		-h | --help)
			usage
			return 0
			;;
		*)
			log_error "Unknown argument: $1"
			usage >&2
			return 1
			;;
		esac
	done

	if [[ -z "${results_dir}" ]]; then
		log_error "--results-dir is required."
		usage >&2
		return 1
	fi
	if ! validate_jobs "${jobs}"; then
		return 1
	fi
	if [[ ! -d "${results_dir}" ]]; then
		log_error "Missing results directory: ${results_dir}"
		return 1
	fi

	input="${input:-${results_dir}/cohort/assembly_stats/assembly_stats.tsv}"
	output="${output:-${results_dir}/cohort/assembly_stats/assembly_stats.with_gc.tsv}"

	if [[ ! -f "${input}" ]]; then
		log_error "Missing input assembly stats table: ${input}"
		return 1
	fi

	log_info "Backfill configuration: results_dir=${results_dir} input=${input} output=${output} jobs=${jobs} seqtk_binary=${seqtk_binary}"

	IFS=$'\t' read -r -a header_array < <(head -n 1 "${input}" | tr -d '\r')
	if ((${#header_array[@]} == 0)); then
		log_error "Input assembly stats table is empty: ${input}"
		return 1
	fi

	for required_column in accession n50 scaffolds genome_size; do
		if ! find_column_index header_array "${required_column}" >/dev/null; then
			log_error "Input assembly stats table is missing the ${required_column} column."
			return 1
		fi
	done
	if find_column_index header_array "gc_content" >/dev/null 2>&1; then
		log_error "Input assembly stats table already contains gc_content: ${input}"
		return 1
	fi
	HEADER_COLUMNS=("${header_array[@]}")

	accession_index="$(find_column_index header_array "accession")"
	duplicate_accessions="$(
		awk -F '\t' -v accession_index="$((accession_index + 1))" 'NR > 1 { print $accession_index }' "${input}" |
			sed '/^[[:space:]]*$/d' |
			sort |
			uniq -d
	)"
	if [[ -n "${duplicate_accessions}" ]]; then
		log_error "Input assembly stats table contains duplicate accession values: ${duplicate_accessions}"
		return 1
	fi

	if ! bash -c 'command -v "$1" >/dev/null 2>&1' _ "${seqtk_binary%% *}"; then
		if [[ "${seqtk_binary}" != *" "* ]]; then
			log_error "seqtk is required to backfill gc_content."
			return 1
		fi
	fi

	output_dir="$(dirname "${output}")"
	mkdir -p "${output_dir}"
	OUTPUT_TMP="$(mktemp "${output_dir}/.assembly_stats_gc.XXXXXX")"
	worker_dir="$(mktemp -d)"
	CLEANUP_DIRS+=("${worker_dir}")
	mkdir -p "${worker_dir}/results" "${worker_dir}/errors"
	task_file="${worker_dir}/tasks.nul"
	worker_script="$(resolve_script_path "$0")"
	task_count=0

	while IFS= read -r input_line; do
		input_line="${input_line%$'\r'}"
		if [[ -z "${input_line}" ]]; then
			continue
		fi

		IFS=$'\t' read -r -a fields <<<"${input_line}"
		if ((${#fields[@]} <= accession_index)); then
			log_error "Assembly stats row is missing the accession field."
			return 1
		fi

		accession="${fields[accession_index]}"
		if [[ -z "${accession}" ]]; then
			log_error "Assembly stats row contains an empty accession."
			return 1
		fi

		staged_fasta="$(resolve_staged_fasta "${results_dir}" "${accession}")" || return 1
		task_count=$((task_count + 1))
		log_info "Queued accession ${accession}: row=${task_count} staged_fasta=${staged_fasta}"
		printf '%s\0%s\0%s\0%s\0%s\0' \
			"${task_count}" \
			"${accession}" \
			"${staged_fasta}" \
			"${worker_dir}" \
			"${seqtk_binary}" \
			>>"${task_file}"
	done < <(tail -n +2 "${input}")

	log_info "Detected ${task_count} assembly-stats row(s) for GC backfill."

	if [[ -s "${task_file}" ]]; then
		log_info "Starting GC backfill workers: rows=${task_count} jobs=${jobs}"
		if ! xargs -0 -n 5 -P "${jobs}" bash "${worker_script}" --worker <"${task_file}"; then
			local worker_error
			worker_error="$(collect_worker_error "${worker_dir}" "${task_count}" || true)"
			if [[ -n "${worker_error}" ]]; then
				log_error "${worker_error}"
			else
				log_error "GC backfill worker execution failed."
			fi
			return 1
		fi
	fi

	log_info "Merging GC backfill results into ${output}"
	if ! merge_worker_results "${worker_dir}" "${task_count}" "${input}" "${OUTPUT_TMP}"; then
		return 1
	fi

	log_info "Writing merged GC backfill output to ${output}"
	mv "${OUTPUT_TMP}" "${output}"
	OUTPUT_TMP=""
	log_info "GC backfill complete: rows=${task_count} output=${output}"
}

seqtk_binary=""
OUTPUT_TMP=""
CLEANUP_DIRS=()
HEADER_COLUMNS=()
trap cleanup_temp_files EXIT
main "$@"
