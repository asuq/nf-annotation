#!/bin/bash
# Producer-owned storage handling. Uses only Bash and core filesystem utilities,
# which are present in all pinned native images (including PADLOC).
set -euo pipefail
if stat -c '%d %s %b' / >/dev/null 2>&1; then stat_style=gnu; else stat_style=bsd; fi
if find /dev/null -printf '' >/dev/null 2>&1; then find_style=gnu; else find_style=bsd; fi

die() { printf 'ERROR: %s\n' "$*" >&2; exit 75; }
quote_json() {
    local value=$1
    value=${value//\\/\\\\}; value=${value//\"/\\\"}
    value=${value//$'\n'/\\n}; value=${value//$'\r'/\\r}; value=${value//$'\t'/\\t}
    printf '"%s"' "$value"
}
ordinary_dir() { [[ -d $1 && ! -L $1 ]]; }
file_stat() {
    if [[ $stat_style == gnu ]]; then stat -c '%d %s %b' -- "$1"
    else stat -f '%d %z %b' "$1"; fi
}
# All traversals are physical and bounded to the root filesystem. Reject any
# linked, special, or mounted entry before allowing removal/copying.
measure() {
    local root=$1 path device size blocks root_device files=0 dirs=0 apparent=0 allocated=0 listing invalid=0 mount_path rest encoded
    ordinary_dir "$root" || return 1
    if [[ -r /proc/self/mountinfo ]]; then
        encoded=${root// /\\040}
        while read -r device size blocks path mount_path rest; do
            [[ $mount_path != "$encoded" && $mount_path != "$encoded/"* ]] || return 1
        done < /proc/self/mountinfo
    fi
    listing=$(mktemp "$execution/.storage-inventory.XXXXXXXX") || return 1
    if [[ $find_style == gnu ]]; then
        if ! find -P "$root" -xdev -printf '%D %s %b %y\n' > "$listing"; then rm -f "$listing"; return 1; fi
        read -r root_device size blocks < <(file_stat "$root")
        while read -r device size blocks path; do
            if [[ $device != "$root_device" ]]; then invalid=1; break; fi
            case $path in
                d) dirs=$((dirs+1));;
                f) files=$((files+1)); apparent=$((apparent+size));;
                *) invalid=1; break;;
            esac
            allocated=$((allocated+blocks*512))
        done < "$listing"
        rm -f "$listing"
        [[ $invalid == 0 ]] || return 1
        printf '{"apparent_bytes":%s,"allocated_bytes":%s,"files":%s,"directories":%s}' "$apparent" "$allocated" "$files" "$dirs"
        return
    fi
    if ! find -P "$root" -xdev -print0 > "$listing"; then rm -f "$listing"; return 1; fi
    read -r root_device size blocks < <(file_stat "$root")
    while IFS= read -r -d '' path; do
        if [[ -L $path ]]; then invalid=1; break; fi
        if ! read -r device size blocks < <(file_stat "$path"); then invalid=1; break; fi
        if [[ $device != "$root_device" ]]; then invalid=1; break; fi
        if [[ -d $path ]]; then dirs=$((dirs+1))
        elif [[ -f $path ]]; then files=$((files+1)); apparent=$((apparent+size))
        else invalid=1; break; fi
        allocated=$((allocated+blocks*512))
    done < "$listing"
    rm -f "$listing"
    [[ $invalid == 0 ]] || return 1
    printf '{"apparent_bytes":%s,"allocated_bytes":%s,"files":%s,"directories":%s}' "$apparent" "$allocated" "$files" "$dirs"
}
equal_tree() {
    local source=$1 target=$2 path relative source_count=0 target_count=0
    measure "$source" >/dev/null && measure "$target" >/dev/null || return 1
    while IFS= read -r -d '' path; do
        relative=${path#"$source"/}; [[ $path != "$source" ]] || continue
        source_count=$((source_count+1))
        if [[ -d $path ]]; then ordinary_dir "$target/$relative" || return 1
        else [[ -f $target/$relative && ! -L $target/$relative ]] && cmp -s "$path" "$target/$relative" || return 1; fi
    done < <(find -P "$source" -xdev -print0)
    while IFS= read -r -d '' path; do
        [[ $path != "$target" ]] && target_count=$((target_count+1))
    done < <(find -P "$target" -xdev -print0)
    [[ $source_count == "$target_count" ]]
}
export_raw() {
    local target=$durable/raw temporary
    if [[ $execution -ef $durable ]]; then failure_export=already_durable; return; fi
    if [[ -e $target || -L $target ]]; then
        equal_tree "$raw" "$target" || return 1
    else
        temporary=$(mktemp -d "$durable/.annotation-export.XXXXXXXX") || return 1
        if cp -R "$raw" "$temporary/raw" && equal_tree "$raw" "$temporary/raw" && mv "$temporary/raw" "$target"; then
            rmdir "$temporary" || true
        else
            # Retain an incomplete export for diagnosis; never advertise it as raw.
            return 1
        fi
    fi
    failure_export=exported
}

command=${1:-}
[[ $command == finalize || $command == cleanup-legacy ]] || die 'Invalid storage command'
shift
execution= durable= key= tool= batch=NA native_exit= policy=
while [[ $# -gt 0 ]]; do
    [[ $# -ge 2 ]] || die 'Missing option value'
    case $1 in
        --execution-root) execution=$2;; --durable-work) durable=$2;;
        --task-directory) key=$2;; --tool) tool=$2;; --batch-id) batch=$2;;
        --native-exit) native_exit=$2;; --policy) policy=$2;;
        *) die "Unknown option: $1";;
    esac
    shift 2
done
[[ $policy == success || $policy == off ]] || die 'Invalid cleanup policy'
if [[ $command == cleanup-legacy ]]; then
    ordinary_dir "$execution" || die 'Invalid execution directory'
    execution=$(cd "$execution" && pwd -P)
    case $tool in
        ccfinder) roots=(ccfinder_raw ccfinder_run ccfinder_tmp);;
        prokka) roots=(prokka_tmp);;
        *) die 'Invalid legacy tool';;
    esac
    report=$execution/${tool}_storage.json
    [[ ! -L $report && ! -L $report.tmp ]] || die 'Linked legacy storage report'
    if [[ $tool == ccfinder ]]; then
        [[ -f $execution/ccfinder_generated_contigs.txt && ! -L $execution/ccfinder_generated_contigs.txt ]] || die 'Missing or linked generated-contig inventory'
    fi
    legacy_status=disabled
    if [[ $policy == success ]]; then
        legacy_status=removed
        for name in "${roots[@]}"; do
            if [[ -e $execution/$name || -L $execution/$name ]]; then
                if ! measure "$execution/$name" >/dev/null || ! rm -rf -- "$execution/$name"; then
                    legacy_status=failed
                    printf 'WARNING: Refused or failed cleanup: %s\n' "$name" >&2
                fi
            fi
        done
        if [[ $tool == ccfinder ]]; then
            while IFS= read -r name; do
                if [[ $name == */* || $name != *.fna || -L $execution/$name || ! -f $execution/$name ]]; then
                    legacy_status=failed
                elif ! rm -f -- "$execution/$name"; then legacy_status=failed; fi
            done < "$execution/ccfinder_generated_contigs.txt"
        fi
    fi
    printf '{"policy":"%s","cleanup_status":"%s"}\n' "$policy" "$legacy_status" > "$report.tmp"
    mv "$report.tmp" "$report"
    exit 0
fi
[[ $native_exit =~ ^[0-9]+$ && -n $key ]] || die 'Invalid execution identity/status'
ordinary_dir "$execution" && ordinary_dir "$durable" || die 'Execution and durable roots must be ordinary directories'
execution=$(cd "$execution" && pwd -P); durable=$(cd "$durable" && pwd -P)
raw=$execution/raw; scratch=$execution/scratch
raw_bytes=null; before=null; after=null; cleanup_status=not_applicable; failure_export=not_applicable
cleanup_error=; recorded=null; diagnostics_complete=false
empty='{"apparent_bytes":0,"allocated_bytes":0,"files":0,"directories":0}'
if [[ -f $raw/exit_code.txt && ! -L $raw/exit_code.txt ]]; then
    read -r recorded < "$raw/exit_code.txt" || true
    [[ $recorded =~ ^[0-9]+$ ]] || recorded=null
fi
[[ $recorded != "$native_exit" ]] || diagnostics_complete=true
if [[ $native_exit == 0 && $recorded != 0 ]]; then die 'Native exit status disagrees with retained evidence'; fi
if raw_bytes=$(measure "$raw"); then
    if [[ $native_exit != 0 ]]; then
        cleanup_status=skipped_native_failure
        export_raw || failure_export=failed
    else
        required=(tool.log versions.txt)
        case $tool in
            kofam) required+=(kofam.tsv);;
            cogclassifier) required+=(rpsblast.tsv cogclassifier.native.tsv);;
            eggnog) required+=(execution.json);;
            padloc) required+=(input.domtblout);;
            pfam) required+=(evidence);;
            *) die "Unknown annotation tool: $tool";;
        esac
        for name in "${required[@]}"; do
            [[ -e $raw/$name && ! -L $raw/$name ]] || cleanup_error="Required retained artifact missing: $name"
        done
        if [[ -e $scratch || -L $scratch ]]; then
            before=$(measure "$scratch") || { before=null; cleanup_error='Unsafe or unreadable scratch tree'; }
            if [[ -z $cleanup_error ]]; then
                if [[ $policy == off ]]; then cleanup_status=disabled; after=$before
                elif rm -rf -- "$scratch"; then cleanup_status=removed; after=$empty
                else cleanup_error='Unable to remove owned scratch'; after=$(measure "$scratch") || after=null; fi
            fi
        else before=$empty; after=$empty; cleanup_status=absent; fi
        if [[ -n $cleanup_error ]]; then cleanup_status=failed; printf 'WARNING: %s\n' "$cleanup_error" >&2; fi
    fi
else
    raw_bytes=null; diagnostics_complete=false
    if [[ $native_exit == 0 ]]; then die 'Unsafe or missing retained raw evidence'; fi
    cleanup_status=skipped_native_failure; failure_export=failed
fi
report=$execution/annotation_storage.json
[[ ! -L $report && ! -L $report.tmp ]] || die 'Linked storage report'
{
    printf '{"policy_version":"annotation-storage-v1","policy":'; quote_json "$policy"
    printf ',"task_directory":'; quote_json "$key"
    printf ',"tool":'; quote_json "$tool"
    printf ',"batch_id":'; quote_json "$batch"
    printf ',"native_exit_status":%s,"recorded_native_exit_status":%s,"diagnostics_complete":%s' "$native_exit" "$recorded" "$diagnostics_complete"
    printf ',"raw":%s,"scratch_before":%s,"scratch_after":%s' "$raw_bytes" "$before" "$after"
    printf ',"cleanup_status":'; quote_json "$cleanup_status"
    printf ',"cleanup_error":'; quote_json "$cleanup_error"
    printf ',"failure_export":'; quote_json "$failure_export"
    printf '}\n'
} > "$report.tmp"
mv "$report.tmp" "$report"
if [[ $native_exit != 0 && ! $execution -ef $durable ]]; then
    [[ ! -L $durable/annotation_storage.json ]] || die 'Linked durable report'
    cp "$report" "$durable/annotation_storage.json" || die 'Unable to export diagnostic status'
fi
