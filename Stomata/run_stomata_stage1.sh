#!/usr/bin/env bash
set -euo pipefail

# GPU stage-1 benchmark matching run_com_final.sh and run_com_final_chm13.sh.
# Usage:
#   ./run_stomata_stage1.sh [hg38|chm13|all]
#
# Environment overrides:
#   STOMATA_DIR=/path/to/stomata
#   STOMATA_BIN=/path/to/stomata-binary
#   GUIDE_DIR=/path/to/benchmark-guides
#   HG38_FORWARD_DIR=/path/to/hg38-split
#   HG38_REVERSE_DIR=/path/to/hg38-split-rc
#   CHM13_FORWARD_DIR=/path/to/chm13-split
#   CHM13_REVERSE_DIR=/path/to/chm13-split-rc
#   STOMATA_OUTPUT_PARENT=/path/to/results
#   BUILD_JOBS=16
#   KEEP_RAW=1
#   SLEEP_BETWEEN_RUNS=5
#   STOMATA_OUTPUT_MODE=direct|bed

TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STOMATA_DIR="${STOMATA_DIR:-${TOOL_DIR}}"
STOMATA_BIN="${STOMATA_BIN:-${STOMATA_DIR}/build/src/stomata}"
GUIDE_DIR="${GUIDE_DIR:-${TOOL_DIR}/benchmark_guides}"
STOMATA_OUTPUT_PARENT="${STOMATA_OUTPUT_PARENT:-${TOOL_DIR}}"
HG38_FORWARD_DIR="${HG38_FORWARD_DIR:-}"
HG38_REVERSE_DIR="${HG38_REVERSE_DIR:-}"
CHM13_FORWARD_DIR="${CHM13_FORWARD_DIR:-}"
CHM13_REVERSE_DIR="${CHM13_REVERSE_DIR:-}"
BUILD_JOBS="${BUILD_JOBS:-$(nproc)}"
KEEP_RAW="${KEEP_RAW:-0}"
SLEEP_BETWEEN_RUNS="${SLEEP_BETWEEN_RUNS:-5}"
STOMATA_OUTPUT_MODE="${STOMATA_OUTPUT_MODE:-direct}"
STOMATA_UPSTREAM_COMMIT="f63d2cc536fe166728e4108cb4d57537c06dac56"
DATASET="${1:-all}"
MAX_DISTANCE=6
PAD_LENGTH=$((MAX_DISTANCE + 1))
PAD_SEQUENCE="$(printf '%*s' "$PAD_LENGTH" '' | tr ' ' 'N')"

GUIDE_1="${GUIDE_DIR}/new_tested_sgRNA_1.txt"
GUIDE_16="${GUIDE_DIR}/new_tested_sgRNA_16.txt"
GUIDE_32="${GUIDE_DIR}/new_tested_sgRNA_32.txt"
GUIDE_64="${GUIDE_DIR}/new_tested_sgRNA_64.txt"
GUIDE_128="${GUIDE_DIR}/new_tested_sgRNA_128.txt"

declare -A GUIDE_FILES=(
    ["1"]="$GUIDE_1"
    ["16"]="$GUIDE_16"
    ["32"]="$GUIDE_32"
    ["64"]="$GUIDE_64"
    ["128"]="$GUIDE_128"
)

die() {
    echo "Error: $*" >&2
    exit 1
}

elapsed_seconds() {
    local start_ns="$1"
    local end_ns="$2"
    awk -v start="$start_ns" -v end="$end_ns" 'BEGIN {printf "%.6f", (end - start) / 1000000000}'
}

configure_cuda_path() {
    if command -v nvcc >/dev/null; then
        return
    fi

    local nvcc_path=""
    if [[ -x /usr/local/cuda/bin/nvcc ]]; then
        nvcc_path="/usr/local/cuda/bin/nvcc"
    else
        nvcc_path="$(find /usr/local -maxdepth 3 -type f -path '*/cuda-*/bin/nvcc' \
            -executable -print 2>/dev/null | sort -V | tail -n 1)"
    fi

    if [[ -n "$nvcc_path" ]]; then
        export PATH="$(dirname "$nvcc_path"):${PATH}"
        export CUDACXX="$nvcc_path"
        echo "Using CUDA compiler: $nvcc_path"
    fi
}

require_inputs() {
    configure_cuda_path
    command -v cmake >/dev/null || die "cmake is required"
    command -v nvcc >/dev/null ||
        die "CUDA toolkit compiler nvcc was not found. A GPU driver alone is not enough; install CUDA 11.8+ or add its bin directory (for example /usr/local/cuda/bin) to PATH."
    command -v nvidia-smi >/dev/null || die "nvidia-smi is required"
    command -v /usr/bin/time >/dev/null || die "GNU /usr/bin/time is required"
    nvidia-smi -L >/dev/null || die "no working NVIDIA GPU was detected"
    [[ -d "$STOMATA_DIR" ]] || die "Stomata source directory not found: $STOMATA_DIR"
    [[ "$STOMATA_OUTPUT_MODE" == "direct" || "$STOMATA_OUTPUT_MODE" == "bed" ]] ||
        die "STOMATA_OUTPUT_MODE must be direct or bed"

    local count
    for count in 1 16 32 64 128; do
        [[ -f "${GUIDE_FILES[$count]}" ]] ||
            die "guide file not found: ${GUIDE_FILES[$count]}"
    done
}

build_stomata() {
    local cuda_architectures="${CUDA_ARCHITECTURES:-}"
    if [[ -z "$cuda_architectures" ]]; then
        cuda_architectures="$(nvidia-smi --query-gpu=compute_cap \
            --format=csv,noheader | head -n 1 | tr -d '.[:space:]')"
    fi
    [[ "$cuda_architectures" =~ ^[0-9]+$ ]] ||
        die "could not determine the GPU compute architecture; set CUDA_ARCHITECTURES explicitly"

    echo "Configuring and building Stomata..."
    echo "CUDA architecture: ${cuda_architectures}"
    cmake -S "$STOMATA_DIR" -B "${STOMATA_DIR}/build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CUDA_ARCHITECTURES="$cuda_architectures"
    cmake --build "${STOMATA_DIR}/build" -j "$BUILD_JOBS"
    [[ -x "$STOMATA_BIN" ]] || die "Stomata binary not found: $STOMATA_BIN"
}

write_environment_metadata() {
    local output_file="$1"
    {
        echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "hostname=$(hostname)"
        echo "stomata_dir=$STOMATA_DIR"
        echo "stomata_binary=$STOMATA_BIN"
        echo "stomata_upstream_commit=$STOMATA_UPSTREAM_COMMIT"
        echo "stomata_source=AutoFFinder adapted source"
        echo "max_distance=$MAX_DISTANCE"
        echo "padding_length=$PAD_LENGTH"
        echo "distance_mode=levenshtein"
        echo "strand_mode=plus_on_separate_forward_and_reverse_complement_indexes"
        echo "index_layout=two_indexes_to_keep_each_gpu_search_below_uint32_position_limit"
        echo "deduplication=disabled"
        echo "alignment_annotation=disabled"
        echo "scores=disabled"
        echo "output_mode=$STOMATA_OUTPUT_MODE"
        echo
        nvidia-smi --query-gpu=name,uuid,driver_version,memory.total,compute_cap \
            --format=csv,noheader
        echo
        nvcc --version
        echo
        "$STOMATA_BIN" --version
        echo
        echo "AutoFFinder-specific changes are documented in LOCAL_CHANGES.md."
    } > "$output_file"
}

append_fasta_record() {
    local sequence_file="$1"
    local record_name="$2"
    local base_name="$3"
    local strand_suffix="$4"
    local fasta_file="$5"
    local map_file="$6"

    awk \
        -v record="$record_name" \
        -v base="$base_name" \
        -v strand="$strand_suffix" \
        -v padding="$PAD_SEQUENCE" \
        -v map_file="$map_file" '
        BEGIN {
            print ">" record
        }
        {
            sequence = toupper($0)
            gsub(/[[:space:]]/, "", sequence)
            if (sequence == "") {
                next
            }
            if (sequence !~ /^[ACGTN]+$/) {
                printf "Invalid sequence in %s at line %d\n", FILENAME, FNR > "/dev/stderr"
                invalid = 1
                exit 2
            }
            print sequence
            sequence_length += length(sequence)
        }
        END {
            if (!invalid) {
                print padding
                print record "\t" base "\t" strand "\t" sequence_length >> map_file
            }
        }
    ' "$sequence_file" >> "$fasta_file"
}

prepare_strand_genome() {
    local sequence_dir="$1"
    local strand_suffix="$2"
    local fasta_file="$3"
    local map_file="$4"

    [[ -d "$sequence_dir" ]] || die "chromosome directory not found: $sequence_dir"
    : > "$fasta_file"

    local sequence_file filename base
    while IFS= read -r sequence_file; do
        filename="$(basename "$sequence_file")"
        base="${filename%.*}"
        append_fasta_record \
            "$sequence_file" "${base}_${strand_suffix}" "$base" "$strand_suffix" \
            "$fasta_file" "$map_file"
    done < <(printf '%s\n' "${sequence_dir}"/*.txt | sort -V)

    [[ -s "$fasta_file" ]] || die "strand FASTA is empty: $fasta_file"
}

prepare_spacer_file() {
    local guide_file="$1"
    local expected_count="$2"
    local output_file="$3"

    awk -v expected="$expected_count" '
        {
            sequence = toupper($0)
            gsub(/[[:space:]]/, "", sequence)
            if (sequence == "") {
                next
            }
            if (sequence !~ /^[ACGTN]+$/) {
                printf "Invalid guide in %s at line %d\n", FILENAME, FNR > "/dev/stderr"
                invalid = 1
                exit 2
            }
            print count "\t" sequence
            count++
        }
        END {
            if (!invalid && count != expected) {
                printf "Expected %d guides in %s, found %d\n", expected, FILENAME, count > "/dev/stderr"
                exit 3
            }
        }
    ' "$guide_file" > "$output_file"
}

adapt_bed_to_candidates() {
    local forward_bed="$1"
    local reverse_bed="$2"
    local chromosome_map="$3"
    local candidate_dir="$4"
    local expected_guides="$5"
    local stats_file="$6"

    mkdir -p "$candidate_dir"

    local record base strand length
    while IFS=$'\t' read -r record base strand length; do
        : > "${candidate_dir}/${base}_${strand}.txt"
    done < "$chromosome_map"

    # Stomata sorts each spacer's SearchHit vector before batch BED formatting.
    # Batch output is emitted one spacer at a time, so positions are already
    # monotonic for each chromosome/spacer pair; no global sort is needed.
    awk \
        -F '\t' \
        -v map_file="$chromosome_map" \
        -v output_dir="$candidate_dir" \
        -v expected_guides="$expected_guides" \
        -v stats_file="$stats_file" '
        BEGIN {
            while ((getline line < map_file) > 0) {
                split(line, fields, "\t")
                record = fields[1]
                base[record] = fields[2]
                strand[record] = fields[3]
                chromosome_length[record] = fields[4] + 0
            }
            close(map_file)
        }
        {
            raw_hits++
            record = $1
            end = $3
            spacer_name = $7

            if (spacer_name ~ /^[0-9]+$/) {
                target_id = spacer_name + 0
            } else if (spacer_name ~ /^spacer_[0-9]+$/) {
                target_id = spacer_name
                sub(/^spacer_/, "", target_id)
                target_id = target_id - 1
            } else {
                printf "Malformed Stomata spacer name at BED row %d: %s\n", \
                    FNR, spacer_name > "/dev/stderr"
                invalid = 1
                exit 2
            }

            if (!(record in chromosome_length)) {
                printf "Unknown chromosome record in Stomata output: %s\n", record > "/dev/stderr"
                invalid = 1
                exit 2
            }
            if (end !~ /^[0-9]+$/ || target_id < 0 ||
                target_id >= expected_guides) {
                printf "Malformed Stomata BED row %d\n", FNR > "/dev/stderr"
                invalid = 1
                exit 2
            }

            # Discard endpoints inside the synthetic N padding.
            if (end + 0 > chromosome_length[record]) {
                padding_hits++
                next
            }

            key = record SUBSEP target_id
            if (seen[key] && end + 0 < last_end[key]) {
                printf "Stomata endpoints are not sorted at BED row %d\n", FNR > "/dev/stderr"
                invalid = 1
                exit 2
            }
            if (seen[key] && end + 0 == last_end[key]) {
                duplicate_hits++
                next
            }

            output_file = output_dir "/" base[record] "_" strand[record] ".txt"
            print (end + 0) ":" (target_id + 0) >> output_file
            last_end[key] = end + 0
            seen[key] = 1
            candidate_hits++
        }
        END {
            if (!invalid) {
                print raw_hits + 0 "\t" candidate_hits + 0 "\t" \
                    padding_hits + 0 "\t" duplicate_hits + 0 > stats_file
            }
        }
    ' "$forward_bed" "$reverse_bed"
}

run_stomata_once() {
    local dataset="$1"
    local figure="$2"
    local guide_count="$3"
    local distance="$4"
    local spacer_file="$5"
    local forward_index="$6"
    local reverse_index="$7"
    local chromosome_map="$8"
    local candidate_dir="$9"
    local raw_output_base="${10}"
    local log_file_base="${11}"
    local summary_csv="${12}"
    local reported_candidate_dir="$candidate_dir"

    local forward_raw="${raw_output_base}_fw.bed"
    local reverse_raw="${raw_output_base}_rc.bed"
    local forward_log="${log_file_base}_fw.log"
    local reverse_log="${log_file_base}_rc.log"
    local forward_metrics="${raw_output_base}_fw.time"
    local reverse_metrics="${raw_output_base}_rc.time"
    local adapter_stats="${raw_output_base}.adapter"
    local commit="${STOMATA_UPSTREAM_COMMIT:0:12}"

    local -a forward_output_args reverse_output_args
    if [[ "$STOMATA_OUTPUT_MODE" == "direct" ]]; then
        mkdir -p "$candidate_dir"
        forward_output_args=(
            --autoofftarget-output-dir "$candidate_dir"
            --autoofftarget-trim-trailing "$PAD_LENGTH"
        )
        reverse_output_args=("${forward_output_args[@]}")
    else
        forward_output_args=(--format bed --output "$forward_raw")
        reverse_output_args=(--format bed --output "$reverse_raw")
    fi

    echo "Running Stomata: dataset=$dataset figure=$figure guides=$guide_count distance=$distance"

    /usr/bin/time \
        -f '%e,%U,%S,%P,%M' \
        -o "$forward_metrics" \
        "$STOMATA_BIN" \
            --genome "$forward_index" \
            --spacer-file "$spacer_file" \
            --threshold "$distance" \
            --distance-mode levenshtein \
            --strand plus \
            --no-deduplicate \
            --no-compute-mismatches \
            --no-scores \
            --no-treat-u-as-t \
            --verbose \
            "${forward_output_args[@]}" \
            > /dev/null 2> "$forward_log"

    grep -Fq "[GPU: yes]" "$forward_log" ||
        die "Stomata did not report GPU execution; see $forward_log"

    /usr/bin/time \
        -f '%e,%U,%S,%P,%M' \
        -o "$reverse_metrics" \
        "$STOMATA_BIN" \
            --genome "$reverse_index" \
            --spacer-file "$spacer_file" \
            --threshold "$distance" \
            --distance-mode levenshtein \
            --strand plus \
            --no-deduplicate \
            --no-compute-mismatches \
            --no-scores \
            --no-treat-u-as-t \
            --verbose \
            "${reverse_output_args[@]}" \
            > /dev/null 2> "$reverse_log"

    grep -Fq "[GPU: yes]" "$reverse_log" ||
        die "Stomata did not report GPU execution; see $reverse_log"

    local forward_elapsed forward_user forward_system forward_cpu forward_rss
    local reverse_elapsed reverse_user reverse_system reverse_cpu reverse_rss
    IFS=',' read -r forward_elapsed forward_user forward_system forward_cpu forward_rss < "$forward_metrics"
    IFS=',' read -r reverse_elapsed reverse_user reverse_system reverse_cpu reverse_rss < "$reverse_metrics"

    local elapsed user_seconds system_seconds cpu_percent max_rss_kb
    elapsed="$(awk -v f="$forward_elapsed" -v r="$reverse_elapsed" \
        'BEGIN {printf "%.6f", f + r}')"
    user_seconds="$(awk -v f="$forward_user" -v r="$reverse_user" \
        'BEGIN {printf "%.6f", f + r}')"
    system_seconds="$(awk -v f="$forward_system" -v r="$reverse_system" \
        'BEGIN {printf "%.6f", f + r}')"
    cpu_percent="$(awk -v u="$user_seconds" -v s="$system_seconds" -v e="$elapsed" \
        'BEGIN {if (e > 0) printf "%.2f", 100 * (u + s) / e; else print 0}')"
    max_rss_kb="$(awk -v f="$forward_rss" -v r="$reverse_rss" \
        'BEGIN {print (f > r) ? f : r}')"

    local adapter_seconds=0
    local candidate_write_seconds=0
    local raw_hits candidate_hits padding_hits duplicate_hits raw_bytes
    if [[ "$STOMATA_OUTPUT_MODE" == "direct" ]]; then
        local forward_write_ms reverse_write_ms
        local forward_candidates reverse_candidates
        local forward_padding reverse_padding forward_duplicates reverse_duplicates
        forward_write_ms="$(awk '/^AutoOffTarget write:/ {value=$3} END {print value + 0}' "$forward_log")"
        reverse_write_ms="$(awk '/^AutoOffTarget write:/ {value=$3} END {print value + 0}' "$reverse_log")"
        candidate_write_seconds="$(awk -v f="$forward_write_ms" -v r="$reverse_write_ms" \
            'BEGIN {printf "%.6f", (f + r) / 1000}')"

        forward_candidates="$(awk '/^AutoOffTarget write:/ {split($5,a,"="); value=a[2]} END {print value + 0}' "$forward_log")"
        reverse_candidates="$(awk '/^AutoOffTarget write:/ {split($5,a,"="); value=a[2]} END {print value + 0}' "$reverse_log")"
        forward_padding="$(awk '/^AutoOffTarget write:/ {split($6,a,"="); value=a[2]} END {print value + 0}' "$forward_log")"
        reverse_padding="$(awk '/^AutoOffTarget write:/ {split($6,a,"="); value=a[2]} END {print value + 0}' "$reverse_log")"
        forward_duplicates="$(awk '/^AutoOffTarget write:/ {split($7,a,"="); value=a[2]} END {print value + 0}' "$forward_log")"
        reverse_duplicates="$(awk '/^AutoOffTarget write:/ {split($7,a,"="); value=a[2]} END {print value + 0}' "$reverse_log")"

        candidate_hits=$((forward_candidates + reverse_candidates))
        padding_hits=$((forward_padding + reverse_padding))
        duplicate_hits=$((forward_duplicates + reverse_duplicates))
        raw_hits=$((candidate_hits + padding_hits + duplicate_hits))
        raw_bytes=0
    else
        local adapter_start adapter_end
        adapter_start="$(date +%s%N)"
        adapt_bed_to_candidates \
            "$forward_raw" "$reverse_raw" "$chromosome_map" "$candidate_dir" \
            "$guide_count" "$adapter_stats"
        adapter_end="$(date +%s%N)"
        adapter_seconds="$(elapsed_seconds "$adapter_start" "$adapter_end")"
        IFS=$'\t' read -r raw_hits candidate_hits padding_hits duplicate_hits < "$adapter_stats"
        raw_bytes="$(awk -v f="$(stat -c '%s' "$forward_raw")" \
            -v r="$(stat -c '%s' "$reverse_raw")" 'BEGIN {print f + r}')"
    fi

    local search_without_candidate_write stage1_seconds candidate_bytes
    search_without_candidate_write="$(awk -v total="$elapsed" -v write="$candidate_write_seconds" \
        'BEGIN {value = total - write; if (value < 0) value = 0; printf "%.6f", value}')"
    stage1_seconds="$(awk -v search="$elapsed" -v adapter="$adapter_seconds" \
        'BEGIN {printf "%.6f", search + adapter}')"
    candidate_bytes="$(find "$candidate_dir" -maxdepth 1 -type f -printf '%s\n' |
        awk '{total += $1} END {print total + 0}')"

    echo "${dataset},${figure},${guide_count},${distance},${commit},${STOMATA_OUTPUT_MODE},"\
"${forward_elapsed},${reverse_elapsed},${elapsed},"\
"${candidate_write_seconds},${search_without_candidate_write},"\
"${adapter_seconds},${stage1_seconds},${user_seconds},${system_seconds},"\
"${cpu_percent},${max_rss_kb},${raw_hits},${candidate_hits},${padding_hits},"\
"${duplicate_hits},${raw_bytes},${candidate_bytes},${reported_candidate_dir}" >> "$summary_csv"

    rm -f "$forward_metrics" "$reverse_metrics" "$adapter_stats"
    if [[ "$KEEP_RAW" != "1" ]]; then
        rm -f "$forward_raw" "$reverse_raw"
    fi

    echo "Finished: guides=$guide_count distance=$distance candidates=$candidate_hits mode=$STOMATA_OUTPUT_MODE search=${elapsed}s candidate_write=${candidate_write_seconds}s search_without_write=${search_without_candidate_write}s adapter=${adapter_seconds}s stage1=${stage1_seconds}s"
    sleep "$SLEEP_BETWEEN_RUNS"
}

warm_up_gpu() {
    local genome_index="$1"
    local spacer_file="$2"
    local warmup_output="$3"
    local warmup_log="$4"

    echo "Warming GPU and filesystem caches..."
    "$STOMATA_BIN" \
        --genome "$genome_index" \
        --spacer-file "$spacer_file" \
        --threshold 0 \
        --distance-mode levenshtein \
        --strand plus \
        --format bed \
        --no-deduplicate \
        --no-compute-mismatches \
        --no-scores \
        --no-treat-u-as-t \
        --verbose \
        --output "$warmup_output" \
        > /dev/null 2> "$warmup_log"

    grep -Fq "[GPU: yes]" "$warmup_log" ||
        die "Stomata warm-up did not use the GPU; see $warmup_log"
    rm -f "$warmup_output"
}

# Populates the (already-declared in the caller's scope) forward_dir,
# reverse_dir, output_name variables for a dataset. Shared by the upfront
# validation pass and run_dataset so the two can't drift out of sync.
dataset_dirs() {
    local dataset="$1"

    case "$dataset" in
        hg38)
            forward_dir="$HG38_FORWARD_DIR"
            reverse_dir="$HG38_REVERSE_DIR"
            output_name="stomata_hg38_stage1_outputs_local"
            ;;
        chm13)
            forward_dir="$CHM13_FORWARD_DIR"
            reverse_dir="$CHM13_REVERSE_DIR"
            output_name="stomata_chm13_stage1_outputs_local"
            ;;
        *)
            die "unsupported dataset: $dataset"
            ;;
    esac
}

require_dataset_inputs() {
    local dataset="$1"
    local forward_dir reverse_dir output_name
    dataset_dirs "$dataset"

    [[ -n "$forward_dir" ]] ||
        die "${dataset^^}_FORWARD_DIR must point to the forward chromosome directory"
    [[ -n "$reverse_dir" ]] ||
        die "${dataset^^}_REVERSE_DIR must point to the reverse-complement chromosome directory"
    [[ -d "$forward_dir" ]] || die "forward chromosome directory not found: $forward_dir"
    [[ -d "$reverse_dir" ]] || die "reverse-complement chromosome directory not found: $reverse_dir"

    local output_root="${STOMATA_OUTPUT_PARENT}/${output_name}"
    [[ ! -e "$output_root" ]] || die "output already exists: $output_root"
}

run_dataset() {
    local dataset="$1"
    local forward_dir reverse_dir output_name
    dataset_dirs "$dataset"

    local output_root="${STOMATA_OUTPUT_PARENT}/${output_name}"
    [[ ! -e "$output_root" ]] || die "output already exists: $output_root"

    local prepared_dir="${output_root}/prepared"
    local raw_dir="${output_root}/raw"
    local log_dir="${output_root}/logs"
    local figure1_dir="${output_root}/figure1_runtime_vs_distance_128sgRNAs"
    local figure2_dir="${output_root}/figure2_runtime_vs_numsgRNAs_dist5"
    mkdir -p "$prepared_dir" "$raw_dir" "$log_dir" "$figure1_dir" "$figure2_dir"

    local forward_fasta="${prepared_dir}/${dataset}_forward_padded.fa"
    local reverse_fasta="${prepared_dir}/${dataset}_reverse_complement_padded.fa"
    local forward_index="${forward_fasta}.st"
    local reverse_index="${reverse_fasta}.st"
    local chromosome_map="${prepared_dir}/chromosome_records.tsv"
    local preparation_summary="${prepared_dir}/preparation_summary.csv"

    local prepare_start prepare_end prepare_seconds
    prepare_start="$(date +%s%N)"
    : > "$chromosome_map"
    prepare_strand_genome "$forward_dir" "fw" "$forward_fasta" "$chromosome_map"
    prepare_strand_genome "$reverse_dir" "rc" "$reverse_fasta" "$chromosome_map"
    [[ -s "$chromosome_map" ]] || die "chromosome map is empty: $chromosome_map"
    prepare_end="$(date +%s%N)"
    prepare_seconds="$(elapsed_seconds "$prepare_start" "$prepare_end")"

    local forward_index_metrics="${prepared_dir}/forward_index.time"
    local reverse_index_metrics="${prepared_dir}/reverse_index.time"
    /usr/bin/time -f '%e,%U,%S,%P,%M' -o "$forward_index_metrics" \
        "$STOMATA_BIN" --index-genome "$forward_fasta" \
        > /dev/null 2> "${log_dir}/forward_index.log"
    /usr/bin/time -f '%e,%U,%S,%P,%M' -o "$reverse_index_metrics" \
        "$STOMATA_BIN" --index-genome "$reverse_fasta" \
        > /dev/null 2> "${log_dir}/reverse_index.log"
    [[ -f "$forward_index" ]] || die "Stomata index was not created: $forward_index"
    [[ -f "$reverse_index" ]] || die "Stomata index was not created: $reverse_index"

    local forward_index_elapsed forward_index_user forward_index_system forward_index_cpu forward_index_rss
    local reverse_index_elapsed reverse_index_user reverse_index_system reverse_index_cpu reverse_index_rss
    IFS=',' read -r forward_index_elapsed forward_index_user forward_index_system forward_index_cpu forward_index_rss < "$forward_index_metrics"
    IFS=',' read -r reverse_index_elapsed reverse_index_user reverse_index_system reverse_index_cpu reverse_index_rss < "$reverse_index_metrics"
    {
        echo "dataset,prepare_fastas_seconds,forward_index_elapsed_seconds,reverse_index_elapsed_seconds,forward_index_user_seconds,reverse_index_user_seconds,forward_index_system_seconds,reverse_index_system_seconds,forward_index_cpu_percent,reverse_index_cpu_percent,forward_index_max_resident_kb,reverse_index_max_resident_kb"
        echo "${dataset},${prepare_seconds},${forward_index_elapsed},${reverse_index_elapsed},${forward_index_user},${reverse_index_user},${forward_index_system},${reverse_index_system},${forward_index_cpu%%%},${reverse_index_cpu%%%},${forward_index_rss},${reverse_index_rss}"
    } > "$preparation_summary"
    rm -f "$forward_index_metrics" "$reverse_index_metrics"

    local count
    for count in 1 16 32 64 128; do
        prepare_spacer_file \
            "${GUIDE_FILES[$count]}" "$count" "${prepared_dir}/guides_${count}.tsv"
    done

    write_environment_metadata "${output_root}/environment.txt"
    warm_up_gpu \
        "$forward_index" "${prepared_dir}/guides_128.tsv" \
        "${raw_dir}/warmup_fw.bed" "${log_dir}/warmup_fw.log"
    warm_up_gpu \
        "$reverse_index" "${prepared_dir}/guides_128.tsv" \
        "${raw_dir}/warmup_rc.bed" "${log_dir}/warmup_rc.log"

    local header
    header="dataset,figure,guide_count,edit_distance,stomata_commit,output_mode,forward_search_elapsed_seconds,reverse_search_elapsed_seconds,search_elapsed_seconds,candidate_write_elapsed_seconds,search_without_candidate_write_seconds,adapter_elapsed_seconds,stage1_elapsed_seconds,user_seconds,system_seconds,cpu_percent,max_resident_kb,raw_hits,candidate_hits,padding_hits,duplicate_hits,raw_output_bytes,candidate_output_bytes,candidate_dir"

    local figure1_summary="${figure1_dir}/Stomata_stage1_runtime_summary_fig1.csv"
    local figure2_summary="${figure2_dir}/Stomata_stage1_runtime_summary_fig2.csv"
    echo "$header" > "$figure1_summary"
    echo "$header" > "$figure2_summary"

    local distance candidate_dir run_name
    for distance in 0 1 2 3 4 5 6; do
        run_name="sg128_dist${distance}"
        candidate_dir="${figure1_dir}/output_ed${distance}"
        run_stomata_once \
            "$dataset" "1" "128" "$distance" \
            "${prepared_dir}/guides_128.tsv" "$forward_index" "$reverse_index" \
            "$chromosome_map" "$candidate_dir" "${raw_dir}/${run_name}" \
            "${log_dir}/${run_name}" "$figure1_summary"
    done

    for count in 1 16 32 64 128; do
        run_name="sg${count}_dist5"
        candidate_dir="${figure2_dir}/output_ed5_${count}"
        run_stomata_once \
            "$dataset" "2" "$count" "5" \
            "${prepared_dir}/guides_${count}.tsv" "$forward_index" "$reverse_index" \
            "$chromosome_map" "$candidate_dir" "${raw_dir}/${run_name}" \
            "${log_dir}/${run_name}" "$figure2_summary"
    done

    # The generated FASTAs and indexes are multi-gigabyte reproducible intermediates.
    rm -f "$forward_fasta" "$reverse_fasta" "$forward_index" "$reverse_index"
    rmdir "$raw_dir" 2>/dev/null || true

    echo "Completed $dataset. Results: $output_root"
}

case "$DATASET" in
    hg38|chm13)
        requested_datasets=("$DATASET")
        ;;
    all)
        requested_datasets=("hg38" "chm13")
        ;;
    *)
        die "usage: $0 [hg38|chm13|all]"
        ;;
esac

require_inputs
for requested_dataset in "${requested_datasets[@]}"; do
    require_dataset_inputs "$requested_dataset"
done
build_stomata

for requested_dataset in "${requested_datasets[@]}"; do
    run_dataset "$requested_dataset"
done

echo "All requested Stomata stage-1 benchmarks completed."
