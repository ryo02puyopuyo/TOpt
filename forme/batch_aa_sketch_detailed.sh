#!/bin/bash
set -u

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

PROJECT_DIR="/home/ryo-mtmt/TOpt"
SCRIPT_DIR="$PROJECT_DIR/forme"
INPUT_DIR="$SCRIPT_DIR/tpar_benchmark"
TOPT_BIN="$PROJECT_DIR/bin/TOpt"
OUTPUT_BASE="$SCRIPT_DIR/sketch_experiments/aa_sketch_detail_$TIMESTAMP"
OUT_DIR="$OUTPUT_BASE/logs"
CSV_FILE="$OUTPUT_BASE/summary.csv"

# Edit these lists before large runs.
CIRCUITS=(
    #"barenco_tof_3.tfc"
    #"barenco_tof_4.tfc"
    #"barenco_tof_5.tfc"
    #"barenco_tof_10.tfc"
    #"gf2^4_mult.tfc"
    #"gf2^5_mult.tfc"
    #"gf2^6_mult.tfc"
    #"gf2^7_mult.tfc"
    #"gf2^8_mult.tfc"
    #"grover_5.tfc"
    #"hwb6.tfc"
    #"mod5_4.tfc"
    #"mod_mult_55.tfc"
    #"mod_red_21.tfc"
    #"qcla_com_7.tfc"
    #"qft_4.tfc"
    #"rc_adder_6.tfc"
    #"tof_3.tfc"
    #"tof_4.tfc"
    #"tof_5.tfc"
    #"tof_10.tfc"
    #"vbe_adder_3.tfc"
    #"adder_8.tfc"
    #"gf2^9_mult.tfc"
    #"ham15-low.tfc"
    #"ham15-med.tfc"
    "gf2^10_mult.tfc"
    "gf2^16_mult.tfc"
    "qcla_mod_7.tfc"
)

P_LIST=(
    3
    5
    8
    10
    13
    15
    20
    25
    30
)

# Available variants:
#   reuse      : lx2_dynamic_repair_chi_packed_aa_sketch_p
#   noreuse    : lx2_dynamic_repair_chi_packed_aa_sketch_noreuse_p
#   basisfirst : lx2_dynamic_repair_chi_packed_aa_basis_first_sketch_p
VARIANTS=(
    "reuse"
    "noreuse"
)

RUN_BASELINE=1

mkdir -p "$OUT_DIR"

extract_value() {
    local key="$1"
    local file="$2"
    grep -m 1 "$key" "$file" | awk -F': ' '{print $2}' | awk '{print $1}'
}

extract_number() {
    local key="$1"
    local file="$2"
    grep -m 1 "$key" "$file" | grep -oE '[0-9]+(\.[0-9]+)?' | head -n 1
}

safe_name() {
    echo "$1" | sed 's/[^A-Za-z0-9_.-]/_/g'
}

percent_rate() {
    local num="${1:-0}"
    local den="${2:-0}"
    awk -v n="$num" -v d="$den" 'BEGIN { if (d == 0) print ""; else printf "%.3f", 100.0 * n / d }'
}

sec_from_ms() {
    local ms="${1:-}"
    awk -v ms="$ms" 'BEGIN { if (ms == "") print ""; else printf "%.3f", ms / 1000.0 }'
}

algo_for_variant() {
    local variant="$1"
    local p="$2"
    case "$variant" in
        reuse)
            echo "lx2_dynamic_repair_chi_packed_aa_sketch_$p"
            ;;
        noreuse)
            echo "lx2_dynamic_repair_chi_packed_aa_sketch_noreuse_$p"
            ;;
        basisfirst)
            echo "lx2_dynamic_repair_chi_packed_aa_basis_first_sketch_$p"
            ;;
        *)
            echo ""
            ;;
    esac
}

write_header() {
    cat > "$CSV_FILE" <<'EOF'
circuit,variant,p,algorithm,input_tcount,output_tcount,total_exec_ms,total_exec_s,overall_exec_s,pairs_tested,sketch_filtered,sketch_skip_reuse,full_after_sketch,full_after_sketch_rate_pct,filtered_miss,nullspace_runs,rebuild_count,local_add_count,basis_update_count,basis_reuse_rate_pct,aa_ms,sketch_ms,sketch_build_ms,sketch_basis_ms,sketch_mem_ms,sketch_other_ms,chi_update_ms,basis_work_ms,nullspace_ms,sketch_rej_rate_pct,avg_sketch_rows,avg_sketch_basis,avg_affected_rows,avg_active_rows,reuse_mem_reject,reuse_ns_runs,reuse_no_pair_y,reuse_no_reduce,reuse_reduced,reuse_mem_rej_rate_pct,reuse_no_red_rate_pct,fail_count,outfile
EOF
}

append_summary() {
    local circuit="$1"
    local variant="$2"
    local p="$3"
    local algo="$4"
    local outfile="$5"

    local input_t output_t total_exec_ms total_exec_s overall_exec_s
    local pairs sketch_filtered sketch_skip_reuse full_after_sketch full_after_sketch_rate
    local filtered_miss nullspace_runs rebuild_count local_add_count basis_update_count basis_reuse_rate
    local aa_ms sketch_ms sketch_build_ms sketch_basis_ms sketch_mem_ms sketch_other_ms
    local chi_update_ms basis_work_ms nullspace_ms sketch_rej_rate
    local avg_sketch_rows avg_sketch_basis avg_affected avg_active
    local reuse_mem_reject reuse_ns_runs reuse_no_pair_y reuse_no_reduce reuse_reduced reuse_mem_rej_rate reuse_no_red_rate
    local fail_count

    input_t=$(extract_number "Initial T-count" "$outfile")
    output_t=$(extract_number "Final T-count" "$outfile")
    total_exec_ms=$(extract_value "Execution Time" "$outfile")
    total_exec_s=$(sec_from_ms "$total_exec_ms")
    overall_exec_s=$(extract_number "Execution time:" "$outfile")
    pairs=$(extract_value "Pairs tested" "$outfile")
    sketch_filtered=$(extract_value "Sketch filtered" "$outfile")
    sketch_skip_reuse=$(extract_value "Sketch skip reuse" "$outfile")
    filtered_miss=$(extract_value "Filtered miss" "$outfile")
    nullspace_runs=$(extract_value "Nullspace runs" "$outfile")
    rebuild_count=$(extract_value "Rebuild count" "$outfile")
    local_add_count=$(extract_value "Local add count" "$outfile")
    aa_ms=$(extract_value "AA Table time" "$outfile")
    sketch_ms=$(extract_value "Sketch time" "$outfile")
    sketch_build_ms=$(extract_value "Sketch build time" "$outfile")
    sketch_basis_ms=$(extract_value "Sketch basis time" "$outfile")
    sketch_mem_ms=$(extract_value "Sketch mem time" "$outfile")
    sketch_other_ms=$(extract_value "Sketch other time" "$outfile")
    chi_update_ms=$(extract_value "Chi update time" "$outfile")
    basis_work_ms=$(extract_value "Basis work time" "$outfile")
    nullspace_ms=$(extract_value "Nullspace time" "$outfile")
    sketch_rej_rate=$(extract_value "Sketch rej rate" "$outfile" | tr -d '%')
    avg_sketch_rows=$(extract_value "Avg sketch rows" "$outfile")
    avg_sketch_basis=$(extract_value "Avg sketch basis" "$outfile")
    avg_affected=$(extract_value "Avg affected" "$outfile")
    avg_active=$(extract_value "Avg active rows" "$outfile")
    reuse_mem_reject=$(extract_value "Reuse mem reject" "$outfile")
    reuse_ns_runs=$(extract_value "Reuse NS runs" "$outfile")
    reuse_no_pair_y=$(extract_value "Reuse no pair y" "$outfile")
    reuse_no_reduce=$(extract_value "Reuse no reduce" "$outfile")
    reuse_reduced=$(extract_value "Reuse reduced" "$outfile")
    reuse_mem_rej_rate=$(extract_value "Reuse mem rej rate" "$outfile" | tr -d '%')
    reuse_no_red_rate=$(extract_value "Reuse no-red rate" "$outfile" | tr -d '%')
    fail_count=$(extract_number "Fail count" "$outfile")

    basis_update_count=$(awk -v r="${rebuild_count:-0}" -v l="${local_add_count:-0}" 'BEGIN { printf "%d", r + l }')
    basis_reuse_rate=$(percent_rate "${local_add_count:-0}" "$basis_update_count")
    full_after_sketch=$(awk -v pairs="${pairs:-0}" -v filt="${sketch_filtered:-0}" 'BEGIN { printf "%d", pairs - filt }')
    full_after_sketch_rate=$(percent_rate "$full_after_sketch" "${pairs:-0}")

    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$circuit" "$variant" "$p" "$algo" \
        "$input_t" "$output_t" "$total_exec_ms" "$total_exec_s" "$overall_exec_s" \
        "$pairs" "$sketch_filtered" "$sketch_skip_reuse" "$full_after_sketch" "$full_after_sketch_rate" \
        "$filtered_miss" "$nullspace_runs" "$rebuild_count" "$local_add_count" "$basis_update_count" "$basis_reuse_rate" \
        "$aa_ms" "$sketch_ms" "$sketch_build_ms" "$sketch_basis_ms" "$sketch_mem_ms" "$sketch_other_ms" \
        "$chi_update_ms" "$basis_work_ms" "$nullspace_ms" "$sketch_rej_rate" \
        "$avg_sketch_rows" "$avg_sketch_basis" "$avg_affected" "$avg_active" \
        "$reuse_mem_reject" "$reuse_ns_runs" "$reuse_no_pair_y" "$reuse_no_reduce" "$reuse_reduced" \
        "$reuse_mem_rej_rate" "$reuse_no_red_rate" "$fail_count" "$outfile" >> "$CSV_FILE"
}

run_one() {
    local circuit="$1"
    local variant="$2"
    local p="$3"
    local algo="$4"
    local circuit_base safe_algo out_tfc outfile

    circuit_base="$(basename "$circuit" .tfc)"
    safe_algo="$(safe_name "$algo")"
    out_tfc="$OUT_DIR/${circuit_base}_${safe_algo}.tfc"
    outfile="$OUT_DIR/${circuit_base}_${safe_algo}.out"

    echo "[$(date '+%F %T')] RUN circuit=$circuit algo=$algo"
    "$TOPT_BIN" circuit "$INPUT_DIR/$circuit" -a "$algo" -o "$out_tfc" > "$outfile" 2>&1
    append_summary "$circuit" "$variant" "$p" "$algo" "$outfile"
}

write_header

echo "Output directory: $OUTPUT_BASE"
echo "CSV: $CSV_FILE"
echo "Circuits: ${CIRCUITS[*]}"
echo "p list: ${P_LIST[*]}"
echo "Variants: ${VARIANTS[*]}"

for circuit in "${CIRCUITS[@]}"; do
    if [ ! -f "$INPUT_DIR/$circuit" ]; then
        echo "Skip missing circuit: $INPUT_DIR/$circuit"
        continue
    fi

    if [ "$RUN_BASELINE" = "1" ]; then
        run_one "$circuit" "baseline" "" "lx2_dynamic_repair_chi_packed_aa"
    fi

    for p in "${P_LIST[@]}"; do
        for variant in "${VARIANTS[@]}"; do
            algo="$(algo_for_variant "$variant" "$p")"
            if [ -z "$algo" ]; then
                echo "Skip unknown variant: $variant"
                continue
            fi
            run_one "$circuit" "$variant" "$p" "$algo"
        done
    done
done

echo "Done: $CSV_FILE"
