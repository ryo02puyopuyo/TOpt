#!/bin/bash

# Multi-circuit AA sketch-before-basis experiment script.
#
# Usage:
#   cd /home/ryo-mtmt/TOpt
#   bash forme/batch_multi_aa_sketch.sh
#
# Slurm:
#   sbatch --wrap="bash forme/batch_multi_aa_sketch.sh"

PROJECT_DIR="/home/ryo-mtmt/TOpt"
SCRIPT_DIR="$PROJECT_DIR/forme"
INPUT_DIR="$SCRIPT_DIR/tpar_benchmark"

TOPT_BIN="$PROJECT_DIR/bin/TOpt"
OUTPUT_BASE="$SCRIPT_DIR/sketch_experiments/aa_sketch_multi_$(date +%Y%m%d_%H%M%S)"

# Default set: tune this list before large Slurm runs.
CIRCUITS=(
    "barenco_tof_3.tfc"
    "barenco_tof_4.tfc"
    "barenco_tof_5.tfc"
    "barenco_tof_10.tfc"
    "gf2^4_mult.tfc"
    "gf2^5_mult.tfc"
    "gf2^6_mult.tfc"
    "gf2^7_mult.tfc"
    "gf2^8_mult.tfc"
    "grover_5.tfc"
    "hwb6.tfc"
    "mod5_4.tfc"
    "mod_mult_55.tfc"
    "mod_red_21.tfc"
    "qcla_com_7.tfc"
    "qft_4.tfc"
    "rc_adder_6.tfc"
    "tof_3.tfc"
    "tof_4.tfc"
    "tof_5.tfc"
    "tof_10.tfc"
    "vbe_adder_3.tfc"
    #"adder_8.tfc"
    "gf2^9_mult.tfc"
    #"ham15-low.tfc"
)

ALGOS=(
    "lx2_dynamic_repair_chi_packed_aa"
    "lx2_dynamic_repair_chi_packed_aa_sketch_5"
    "lx2_dynamic_repair_chi_packed_aa_sketch_8"
    "lx2_dynamic_repair_chi_packed_aa_sketch_10"
    "lx2_dynamic_repair_chi_packed_aa_sketch_13"
    "lx2_dynamic_repair_chi_packed_aa_sketch_15"
)

safe_name() {
    echo "$1" | tr '^' '_' | tr '/' '_' | tr ':' '_' | tr ' ' '_'
}

extract_value() {
    local pattern="$1"
    local file="$2"
    grep -m1 "$pattern" "$file" | awk -F ':' '{gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}'
}

extract_number() {
    local pattern="$1"
    local file="$2"
    grep -m1 "$pattern" "$file" | grep -oP '[0-9.]+' | head -1
}

extract_last_number() {
    local pattern="$1"
    local file="$2"
    grep "$pattern" "$file" | grep -oP '[0-9.]+' | tail -1
}

mkdir -p "$OUTPUT_BASE"

SUMMARY="$OUTPUT_BASE/summary.csv"
echo "circuit,algorithm,input_tcount,output_tcount,total_exec_ms,overall_exec_s,pairs_tested,sketch_filtered,filtered_miss,nullspace_runs,rebuild_count,local_add_count,aa_ms,sketch_ms,chi_update_ms,basis_work_ms,nullspace_ms,sketch_rej_rate_pct,avg_sketch_rows,avg_sketch_basis,avg_affected_rows,avg_active_rows,fail_count,outfile" > "$SUMMARY"

if [ ! -x "$TOPT_BIN" ]; then
    echo "ERROR: $TOPT_BIN not found or not executable"
    exit 1
fi

echo "========================================"
echo " multi-circuit AA sketch-before-basis experiments"
echo " Input dir: $INPUT_DIR"
echo " Output:    $OUTPUT_BASE"
echo " Circuits:  ${#CIRCUITS[@]}"
echo " Algos:     ${#ALGOS[@]}"
echo " Runs:      $(( ${#CIRCUITS[@]} * ${#ALGOS[@]} ))"
echo "========================================"

for circuit in "${CIRCUITS[@]}"; do
    input="$INPUT_DIR/$circuit"
    circuit_name="${circuit%.tfc}"
    circuit_safe="$(safe_name "$circuit_name")"

    if [ ! -f "$input" ]; then
        echo "[SKIP] $circuit not found"
        continue
    fi

    echo ""
    echo "=== $circuit ==="

    for algo in "${ALGOS[@]}"; do
        algo_safe="$(safe_name "$algo")"
        outfile="$OUTPUT_BASE/${circuit_safe}_${algo_safe}.out"

        if [ -f "$outfile" ]; then
            echo "$algo ... SKIP (already exists)"
            continue
        fi

        echo -n "$algo ... "
        "$TOPT_BIN" circuit "$input" -a "$algo" > "$outfile" 2>&1
        status=$?

        if [ $status -ne 0 ]; then
            echo "ERROR"
            echo "$circuit_name,$algo,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,$outfile" >> "$SUMMARY"
            continue
        fi

        in_t=$(grep -m1 "^T count = " "$outfile" | awk '{print $NF}')
        out_t=$(grep "Output T Count (PhasePolynomial)" "$outfile" | awk '{print $NF}' | tail -1)
        exec_ms=$(extract_number "Execution Time" "$outfile")
        exec_s=$(extract_number "Execution time:" "$outfile")
        pairs_tested=$(extract_value "Pairs tested" "$outfile")
        sketch_filtered=$(extract_value "Sketch filtered" "$outfile")
        filtered_miss=$(extract_value "Filtered miss" "$outfile")
        nullspace_runs=$(extract_value "Nullspace runs" "$outfile")
        rebuild_count=$(extract_value "Rebuild count" "$outfile")
        local_add_count=$(extract_value "Local add count" "$outfile")
        aa_ms=$(extract_number "AA Table time" "$outfile")
        sketch_ms=$(extract_number "Sketch time" "$outfile")
        chi_update_ms=$(extract_number "Chi update time" "$outfile")
        basis_work_ms=$(extract_number "Basis work time" "$outfile")
        nullspace_ms=$(extract_number "Nullspace time" "$outfile")
        sketch_rej_rate_pct=$(extract_value "Sketch rej rate" "$outfile" | tr -d '%')
        avg_sketch_rows=$(extract_value "Avg sketch rows" "$outfile")
        avg_sketch_basis=$(extract_value "Avg sketch basis" "$outfile")
        avg_affected_rows=$(extract_value "Avg affected" "$outfile" | awk '{print $1}')
        avg_active_rows=$(extract_value "Avg active rows" "$outfile")
        fail_count=$(extract_last_number "Fail count" "$outfile")

        echo "T: $in_t -> $out_t, time: ${exec_s}s, sketch_filtered: ${sketch_filtered:-}, sketch_reject: ${sketch_rej_rate_pct:-}%, basis: ${basis_work_ms:-}ms"
        echo "$circuit_name,$algo,$in_t,$out_t,$exec_ms,$exec_s,$pairs_tested,$sketch_filtered,$filtered_miss,$nullspace_runs,$rebuild_count,$local_add_count,$aa_ms,$sketch_ms,$chi_update_ms,$basis_work_ms,$nullspace_ms,$sketch_rej_rate_pct,$avg_sketch_rows,$avg_sketch_basis,$avg_affected_rows,$avg_active_rows,$fail_count,$outfile" >> "$SUMMARY"
    done
done

echo ""
echo "========================================"
echo "Done"
echo "Results: $OUTPUT_BASE"
echo "CSV:     $SUMMARY"
echo "========================================"
