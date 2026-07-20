#!/bin/bash

# Multi-circuit sketch experiment script.
#
# Usage:
#   cd /home/ryo-mtmt/TOpt
#   bash forme/batch_multi_sketch.sh
#
# Slurm:
#   sbatch --wrap="bash forme/batch_multi_sketch.sh"

PROJECT_DIR="/home/ryo-mtmt/TOpt"
SCRIPT_DIR="$PROJECT_DIR/forme"
INPUT_DIR="$SCRIPT_DIR/tpar_benchmark"

TOPT_BIN="$PROJECT_DIR/bin/TOpt"
OUTPUT_BASE="$SCRIPT_DIR/sketch_experiments/multi_$(date +%Y%m%d_%H%M%S)"

# Default set: keep this moderate for Slurm experiments.
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
    "adder_8.tfc"
    "gf2^9_mult.tfc"
    "ham15-low.tfc"
)

ALGOS=(
    "todd_exp_packedlocal_diff_1110"
    "todd_exp_packedlocal_diff_1111_p5"
    "todd_exp_packedlocal_diff_1111_p10"
    "todd_exp_packedlocal_diff_1111_p15"
    "todd_exp_packedlocal_diff_1111_p20"
)

# Optional wider sweep. Append these to ALGOS if needed.
MORE_P_ALGOS=(
    "todd_exp_packedlocal_diff_1111_p5"
    "todd_exp_packedlocal_diff_1111_p20"
    "todd_exp_packedlocal_diff_1111_p25"
    "todd_exp_packedlocal_diff_1111_p30"
)

# Optional heavier circuits. Append these to CIRCUITS if needed.
HEAVY_CIRCUITS=(
    "adder_8.tfc"
    "gf2^9_mult.tfc"
    "gf2^10_mult.tfc"
    "ham15-low.tfc"
    "ham15-med.tfc"
    "hwb8.tfc"
    "qcla_adder_10.tfc"
    "qcla_mod_7.tfc"
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
    grep -m1 "$pattern" "$file" | grep -oP '[0-9.]+' | tail -1
}

percent_rate() {
    local numerator="$1"
    local denominator="$2"
    awk -v n="${numerator:-0}" -v d="${denominator:-0}" 'BEGIN {
        if (d + 0 <= 0) print "";
        else printf "%.2f", 100.0 * (n + 0) / (d + 0);
    }'
}

mkdir -p "$OUTPUT_BASE"

SUMMARY="$OUTPUT_BASE/summary.csv"
echo "circuit,algorithm,input_tcount,output_tcount,total_exec_ms,overall_exec_s,chi_ms,aa_ms,nullspace_ms,sketch_ms,reject_rate_pct,sketch_mode,sketch_attempts,sketch_skipped,sketch_rejected,sketch_rank_rejected,sketch_pair_rejected,pair_reject_rate_pct,sketch_passed,sketch_rows_avg,sketch_y_vectors,sketch_pair_y,sketch_y_full_ok,sketch_y_full_ng,sketch_y_reduced,full_ns_skipped,full_ns_skip_rate_pct,full_ns_computed,fail_count,outfile" > "$SUMMARY"

if [ ! -x "$TOPT_BIN" ]; then
    echo "ERROR: $TOPT_BIN not found or not executable"
    exit 1
fi

echo "========================================"
echo " multi-circuit sketch experiments"
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
            echo "$circuit_name,$algo,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,ERR,$outfile" >> "$SUMMARY"
            continue
        fi

        in_t=$(grep -m1 "^T count = " "$outfile" | awk '{print $NF}')
        out_t=$(grep "Output T Count (PhasePolynomial)" "$outfile" | awk '{print $NF}' | tail -1)
        exec_ms=$(extract_number "Execution Time" "$outfile")
        exec_s=$(extract_number "Execution time:" "$outfile")
        chi_ms=$(extract_number "Chi calculation" "$outfile")
        aa_ms=$(extract_number "AA Table calc" "$outfile")
        ns_ms=$(extract_number "Nullspace calc" "$outfile")
        sketch_mode=$(extract_value "Sketch mode" "$outfile")
        sketch_attempts=$(extract_value "Sketch attempts" "$outfile")
        sketch_skipped=$(extract_value "Sketch skipped" "$outfile")
        sketch_rejected=$(extract_value "Sketch rejected" "$outfile")
        sketch_rank_rejected=$(extract_value "Sketch rank rej" "$outfile")
        sketch_pair_rejected=$(extract_value "Sketch pair rej" "$outfile")
        sketch_passed=$(extract_value "Sketch passed" "$outfile")
        sketch_rows_avg=$(extract_value "Sketch rows avg" "$outfile")
        sketch_y_vectors=$(extract_value "Sketch y vectors" "$outfile")
        sketch_pair_y=$(extract_value "Sketch pair y" "$outfile")
        sketch_y_full_ok=$(extract_value "Sketch y full ok" "$outfile")
        sketch_y_full_ng=$(extract_value "Sketch y full ng" "$outfile")
        sketch_y_reduced=$(extract_value "Sketch y reduced" "$outfile")
        full_ns_skipped=$(extract_value "Full NS skipped" "$outfile")
        full_ns_computed=$(extract_value "Full NS computed" "$outfile")
        sketch_ms=$(extract_number "Sketch time" "$outfile")
        fail_count=$(extract_last_number "Fail count" "$outfile")
        reject_rate_pct=$(percent_rate "$sketch_rejected" "$sketch_attempts")
        pair_reject_rate_pct=$(percent_rate "$sketch_pair_rejected" "$sketch_attempts")
        full_ns_skip_rate_pct=$(percent_rate "$full_ns_skipped" "$sketch_attempts")

        echo "T: $in_t -> $out_t, time: ${exec_s}s, nullspace: ${ns_ms:-}ms, sketch: ${sketch_ms:-}ms, reject: ${reject_rate_pct:-}%"
        echo "$circuit_name,$algo,$in_t,$out_t,$exec_ms,$exec_s,$chi_ms,$aa_ms,$ns_ms,$sketch_ms,$reject_rate_pct,$sketch_mode,$sketch_attempts,$sketch_skipped,$sketch_rejected,$sketch_rank_rejected,$sketch_pair_rejected,$pair_reject_rate_pct,$sketch_passed,$sketch_rows_avg,$sketch_y_vectors,$sketch_pair_y,$sketch_y_full_ok,$sketch_y_full_ng,$sketch_y_reduced,$full_ns_skipped,$full_ns_skip_rate_pct,$full_ns_computed,$fail_count,$outfile" >> "$SUMMARY"
    done
done

echo ""
echo "========================================"
echo "Done"
echo "Results: $OUTPUT_BASE"
echo "CSV:     $SUMMARY"
echo "========================================"
