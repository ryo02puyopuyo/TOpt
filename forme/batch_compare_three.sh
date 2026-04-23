#!/bin/bash

# 3手法比較実験スクリプト
#
# 使い方:
#   ./batch_compare_three.sh
#
# スクリプト内で指定した回路に対して，
#   - todd_exp_packedlocal_diff_1110
#   - todd_exp_1110
#   - lx2_dynamic_repair_chi_packed_aa
# の3手法を実行する．

# --- 設定 ---
PROJECT_DIR="/home/rest/forpyzx/TOpt"
SCRIPT_DIR="$PROJECT_DIR/forme"

TOPT_BIN="$PROJECT_DIR/bin/TOpt"
INPUT_DIR="$SCRIPT_DIR/tpar_benchmark"
OUTPUT_BASE="$SCRIPT_DIR/results_compare_three"

# ★ 実行する回路をここに記述 ★
CIRCUITS=(
    "adder_8.tfc"
    "barenco_tof_3.tfc"
    "barenco_tof_4.tfc"
    "barenco_tof_5.tfc"
    "barenco_tof_10.tfc"
    "csla_mux_3.tfc"
    "csum_mux_9.tfc"
    "gf2^4_mult.tfc"
    "gf2^5_mult.tfc"
    "gf2^6_mult.tfc"
    "gf2^7_mult.tfc"
    "gf2^8_mult.tfc"
    "gf2^9_mult.tfc"
    "grover_5.tfc"
    "ham15-high.tfc"
    "ham15-med.tfc"
    "ham15-low.tfc"
    "hwb6.tfc"
    "hwb8.tfc"
    "mod5_4.tfc"
    "mod_adder_1024.tfc"
    "mod_mult_55.tfc"
    "mod_red_21.tfc"
    "qcla_adder_10.tfc"
    "qcla_com_7.tfc"
    "qcla_mod_7.tfc"
    "qft_4.tfc"
    "rc_adder_6.tfc"
    "tof_3.tfc"
    "tof_4.tfc"
    "tof_5.tfc"
    "tof_10.tfc"
    "vbe_adder_3.tfc"
)

# ★ 比較するアルゴリズム ★
ALGOS=(
    "todd_exp_packedlocal_diff_1110"
    "todd_exp_1110"
    "lx2_dynamic_repair_chi_packed_aa"
)

algo_to_dir() {
    local algo="$1"
    case "$algo" in
        "todd_exp_packedlocal_diff_1110")
            echo "todd_exp_diff"
            ;;
        "todd_exp_1110")
            echo "todd_exp_1110"
            ;;
        "lx2_dynamic_repair_chi_packed_aa")
            echo "lx2_repair_chi_packed_aa"
            ;;
        *)
            echo "$algo"
            ;;
    esac
}

# --- 出力ディレクトリ作成 ---
mkdir -p "$OUTPUT_BASE"
for algo in "${ALGOS[@]}"; do
    mkdir -p "$OUTPUT_BASE/$(algo_to_dir "$algo")"
done

# --- サマリーファイル ---
SUMMARY="$OUTPUT_BASE/summary.csv"
if [ ! -f "$SUMMARY" ]; then
    echo "circuit,algorithm,input_tcount,output_tcount,time_s" > "$SUMMARY"
fi

echo "========================================"
echo " Three-Algorithm Comparison"
echo " Circuits: ${#CIRCUITS[@]}"
echo " Algos:    ${ALGOS[*]}"
echo " Total:    $(( ${#CIRCUITS[@]} * ${#ALGOS[@]} )) runs"
echo "========================================"

for circuit in "${CIRCUITS[@]}"; do
    name="${circuit%.tfc}"
    path="$INPUT_DIR/$circuit"

    if [ ! -f "$path" ]; then
        echo "[SKIP] $circuit not found"
        continue
    fi

    echo ""
    echo "=== $circuit ==="

    for algo in "${ALGOS[@]}"; do
        outdir="$OUTPUT_BASE/$(algo_to_dir "$algo")"
        outfile="$outdir/${name}.out"

        if [ -f "$outfile" ]; then
            echo "  $algo ... SKIP (already exists)"
            continue
        fi

        echo -n "  $algo ... "
        "$TOPT_BIN" circuit "$path" -a "$algo" > "$outfile" 2>&1

        if [ $? -eq 0 ]; then
            in_t=$(grep -m1 "^T count = " "$outfile" | awk '{print $NF}')
            out_t=$(grep "Output T Count (PhasePolynomial)" "$outfile" | awk '{print $NF}' | tail -1)

            t_ms=$(grep "Execution Time" "$outfile" | grep -oP '[0-9.]+' | head -1)
            if [ -n "$t_ms" ]; then
                t_s=$(echo "scale=3; $t_ms / 1000" | bc 2>/dev/null || echo "$t_ms")
            else
                t_s=$(grep "Execution time:" "$outfile" | grep -oP '[0-9.]+' | head -1)
            fi

            echo "T: $in_t -> $out_t  (${t_s}s)"
            echo "$name,$algo,$in_t,$out_t,$t_s" >> "$SUMMARY"
        else
            echo "ERROR"
            echo "$name,$algo,ERR,ERR,ERR" >> "$SUMMARY"
        fi
    done
done

echo ""
echo "========================================"
echo "Done! Results: $OUTPUT_BASE"
echo "CSV: $SUMMARY"
echo "========================================"
