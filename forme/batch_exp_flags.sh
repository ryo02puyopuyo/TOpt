#!/bin/bash

# TODD Experimental フラグ比較実験スクリプト
#
# 使い方:
#   ./batch_exp_flags.sh              # 全フラグ構成を実行
#   ./batch_exp_flags.sh todd         # 元の TODD のみ実行
#   ./batch_exp_flags.sh 1000         # Packed Chi のみ実行
#   ./batch_exp_flags.sh 1000 0100    # Packed Chi と AA Table を実行
#   ./batch_exp_flags.sh todd 1110    # 元の TODD と全ON を実行

# --- 設定 ---
PROJECT_DIR="/home/ryo-mtmt/TOpt"
SCRIPT_DIR="$PROJECT_DIR/forme"

TOPT_BIN="$PROJECT_DIR/bin/TOpt"
INPUT_DIR="$SCRIPT_DIR/tpar_benchmark"
OUTPUT_BASE="$SCRIPT_DIR/results_exp_flags"

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

# ★ 全フラグ構成（引数なしで全て実行）★
ALL_FLAGS=(
    "0000"
    "1000"
    "0100"
    "0010"
    "1110"
)

# --- フラグ構成の決定 ---
if [ $# -gt 0 ]; then
    FLAGS=("$@")
else
    FLAGS=("${ALL_FLAGS[@]}")
fi

# --- アルゴリズム名の変換 ---
get_algo() {
    local f="$1"
    if [ "$f" = "todd" ]; then
        echo "todd"
    else
        echo "todd_exp_${f}"
    fi
}

# --- 出力ディレクトリ作成 ---
mkdir -p "$OUTPUT_BASE"
for f in "${FLAGS[@]}"; do
    mkdir -p "$OUTPUT_BASE/$f"
done

# --- サマリーファイル ---
SUMMARY="$OUTPUT_BASE/summary.csv"
# ヘッダーが無い場合のみ作成
if [ ! -f "$SUMMARY" ]; then
    echo "circuit,flags,input_tcount,output_tcount,time_s" > "$SUMMARY"
fi

echo "========================================"
echo " TODD Flag Comparison"
echo " Circuits: ${#CIRCUITS[@]}"
echo " Configs:  ${FLAGS[*]}"
echo " Total:    $(( ${#CIRCUITS[@]} * ${#FLAGS[@]} )) runs"
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

    for flags in "${FLAGS[@]}"; do
        algo=$(get_algo "$flags")
        outfile="$OUTPUT_BASE/$flags/${name}.out"

        # 既に結果がある場合はスキップ（再実行したい場合は .out を削除）
        if [ -f "$outfile" ]; then
            echo "  $flags ... SKIP (already exists)"
            continue
        fi

        echo -n "  $flags ($algo) ... "
        "$TOPT_BIN" circuit "$path" -a "$algo" > "$outfile" 2>&1
        
        if [ $? -eq 0 ]; then
            in_t=$(grep -m1 "^T count = " "$outfile" | awk '{print $NF}')
            out_t=$(grep "Output T Count (PhasePolynomial)" "$outfile" | awk '{print $NF}')

            # Execution Time (experimental) or Execution time (original)
            t_ms=$(grep "Execution Time" "$outfile" | grep -oP '[0-9.]+' | head -1)
            if [ -n "$t_ms" ]; then
                t_s=$(echo "scale=3; $t_ms / 1000" | bc 2>/dev/null || echo "$t_ms")
            else
                t_s=$(grep "Execution time:" "$outfile" | grep -oP '[0-9.]+' | head -1)
            fi

            echo "T: $in_t -> $out_t  (${t_s}s)"
            echo "$name,$flags,$in_t,$out_t,$t_s" >> "$SUMMARY"
        else
            echo "ERROR"
            echo "$name,$flags,ERR,ERR,ERR" >> "$SUMMARY"
        fi
    done
done

echo ""
echo "========================================"
echo "Done! Results: $OUTPUT_BASE"
echo "CSV: $SUMMARY"
echo "========================================"
