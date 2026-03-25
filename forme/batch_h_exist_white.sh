#!/bin/bash

# --- 設定項目 ---
INPUT_DIR="/home/ryo-mtmt/TOpt/forme/tpar_benchmark"      # .tfcファイルが入っているディレクトリ
OUTPUT_DIR="./results_h_既存手法"         # 結果を保存するディレクトリ

# ★ 特定のファイルだけ実行したい場合は、ここにスペース区切りで入力してください
# （例: "adder_8.tfc" "hwb8.tfc"）
# ここが空（()）の場合は、下の除外リストが適用されます。
INCLUDE_LIST=(
    "gf2^8_mult.tfc"
    "gf2^9_mult.tfc"
    "ham15-med.tfc"
    "adder_8.tfc"
    "qcla_mod_7.tfc"
    "qcla_adder_10.tfc"
    "mod_adder_1024.tfc"
    "mod_red_21.tfc"
#    "ham15-med.tfc"
 #   "ham15-low.tfc"
  #  "hwb6.tfc"
   # "mod5_4.tfc"
   # "mod_red_21.tfc"
   # "mod_mult_55.tfc"
   # "qcla_com_7.tfc"
   # "qcla_mod_7.tfc"
   # "qft_4.tfc"
   # "rc_adder_6.tfc"
   # "tof_3.tfc"
   # "tof_4.tfc"
   # "tof_5.tfc"
   # "vbe_adder_3.tfc"
)

# ★ 実行したくないファイル名（INCLUDE_LISTが空の時のみ有効）
EXCLUDE_LIST=(
    "adder_8.tfc"
    "gf2^8_mult.tfc"
    "gf2^9_mult.tfc"
    "gf2^10_mult.tfc"
    "gf2^16_mult.tfc"
    "gf2^32_mult.tfc"
    "gf2^64_mult.tfc"
    "gf2^128_mult.tfc"
    "gf2^256_mult.tfc"
    "cycle_17_3.tfc"
    "hwb8.tfc"
    "hwb10.tfc"
    "hwb11.tfc"
    "hwb12.tfc"
    "mod_adder_1048576.tfc"
    "mod_adder_1024.tfc"
    "ham15-high.tfc"
    "ham15-med.tfc"
    "qcla_adder_10.tfc"
    "qcla_mod_7.tfc"
)

# ライブラリパスの設定（前回のlibm4ri対策）
#export LD_LIBRARY_PATH=/home/ryo-mtmt/.local/lib:$LD_LIBRARY_PATH

# ----------------

# 出力ディレクトリの作成
mkdir -p "$OUTPUT_DIR"

# モード判定
if [ ${#INCLUDE_LIST[@]} -gt 0 ]; then
    RUN_MODE="INCLUDE (Only specified files)"
    TARGET_LIST=("${INCLUDE_LIST[@]}")
else
    RUN_MODE="EXCLUDE (All except excluded files)"
fi

echo "Starting optimization..."
echo "Mode:             $RUN_MODE"
echo "Input Directory:  $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "---------------------------------------"

count=0
skip_count=0

# ファイル処理のメインループ
for tfc_path in "$INPUT_DIR"/*.tfc; do
    [ -e "$tfc_path" ] || continue
    filename=$(basename "$tfc_path")

    # 実行判定フラグ
    should_run=false

    if [ ${#INCLUDE_LIST[@]} -gt 0 ]; then
        # ホワイトリストモード：INCLUDE_LISTに含まれているか？
        for target in "${INCLUDE_LIST[@]}"; do
            if [[ "$filename" == "$target" ]]; then
                should_run=true
                break
            fi
        done
    else
        # ブラックリストモード：EXCLUDE_LISTに含まれていないか？
        is_excluded=false
        for excluded in "${EXCLUDE_LIST[@]}"; do
            if [[ "$filename" == "$excluded" ]]; then
                is_excluded=true
                break
            fi
        done
        if ! $is_excluded; then
            should_run=true
        fi
    fi

    # 実行
    if $should_run; then
        output_file="$OUTPUT_DIR/${filename%.tfc}.out"
        echo "Running: $filename -> $(basename "$output_file")"
        
        /home/ryo-mtmt/TOpt/bin/TOpt circuit "$tfc_path" -a todd > "$output_file"

        if [ $? -eq 0 ]; then
            ((count++))
        else
            echo "Error occurred during: $filename"
        fi
    else
        ((skip_count++))
    fi
done

echo "---------------------------------------"
echo "Done!"
echo "Optimized: $count files"
echo "Skipped:   $skip_count files"