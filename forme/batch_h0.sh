#!/bin/bash

#回路分割手法

# --- 設定項目 ---
INPUT_DIR="/home/ryo-mtmt/TOpt/forme/tpar_benchmark"      # .tfcファイルが入っているディレクトリ
OUTPUT_DIR="./results"        # 結果を保存するディレクトリ

# ★ ここに実行したくないファイル名をスペース区切りで入力してください
EXCLUDE_LIST=(

    "gf2^16_mult.tfc"
    "gf2^32_mult.tfc"
    "gf2^64_mult.tfc"
    "gf2^128_mult.tfc"
    "gf2^256_mult.tfc"
    "cycle_17_3.tfc"
    "hwb11.tfc"
    "hwb12.tfc"
    "mod_adder_1048576.tfc"

)
# ----------------

# 出力ディレクトリの作成
mkdir -p "$OUTPUT_DIR"

echo "Starting optimization..."
echo "Input Directory:  $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "Exclude List:     ${EXCLUDE_LIST[*]}"
echo "---------------------------------------"

count=0
skip_count=0

# .tfc ファイルをループ処理
for tfc_path in "$INPUT_DIR"/*.tfc; do
    # ファイルが存在しない場合のチェック
    [ -e "$tfc_path" ] || continue

    # パスからファイル名のみを取得
    filename=$(basename "$tfc_path")

    # 除外リストに含まれているかチェック
    is_excluded=false
    for excluded in "${EXCLUDE_LIST[@]}"; do
        if [[ "$filename" == "$excluded" ]]; then
            is_excluded=true
            break
        fi
    done

    if $is_excluded; then
        echo "Skipping: $filename (Found in exclude list)"
        ((skip_count++))
        continue
    fi

    # 出力ファイル名の決定 (例: results/adder.out)
    output_file="$OUTPUT_DIR/${filename%.tfc}.out"

    echo "Running: $filename -> $(basename "$output_file")"
    
    # 実行コマンド
    /home/ryo-mtmt/TOpt/bin/TOpt circuit "$tfc_path" -a todd -h 0 > "$output_file"

    if [ $? -eq 0 ]; then
        ((count++))
    else
        echo "Error occurred during: $filename"
    fi
done

echo "---------------------------------------"
echo "Done!"
echo "Optimized: $count files"
echo "Excluded : $skip_count files"