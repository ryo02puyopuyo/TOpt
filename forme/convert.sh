#!/bin/bash

# このスクリプトがあるディレクトリの絶対パスを取得し、Pythonスクリプトを特定
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PY_CONVERTER="$SCRIPT_DIR/convert.py"

if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_directory_path>"
    exit 1
fi

TARGET_DIR=$1

# 入力ディレクトリの存在確認
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Input directory '$TARGET_DIR' not found."
    exit 1
fi

# Pythonスクリプトの存在確認
if [ ! -f "$PY_CONVERTER" ]; then
    echo "Error: $PY_CONVERTER not found."
    exit 1
fi

count=0
echo "Converting files from: $TARGET_DIR"
echo "Output directory: $(pwd)"
echo "---------------------------------------"

for qc_file in "$TARGET_DIR"/*.qc; do
    [ -e "$qc_file" ] || continue

    # ファイル名のみを抽出 (例: /path/to/adder.qc -> adder)
    base_name=$(basename "$qc_file" .qc)
    
    # カレントディレクトリに出力ファイルパスを設定
    tfc_file="./${base_name}.tfc"

    echo "Processing: $(basename "$qc_file") -> $tfc_file"
    python3 "$PY_CONVERTER" "$qc_file" "$tfc_file"

    if [ $? -eq 0 ]; then
        ((count++))
    else
        echo "Failed to convert: $qc_file"
    fi
done

echo "---------------------------------------"
echo "Done! $count files have been saved in the current directory."