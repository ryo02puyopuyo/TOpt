#!/bin/bash

PROGRAM="/home/rest/forpyzx/TOpt/bin/TOpt gsm /home/rest/forpyzx/TOpt/data_for_TOPT/gf27/topt_input_0007.gsm"   # 実行したいプログラム
OUTFILE="1211_3.txt"       # 出力先ファイル


for i in {1..5}; do
    echo "===== Run $i =====" >> "$OUTFILE"
    $PROGRAM >> "$OUTFILE"
    echo "" >> "$OUTFILE"
done