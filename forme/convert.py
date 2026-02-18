import sys
import os

def convert_qc_to_tfc(input_path, output_path):
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    with open(input_path, 'r') as f:
        lines = f.readlines()

    tfc_lines = []
    variables = []
    in_begin_block = False

    for line in lines:
        line = line.strip()
        # コメント行の処理
        if not line or line.startswith('#'):
            continue
        
        # ヘッダー処理 (.v)
        if line.startswith('.v'):
            variables = line.split()[1:]
            var_str = ",".join(variables)
            tfc_lines.append(f".v {var_str}")
            # .i と .o も同じ変数リストで生成（パーサーの期待値）
            tfc_lines.append(f".i {var_str}")
            tfc_lines.append(f".o {var_str}")
            continue

        if line == "BEGIN":
            in_begin_block = True
            tfc_lines.append("BEGIN")
            continue
        
        if line == "END":
            tfc_lines.append("END")
            break

        # ゲートの変換（BEGIN〜END間）
        if in_begin_block:
            parts = line.split()
            gate_name = parts[0].upper() # パーサーに合わせて大文字化
            args = parts[1:]
            
            # 1. Hadamard 
            if gate_name == 'H':
                # ロジック: H <target>
                tfc_lines.append(f"H {args[0]}")
            
            # 2. Toffoli (tof / t3 / t2)
            elif gate_name in ['TOF', 'TOFFOLI', 'T3', 'T2']:
                if len(args) == 3:
                    # 3入力: t3 <control1>,<control2>,<target>
                    tfc_lines.append(f"t3 {','.join(args)}")
                elif len(args) == 2:
                    # 2入力: CNOT または t2
                    tfc_lines.append(f"CNOT {','.join(args)}")
            
            # 3. Z / CCZ
            elif gate_name == 'Z':
                # 引数の数でパーサーが Z, CZ, CCZ を自動判別
                tfc_lines.append(f"Z {','.join(args)}")

            # 4. CNOT
            elif gate_name == 'CNOT':
                tfc_lines.append(f"CNOT {','.join(args)}")

    # ファイル書き出し
    with open(output_path, 'w') as f:
        f.write("\n".join(tfc_lines))
    
    print(f"Done: {output_path} has been generated.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python convert2tfc.py input.qc output.tfc")
    else:
        convert_qc_to_tfc(sys.argv[1], sys.argv[2])