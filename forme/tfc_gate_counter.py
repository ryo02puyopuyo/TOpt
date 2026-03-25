import os
import re

def count_tfc_gates(file_path):
    """
    TFCファイルを解析し、CCZを 7T + 6CNOT に分解した後のゲート数をカウントする
    """
    total_gates = 0
    one_qubit_gates = 0
    two_qubit_gates = 0
    t_gates = 0

    # ゲート判定用設定
    # 1量子ビットゲートの定義
    single_qubit_names = {'h', 't', 't*', 's', 's*', 'z', 'x', 'y', 'p', 'p*'}
    # 2量子ビットゲートの定義
    double_qubit_names = {'cnot', 'cz', 'cs', 'cs*', 't2'}

    try:
        with open(file_path, 'r') as f:
            content = f.read()
            # BEGINとENDの間のブロックを抽出
            body_match = re.search(r'BEGIN\s+(.*?)\s+END', content, re.DOTALL | re.IGNORECASE)
            if not body_match:
                return None
            
            lines = body_match.group(1).splitlines()
            
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # ゲート名と引数を分離（スペースまたはカンマ区切りに対応）
                parts = re.split(r'[\s,]+', line)
                gate_name = parts[0].lower()
                args = parts[1:]
                arg_count = len(args)

                # --- ゲート判定ロジック ---

                # 1. CCZゲート (Zの3引数, CCZ, t3, または3引数のtof)
                if (gate_name == 'z' and arg_count == 3) or \
                   (gate_name in ['ccz', 't3']) or \
                   (gate_name in ['tof', 't'] and arg_count == 3):
                    # 分解ルール: 7T + 6CNOT
                    t_gates += 7
                    one_qubit_gates += 7  # Tゲートは1量子ビットゲート
                    two_qubit_gates += 6  # CNOTは2量子ビットゲート
                    total_gates += (7 + 6)

                # 2. 2量子ビットゲート
                elif gate_name in double_qubit_names or \
                     (gate_name in ['tof', 't'] and arg_count == 2):
                    two_qubit_gates += 1
                    total_gates += 1

                # 3. 1量子ビットゲート
                elif gate_name in single_qubit_names:
                    one_qubit_gates += 1
                    total_gates += 1
                    # Tゲート単体の場合のカウント
                    if gate_name in ['t', 't*']:
                        t_gates += 1
                
                # 4. その他（多量子ビットToffoliなど）
                else:
                    # 必要に応じてここに追加の分解ロジックを記述
                    pass

        return {
            "total": total_gates,
            "one_qubit": one_qubit_gates,
            "two_qubit": two_qubit_gates,
            "t_count": t_gates
        }
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None

def main(directory_path):
    if not os.path.isdir(directory_path):
        print(f"Error: {directory_path} は有効なディレクトリではありません。")
        return

    print(f"{'File Name':<25} | {'Total':<6} | {'1-Qbt':<6} | {'2-Qbt':<6} | {'T-Gate':<6}")
    print("-" * 65)

    files = [f for f in os.listdir(directory_path) if f.endswith('.tfc')]
    if not files:
        print("指定されたディレクトリに .tfc ファイルは見つかりませんでした。")
        return

    for filename in sorted(files):
        path = os.path.join(directory_path, filename)
        res = count_tfc_gates(path)
        
        if res:
            print(f"{filename:<25} | {res['total']:<6} | {res['one_qubit']:<6} | {res['two_qubit']:<6} | {res['t_count']:<6}")

if __name__ == "__main__":
    # 解析したいディレクトリパスを指定してください
    target_dir = "/home/ryo-mtmt/TOpt/forme/tpar_benchmark" 
    main(target_dir)