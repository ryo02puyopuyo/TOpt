import os
import glob
import subprocess
import re

# 対象とする回路のリスト（コメントアウトや追加で実験対象を自由に変更できます）
TARGET_CIRCUITS = [
    "adder_8.tfc",
    "barenco_tof_3.tfc",
    "barenco_tof_4.tfc",
    "barenco_tof_5.tfc",
    "barenco_tof_10.tfc",
    "csla_mux_3.tfc",
    "csum_mux_9.tfc",
    "gf2^4_mult.tfc",
    "gf2^5_mult.tfc",
    "gf2^6_mult.tfc",
    "gf2^7_mult.tfc",
    "gf2^8_mult.tfc",
    "gf2^9_mult.tfc",
    "grover_5.tfc",
    "ham15-high.tfc",
    "ham15-med.tfc",
    "ham15-low.tfc",
    "hwb6.tfc",
    "hwb8.tfc",
    "mod5_4.tfc",
    "mod_adder_1024.tfc",
    "mod_mult_55.tfc",
    "mod_red_21.tfc",
    "qcla_adder_10.tfc",
    "qcla_com_7.tfc",
    "qcla_mod_7.tfc",
    "qft_4.tfc",
    "rc_adder_6.tfc",
    "tof_3.tfc",
    "tof_4.tfc",
    "tof_5.tfc",
    "tof_10.tfc",
    "vbe_adder_3.tfc"
]

def run_experiment():
    benchmarks_dir = "t-par/Benchmarks"
    topt_exec = "./bin/TOpt"

    if not os.path.exists(topt_exec):
        print(f"Error: {topt_exec} not found. バイナリをコンパイルしてください。")
        return

    # タイムアウト設定（秒）
    TIMEOUT_SECONDS = 600
    
    # アンシラ設定オプション (-h 0 を付ける場合は True)
    USE_H_ZERO = False

    print("=" * 110)
    print(f"{'Circuit Name':<20} | {'Blk':<3} | {'Init 1s':<8} | {'Fin 1s':<8} | {'Init T':<8} | {'Fin T':<8} | {'T-Red':<6} | {'T-Red %':<7}")
    print("-" * 110)

    # 1の数の推移を抽出する正規表現
    weight_pattern = re.compile(r"Total weight: (\d+) -> (\d+) \((\d+) ones removed\)")
    # Tゲート数（列数）の削減を抽出する正規表現
    init_t_pattern = re.compile(r"Initial T-count\s*:\s*(\d+)")
    final_t_pattern = re.compile(r"Final T-count\s*:\s*(\d+)")

    for name in TARGET_CIRCUITS:
        file = os.path.join(benchmarks_dir, name)
        if not os.path.exists(file):
            print(f"{name:<20} | NOT FOUND")
            continue
        
        args = [topt_exec, "circuit", file, "-a", "todd_greedy_preprocess"]
        if USE_H_ZERO:
            args.extend(["-h", "0"])
            
        try:
            result = subprocess.run(
                args,
                capture_output=True,
                text=True,
                timeout=TIMEOUT_SECONDS
            )
            
            # 1の数（密度）の変化
            matches = weight_pattern.findall(result.stdout)
            blocks = len(matches)
            
            if blocks == 0:
                print(f"{name:<20} | {'0':<3} | {'-':<8} | {'-':<8} | {'-':<8} | {'-':<8} | {'-':<6} | -")
                continue
                
            total_init_1s = sum(int(m[0]) for m in matches)
            total_fin_1s = sum(int(m[1]) for m in matches)
            
            # Tゲート数（Phase Polynomial内の列数）の変化
            init_ts = init_t_pattern.findall(result.stdout)
            final_ts = final_t_pattern.findall(result.stdout)
            
            total_init_t = sum(int(t) for t in init_ts)
            total_fin_t = sum(int(t) for t in final_ts)
            
            removed_t = total_init_t - total_fin_t
            t_reduction_pct = (removed_t / total_init_t * 100) if total_init_t > 0 else 0.0
            
            print(f"{name:<20} | {blocks:<3} | {total_init_1s:<8} | {total_fin_1s:<8} | {total_init_t:<8} | {total_fin_t:<8} | {removed_t:<6} | {t_reduction_pct:.2f}%")
            
        except subprocess.TimeoutExpired:
            print(f"{name:<20} | {'T/O':<3} | {'-':<8} | {'-':<8} | {'-':<8} | {'-':<8} | {'-':<6} | Timeout (> {TIMEOUT_SECONDS}s)")
        except Exception as e:
            print(f"{name:<20} | {'ERR':<3} | {'-':<8} | {'-':<8} | {'-':<8} | {'-':<8} | {'-':<6} | Error: {e}")

    print("=" * 110)
    print("Experiment completed.")

if __name__ == "__main__":
    run_experiment()
