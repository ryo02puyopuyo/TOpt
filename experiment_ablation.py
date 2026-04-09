import os
import subprocess
import re
import csv
from datetime import datetime

TARGET_CIRCUITS = [
    "adder_8.tfc"
]

ALGORITHMS = [
    ("Baseline (Hamming+M4RI)", "todd_m4ri_hamming"),
    ("Exp 000 (No Opts)", "todd_exp_000"),
    ("Exp 100 (Packed Chi)", "todd_exp_100"),
    ("Exp 010 (AA Table)", "todd_exp_010"),
    ("Exp 001 (Memoization)", "todd_exp_001"),
    ("Exp 110 (Packed + AA)", "todd_exp_110"),
    ("Exp 101 (Packed + Memo)", "todd_exp_101"),
    ("Exp 011 (AA + Memo)", "todd_exp_011"),
    ("Exp 111 (All Opts)", "todd_exp_111")
]

def run_experiment():
    benchmarks_dir = "t-par/Benchmarks"
    topt_exec = "./bin/TOpt"
    
    if not os.path.exists(topt_exec):
        print(f"Error: {topt_exec} not found.")
        return

    TIMEOUT_SECONDS = 600
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_file = f"ablation_results_{timestamp}.csv"
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Circuit", "Algorithm", "Init_T", "Fin_T", "Reduced_T", 
                         "Total_Time_ms", "Chi_Time_ms", "AA_Time_ms", "NS_Time_ms", "Status"])
        
        print("=" * 115)
        print(f"{'Circuit Name':<18} | {'Algorithm':<24} | {'FinT':<4} | {'Total(ms)':<10} | {'Chi(ms)':<8} | {'AA(ms)':<6} | {'NS(ms)':<8}")
        print("-" * 115)

        # Regex patterns
        pat_init = re.compile(r"Initial T-count\s*:\s*(\d+)")
        pat_fin = re.compile(r"Final T-count\s*:\s*(\d+)")
        pat_time = re.compile(r"Execution Time\s*:\s*(\d+)\s*ms")
        pat_chi = re.compile(r"Chi calculation\s*:\s*([0-9.]+)\s*ms")
        pat_aa = re.compile(r"AA Table calc\s*:\s*([0-9.]+)\s*ms")
        pat_ns = re.compile(r"Nullspace calc\s*:\s*([0-9.]+)\s*ms")

        # For baseline parsing
        pat_base_ns = re.compile(r"Nullspace time:\s*(\d+)\s*ms")
        pat_base_chi = re.compile(r"Chi time:\s*(\d+)\s*ms")

        for name in TARGET_CIRCUITS:
            file = os.path.join(benchmarks_dir, name)
            if not os.path.exists(file):
                print(f"{name:<18} | NOT FOUND")
                continue
                
            for alg_name, alg_tag in ALGORITHMS:
                args = [topt_exec, "circuit", file, "-a", alg_tag, "-h", "8"]
                
                init_t, fin_t, time_ms, chi_ms, aa_ms, ns_ms = 0, 0, 0, 0.0, 0.0, 0.0
                status = "Success"
                
                try:
                    result = subprocess.run(args, capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
                    
                    # 共通: 初期・最終Tカウント
                    m_init = pat_init.search(result.stdout)
                    init_t = int(m_init.group(1)) if m_init else 0
                    
                    m_fin = pat_fin.search(result.stdout)
                    fin_t = int(m_fin.group(1)) if m_fin else 0
                    
                    if alg_tag.startswith("todd_exp_"):
                        m_time = pat_time.search(result.stdout)
                        time_ms = int(m_time.group(1)) if m_time else 0
                        m_chi = pat_chi.search(result.stdout)
                        chi_ms = float(m_chi.group(1)) if m_chi else 0.0
                        m_aa = pat_aa.search(result.stdout)
                        aa_ms = float(m_aa.group(1)) if m_aa else 0.0
                        m_ns = pat_ns.search(result.stdout)
                        ns_ms = float(m_ns.group(1)) if m_ns else 0.0
                    else:
                        m_base_ns = pat_base_ns.search(result.stdout)
                        ns_ms = float(m_base_ns.group(1)) if m_base_ns else 0.0
                        m_base_chi = pat_base_chi.search(result.stdout)
                        chi_ms = float(m_base_chi.group(1)) if m_base_chi else 0.0
                        
                        m_fail = re.search(r"Fail count = 1", result.stdout)
                        if m_fail:
                            status = "Verification Failed"

                    print(f"{name:<18} | {alg_name:<24} | {fin_t:<4} | {time_ms:<10} | {chi_ms:<8.1f} | {aa_ms:<6.1f} | {ns_ms:<8.1f}")
                    
                except subprocess.TimeoutExpired:
                    status = "Timeout"
                    print(f"{name:<18} | {alg_name:<24} | {'-':<4} | {'>'+str(TIMEOUT_SECONDS*1000):<10} | {'-':<8} | {'-':<6} | {'-':<8}")
                except Exception as e:
                    status = f"Error: {str(e)}"
                    print(f"{name:<18} | {alg_name:<24} | ERROR")
                    
                writer.writerow([name, alg_name, init_t, fin_t, init_t - fin_t, time_ms, chi_ms, aa_ms, ns_ms, status])
                f.flush()
            print("-" * 115)
            
    print(f"Results saved to {csv_file}")

if __name__ == "__main__":
    run_experiment()
