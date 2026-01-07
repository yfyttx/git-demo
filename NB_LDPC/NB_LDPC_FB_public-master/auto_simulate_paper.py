import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# ==========================================
# 1. 配置区域
# ==========================================
ESSAI_EXE = "./essai"
MATRIX_FILE = "../../build/H_gamma.txt"

# 仿真参数
NB_MONTE_CARLO = 2000   # 跑多一点，曲线更平滑
NB_ITER_MAX    = 50
N_VC           = 2
N_CV           = 6      # 你的 L=6
OFFSET         = 0.5
NB_OPER        = 64

# --- [关键] 扫描概率点 (Flip Probability) ---
# 我们扫描从 p=0.10 到 p=0.01
# 注意：必须写成负数，传给 C 程序
PROBABILITY_POINTS = [-0.10, -0.08, -0.06, -0.05, -0.04, -0.03, -0.02, -0.015, -0.01]

# ==========================================
# 2. 仿真函数
# ==========================================
def run_simulation(val):
    cmd = [
        ESSAI_EXE, str(NB_MONTE_CARLO), str(NB_ITER_MAX), MATRIX_FILE,
        str(val), str(N_VC), str(N_CV), str(OFFSET), str(NB_OPER)
    ]
    
    p_real = abs(val)
    print(f"\n>>> Simulating p = {p_real:.3f} ...")
    
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=0)
        
        full_output = ""
        while True:
            char = process.stdout.read(1)
            if not char and process.poll() is not None: break
            if char: full_output += char
                
        # 提取 FER
        fer_matches = re.findall(r"FER=\s*[\d]+\s*/\s*[\d]+\s*=\s*([\d\.]+)", full_output)
        if fer_matches:
            fer = float(fer_matches[-1])
            print(f"    Done. FER = {fer:.6f}")
            return fer
        return None
            
    except Exception as e:
        print(f"Error: {e}")
        return None

# ==========================================
# 3. 主程序
# ==========================================
def main():
    if not os.path.exists(ESSAI_EXE): return

    results_fer = []
    valid_prob  = []

    print("=== Starting Paper Replication (N=42) ===")

    for p_val in PROBABILITY_POINTS:
        fer = run_simulation(p_val)
        if fer is not None:
            results_fer.append(fer)
            valid_prob.append(abs(p_val))

    if not valid_prob: return

    # ==========================================
    # 4. 绘图 (Paper Style)
    # ==========================================
    print("\nGenerating Plot...")
    plt.figure(figsize=(10, 8))
    
    # 绘制 FER vs Probability
    plt.loglog(valid_prob, results_fer, 'bo-', linewidth=2, markersize=8, label=f'Block Error Rate (N={N_VC*21})')

    plt.grid(True, which="both", linestyle='--', alpha=0.6)
    plt.xlabel('Flip Probability (p)', fontsize=14)
    plt.ylabel('Block Error Rate', fontsize=14)
    
    # 反转 X 轴 (让 0.01 在右边，或者根据论文习惯调整)
    # plt.gca().invert_xaxis() 

    plt.title(f'Performance over Depolarizing Channel', fontsize=16)
    plt.legend(loc='lower left', fontsize=12)
    
    output_name = "paper_figure_reproduction.png"
    plt.savefig(output_name, dpi=300)
    print(f"Done! Saved as {output_name}")

if __name__ == "__main__":
    main()
