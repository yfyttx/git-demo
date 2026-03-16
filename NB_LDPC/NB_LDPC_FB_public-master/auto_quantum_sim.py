import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np

# --- 配置参数 ---
# 错误率扫描点 (建议包含一个较大的点如 0.01 确保能看到 FER 下降)
error_rates = [0.001, 0.0008, 0.0005, 0.0003, 0.0001, 0.00005]
nb_monte_carlo = 1000
max_iter = 50
n_vc, n_cv = 12, 12  
offset = 0.5
nb_oper = 15
k2_hz = 20

# 矩阵路径 (保持你原始能跑通的路径)
hx_path = "../../build/HX_Safe_GF64.txt"
h3_path = "../../build/H3_sparse.txt" 

def run_cmd(cmd):
    try:
        print(f"Executing: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        output = result.stdout
        
        # 提取 FER，兼容崩溃后的输出捕获
        fer_match = re.findall(r"FER=([\d.]+)", output)
        if not fer_match:
            fer_match = re.findall(r"Final Quantum FER.*: ([\d.]+)", output)
            
        if fer_match:
            return float(fer_match[-1])
        else:
            print(f"Warning: No FER found! Program output snippet:\n{output[-300:]}")
            return 1.0
            
    except Exception as e:
        print(f"Error running command: {e}")
        return 1.0

results_total = []

print("=== Quantum Performance Simulation Starting ===")

for p in error_rates:
    print(f"\n[Testing Physical Error Rate: {p}]")
    
    # 1. 运行 HX 仿真 (纠正 Z 错误)
    cmd_x = f"./essai {nb_monte_carlo} {max_iter} {hx_path} -{p} {n_vc} {n_cv} {offset} {nb_oper} 1"
    fer_x = run_cmd(cmd_x)
    
    # 2. 运行 HZ 仿真 (纠正 X 错误)
    cmd_z = f"./essai {nb_monte_carlo} {max_iter} {h3_path} -{p} {n_vc} {n_cv} {offset} {nb_oper} {k2_hz}"
    fer_z = run_cmd(cmd_z)
    
    # 3. 计算综合量子 FER
    fer_total = 1 - (1 - fer_x) * (1 - fer_z)
    results_total.append(fer_total)
    
    print(f"Result for p={p}: FER_X={fer_x:.4f}, FER_Z={fer_z:.4f} -> Total_FER={fer_total:.4f}")

# --- 绘图部分 (仅保留 Total FER) ---
plt.figure(figsize=(10, 6))

# 只有在有非零数据时才画 log 曲线，防止报错
if any(v > 0 for v in results_total):
    plt.loglog(error_rates, results_total, 'r-^', linewidth=2.5, markersize=8, label='Total Quantum FER')
    plt.grid(True, which="both", ls="-", alpha=0.5)
else:
    # 如果全是0，画普通的 plot 避免 log 报错
    plt.plot(error_rates, results_total, 'r-^', linewidth=2.5, markersize=8, label='Total Quantum FER')
    plt.grid(True)

plt.xlabel('Physical Error Rate (p)')
plt.ylabel('Frame Error Rate (FER)')
plt.title('Quantum LDPC Total Performance (N=6510, GF64)')
plt.legend()
plt.gca().invert_xaxis() 

plt.savefig('quantum_performance_total_only.png')
print("\nSimulation Finished. Plot saved as 'quantum_performance_total_only.png'")
plt.show()
