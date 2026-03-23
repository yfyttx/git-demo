import os, subprocess, re
import matplotlib.pyplot as plt

# 重点探测门限区域
p_list = [0.15, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
fer_results = []

print("=== Starting TIT Split Decoding Simulation ===")
os.system("make clean && make essai")

for p in p_list:
    # 模拟 1000 帧以获取稳定数据
    cmd = f"./essai 1000 50 ../../build/H1_sparse.txt ../../build/G2_sparse.txt {p} 0.5 50"
    print(f"Testing p={p}...", end="", flush=True)
    
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    match = re.search(r"Final Quantum FER: ([\d.]+)", res.stdout)
    
    fer = float(match.group(1)) if match else 1.0
    fer_results.append(fer)
    print(f" FER={fer}")

# 绘制瀑布图
plt.figure(figsize=(10, 7))
plt.semilogy(p_list, fer_results, 'r-o', linewidth=2, label='TIT Split Decoder (GF64)')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.xlabel('Physical Error Rate (p)')
plt.ylabel('Frame Error Rate (FER)')
plt.title('Quantum LDPC Performance: TIT Reproduction')
plt.legend()
plt.savefig('tit_final_waterfall.png')
print("\n[OK] Plot saved as tit_final_waterfall.png")
