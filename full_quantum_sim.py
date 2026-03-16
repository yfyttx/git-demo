import subprocess
import re
import matplotlib.pyplot as plt
import os

# ==========================================
# 1. 实验配置 - 必须与 batch_gen.py 的参数严格对应
# ==========================================
# 格式: (显示标签, P参数, n2参数, k2子块数)
# 注意：k2 = n2 - m2 (根据你的 main.cpp 逻辑)
CONFIGS = [
    ("N=1820 (P=13, n2=20)", 13, 20, 10),
    ("N=6510 (P=31, n2=30)", 31, 30, 20),
    ("N=12810 (P=61, n2=30)", 61, 30, 20),
    ("N=6510_HighRate (m2=5)", 31, 30, 25)
]

# 之前跑通的 4 个错误率点
ERROR_RATES = [0.02, 0.015, 0.01, 0.005, 0.001]

# 译码器路径（假设 essai 在根目录，如果在 build 里请改为 ./build/essai）
DECODER = "./essai"

def run_decode(file_name, p, k2_val):
    """调用译码器并提取 FER，即使程序崩溃也尝试读取输出"""
    cmd = f"{DECODER} 100 10 {file_name} {p} 12 12 0.5 15 {k2_val}"
    try:
        # 使用 check_output 但捕获异常，因为程序崩溃会抛出 CalledProcessError
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        # 无论程序是否正常退出，都从 stdout 中找 FER
        output = process.stdout
        match = re.search(r"Final Quantum FER.*?([\d\.]+)", output)
        if match:
            return float(match.group(1))
        
        # 如果第一种没匹配到，试试匹配普通的 FER=
        match = re.search(r"FER=([\d\.]+)", output)
        return float(match.group(1)) if match else 1.0
    except Exception as e:
        print(f"执行出错: {e}")
        return 1.0



def main():
    plt.figure(figsize=(10, 7))
    
    for label, p_val, n2_val, k2 in CONFIGS:
        print(f"\n>>> 正在仿真矩阵: {label}")
        
        # 对应 batch_gen.py 生成的文件名格式
        hx_file = f"HX_P{p_val}_n{n2_val}.txt"
        h3_file = f"H3_P{p_val}_n{n2_val}.txt"
        
        if not os.path.exists(hx_file) or not os.path.exists(h3_file):
            print(f"  [跳过] 找不到文件: {hx_file} 或 {h3_file}")
            continue

        total_fers = []
        for p in ERROR_RATES:
            print(f"  正在测试 p={p}...", end=" ", flush=True)
            
            # 1. 纠正 X 错误 (利用 H3 拆分译码)
            fer_x = run_decode(h3_file, p, k2)
            
            # 2. 纠正 Z 错误 (利用 HX 整体译码)
            fer_z = run_decode(hx_file, p, 1)
            
            # 3. 计算量子总误帧率
            total_fer = 1 - (1 - fer_x) * (1 - fer_z)
            total_fers.append(total_fer)
            print(f"FER_total={total_fer:.4f}")

        plt.plot(ERROR_RATES, total_fers, marker='s', label=label)

    # 图表美化
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Physical Error Rate (p)')
    plt.ylabel('Total Quantum FER')
    plt.title('Performance of Non-binary QLDPC Codes (Various Sizes)')
    plt.grid(True, which="both", ls="-", alpha=0.3)
    plt.legend()
    plt.gca().invert_xaxis()
    
    plt.savefig('quantum_final_comparison.png')
    print("\n[成功] 最终对比图已保存为 quantum_final_comparison.png")

if __name__ == "__main__":
    main()
