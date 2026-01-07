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
# 确保这里指向你正确生成的矩阵文件
MATRIX_FILE = "../../build/H_gamma.txt"

# 仿真参数
NB_MONTE_CARLO = 1000   # 目标帧数 (高信噪比下很快，1000帧没问题)
NB_ITER_MAX    = 50     
N_VC           = 2      
N_CV           = 6      # 你的矩阵行重是 L=6
OFFSET         = 0.5    
NB_OPER        = 64     

# --- [关键修改] 扫描 AWGN 信噪比 (3~18 dB, 步长 1) ---
# range(3, 19) 会生成 3 到 18 的整数列表
SNR_POINTS = list(range(3, 19)) 

# ==========================================
# 2. 仿真函数 (支持实时输出)
# ==========================================
def run_simulation(snr):
    # 构建命令 (AWGN 传入正数 dB)
    cmd = [
        ESSAI_EXE,
        str(NB_MONTE_CARLO),
        str(NB_ITER_MAX),
        MATRIX_FILE,
        str(snr),  # 正数表示 dB
        str(N_VC),
        str(N_CV),
        str(OFFSET),
        str(NB_OPER)
    ]
    
    print(f"\n>>> Simulating AWGN: Eb/No = {snr} dB ...")
    
    fer = None
    full_output = ""
    
    try:
        # 启动进程
        process = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            text=True, 
            bufsize=0
        )
        
        # 实时读取输出 (防止死锁)
        while True:
            char = process.stdout.read(1)
            if not char and process.poll() is not None:
                break
            if char:
                # 只打印关键信息行，避免刷屏太快
                full_output += char
                # 如果遇到回车符或换行符，且行内包含 "FER=", 则打印该行
                if char in ['\r', '\n'] and "FER=" in full_output.split('\n')[-1]:
                     sys.stdout.write(full_output.split('\n')[-1] + '\r')
                     sys.stdout.flush()

        # 提取最终 FER
        fer_matches = re.findall(r"FER=\s*[\d]+\s*/\s*[\d]+\s*=\s*([\d\.]+)", full_output)
        
        if fer_matches:
            fer = float(fer_matches[-1])
            print(f"\nCompleted {snr} dB. FER = {fer}")
            return fer
        else:
            print(f"\nError: Could not parse FER for {snr} dB.")
            return None
            
    except Exception as e:
        print(f"\nSimulation failed: {e}")
        return None

# ==========================================
# 3. 主程序
# ==========================================
def main():
    if not os.path.exists(ESSAI_EXE):
        print("Error: Executable not found.")
        return

    results_fer = []
    valid_snr   = []
    zero_error_count = 0

    print("========================================")
    print("   AWGN Waterfall Test (3 ~ 18 dB)      ")
    print("========================================")

    for snr in SNR_POINTS:
        fer = run_simulation(snr)
        
        if fer is not None:
            # 记录数据
            # 如果 FER 为 0，为了画对数图不报错，我们记录一个极小值 (如 1e-7) 或者在绘图时处理
            # 这里我们原样记录，画图时过滤
            results_fer.append(fer)
            valid_snr.append(snr)
            
            # --- 智能停止策略 ---
            # 如果连续两个点 FER 都是 0 (完美解码)，就没必要跑更离谱的高信噪比了
            if fer == 0.0:
                zero_error_count += 1
                if zero_error_count >= 2:
                    print("Performance is perfect (0 errors). Stopping early to save time.")
                    break
            else:
                zero_error_count = 0

    if not valid_snr:
        print("No data collected.")
        return

    # ==========================================
    # 4. 绘图 (AWGN 风格)
    # ==========================================
    print("\nGenerating AWGN Plot...")
    plt.figure(figsize=(10, 8))
    
    # 预处理数据：把 0 替换成非 0 的极小值以便 log 坐标显示，或者直接不画 0 点
    plot_snr = []
    plot_fer = []
    for s, f in zip(valid_snr, results_fer):
        if f > 0:
            plot_snr.append(s)
            plot_fer.append(f)
        else:
            # 可选：画一个极低的点表示0错误
            plot_snr.append(s)
            plot_fer.append(1e-6) 

    plt.semilogy(plot_snr, plot_fer, 'bo-', linewidth=2, markersize=8, label='Block Error Rate (AWGN)')

    plt.grid(True, which="both", linestyle='--', alpha=0.6)
    
    # 设置标签
    plt.xlabel('Eb/No (dB)', fontsize=14)
    plt.ylabel('Block Error Rate (FER)', fontsize=14)
    plt.title('NB-LDPC Performance over AWGN Channel', fontsize=16)
    plt.legend(loc='lower left', fontsize=12)
    
    # 设置坐标轴
    plt.ylim([1e-6, 1.1])
    plt.xlim([2.5, max(valid_snr) + 0.5])

    output_png = "awgn_waterfall.png"
    plt.savefig(output_png, dpi=300)
    print(f"Success! Plot saved to: {output_png}")

if __name__ == "__main__":
    main()
