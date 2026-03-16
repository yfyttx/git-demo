import subprocess
import re
import os

# --- 实验配置 ---
# 矩阵文件（确保这些文件在当前目录下）
hx_path = "../../build/HX_Safe_GF64.txt"
h3_path = "../../build/H3_sparse.txt" # 用于仿真 HZ 的子块

# 物理错误率扫点
error_rates = [0.001, 0.0005, 0.0001, 0.00005]
# 蒙特卡洛帧数（建议先跑 100 帧看趋势，正式汇报跑 500-1000）
NB_MONTE = 100 
MAX_ITER = 100

# 子块数设置 (由 G2 的行数决定，之前确定是 20)
K2_SUBBLOCKS = 20

def run_sim():
    results = []
    print(f"{'Prob':<10} | {'FER_X (Split)':<15} | {'FER_Z (Global)':<15} | {'Total_FER':<15}")
    print("-" * 60)

    for p in error_rates:
        # --- 对调实验逻辑 ---
        
        # 1. 测 X 错误：改用 H3 矩阵 + 拆分译码 (原本给 Z 用的逻辑)
        # 命令格式: ./essai 帧数 迭代 矩阵 错误率 n_vc n_cv offset nb_oper 子块数
        cmd_x = f"./essai {NB_MONTE} {MAX_ITER} {H3_FILE} {p} 12 12 0.5 15 {K2_SUBBLOCKS}"
        
        # 2. 测 Z 错误：改用 HX 矩阵 + 整体译码 (原本给 X 用的逻辑)
        cmd_z = f"./essai {NB_MONTE} {MAX_ITER} {HX_FILE} {p} 12 12 0.5 15 1"

        # 执行仿真并抓取输出
        def get_fer(cmd):
            try:
                out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode()
                # 假设输出格式包含 "FER = 0.xxxx"
                match = re.search(r"FER\s*=\s*([\d\.]+)", out)
                return float(match.group(1)) if match else 1.0
            except:
                return 1.0

        fer_x = get_fer(cmd_x)
        fer_z = get_fer(cmd_z)
        
        # 量子总错误率判定：只要 X 或 Z 有一个错，整帧就错
        # 统计学近似: Total_FER = 1 - (1 - FER_X) * (1 - FER_Z)
        total_fer = 1 - (1 - fer_x) * (1 - fer_z)
        
        print(f"{p:<10} | {fer_x:<15.4f} | {fer_z:<15.4f} | {total_fer:<15.4f}")
        results.append((p, fer_x, fer_z, total_fer))
    
    return results

if __name__ == "__main__":
    if not os.path.exists("./essai"):
        print("错误：未找到编译好的 essai 可执行程序！")
    else:
        run_sim()
