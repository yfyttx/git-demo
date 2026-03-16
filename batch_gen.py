import subprocess
import os

# 指向 build 目录下的可执行文件
executable = "./build/construct" 

# 4 组参数配置: (J, L, P, m2, n2)
configs = [
    (3, 7, 13, 10, 20),
    (3, 7, 31, 10, 30),
    (3, 7, 61, 10, 30),
    (3, 7, 31, 5, 30)
]

if not os.path.exists(executable):
    print(f"错误: 找不到程序 {executable}")
    print("请确认你已经在 build 文件夹内执行过 make，且生成的文件名叫 construct")
    exit(1)

for J, L, P, m2, n2 in configs:
    # 构造命令
    cmd = f"{executable} {J} {L} {P} {m2} {n2}"
    print(f"\n[运行中] {cmd}")
    
    # 执行生成
    subprocess.run(cmd, shell=True)

print("\n--- 所有矩阵对已生成完毕！ ---")
print("生成的 .txt 文件现在应该就在你的 git-demo 根目录下。")
