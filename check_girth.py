import numpy as np

def check_hx_girth(file_path):
    print(f"正在读取矩阵: {file_path}")
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # 解析 Essai 格式头部: N M Q
        header = lines[0].strip().split()
        N, M = int(header[0]), int(header[1])
        
        # 我们只需要记录每个校验节点（行）连接了哪些变量节点（列）
        # 这里的行索引从 0 开始，数据从文件的第 4 行开始
        adj_rows = []
        for line in lines[3:]:
            if not line.strip(): continue
            data = line.strip().split()
            # 提取列索引 (Essai 格式中 col_idx 是每对数据的第一个值，且从 1 开始)
            cols = [int(data[i]) - 1 for i in range(0, len(data), 2)]
            adj_rows.append(set(cols))

        print(f"矩阵维度: {M} x {N}. 开始检测 4 环...")
        
        # 4 环检测逻辑：检查任意两行是否有 2 个及以上的公共列
        found_4_cycle = False
        for i in range(len(adj_rows)):
            for j in range(i + 1, len(adj_rows)):
                # 求交集
                common = adj_rows[i].intersection(adj_rows[j])
                if len(common) >= 2:
                    print(f"发现 4 环！校验行 {i+1} 和 校验行 {j+1} 在列 {list(common)} 处相交。")
                    found_4_cycle = True
                    break
            if found_4_cycle: break
        
        if found_4_cycle:
            print("\n结论: HX 的围长 (Girth) 为 4。")
            print("建议: 这通常是导致 FER 在 p=10^-4 时依然很高 (0.68) 的主因。")
        else:
            print("\n结论: 未发现 4 环，围长 (Girth) 至少为 6。")
            
    except Exception as e:
        print(f"读取或处理失败: {e}")

if __name__ == "__main__":
    # 路径对齐：指向你 build 文件夹下的 HX 文件
    check_hx_girth("build/HX_Safe_GF64.txt")
