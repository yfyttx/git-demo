import sys
import random

def convert_alist_to_kn(input_file, output_file, gf):
    print(f"Converting {input_file} to {output_file} with GF({gf})...")
    
    with open(input_file, 'r') as f:
        # 读取所有 token (整数)，忽略换行符
        all_numbers = [int(x) for line in f for x in line.split()]

    ptr = 0
    
    # 1. Header: N M
    N = all_numbers[ptr]; ptr += 1
    M = all_numbers[ptr]; ptr += 1
    
    # 2. Max Degrees
    max_dv = all_numbers[ptr]; ptr += 1
    max_dc = all_numbers[ptr]; ptr += 1
    
    print(f"Matrix Size: {N}x{M}, max_dv={max_dv}, max_dc={max_dc}")

    # 3. Column Degrees
    col_degrees = []
    for _ in range(N):
        col_degrees.append(all_numbers[ptr])
        ptr += 1
        
    # 4. Row Degrees
    row_degrees = []
    for _ in range(M):
        row_degrees.append(all_numbers[ptr])
        ptr += 1
        
    # 5. Skip Column Connections (跳过列连接部分)
    # 关键修复：Alist 文件通常会被 0 填充到 max_dv 长度
    # 我们必须读取 max_dv 个数，而不是 col_degrees[i] 个数
    for _ in range(N):
        ptr += max_dv 
        
    # 6. Read Row Connections (读取行连接部分)
    # 同样，每一行在文件中都占用了 max_dc 个位置（包含填充的 0）
    kn_rows = []
    
    for m in range(M):
        degree = row_degrees[m]
        indices = []
        
        # 读取当前行的所有数据 (包含填充的 0)
        row_data = []
        for _ in range(max_dc):
            row_data.append(all_numbers[ptr])
            ptr += 1
            
        # 只取前 degree 个有效数据 (非 0)
        # 验证：检查是否确实取到了非 0 值
        valid_indices = row_data[:degree]
        
        # 构建 KN 格式行
        new_line_parts = []
        for idx in valid_indices:
            if idx == 0:
                print(f"Error: Found 0 in valid index region at Row {m+1}!")
            # 随机生成非零系数 [1, GF-1]
            val = random.randint(1, gf - 1) 
            new_line_parts.append(f"{idx} {val}")
        
        kn_rows.append(" ".join(new_line_parts))

    # 7. Write Output
    with open(output_file, 'w') as f:
        f.write(f"{N} {M} {gf}\n")
        f.write(" ".join(map(str, col_degrees)) + "\n")
        f.write(" ".join(map(str, row_degrees)) + "\n")
        for row_str in kn_rows:
            f.write(row_str + "\n")
            
    print("Conversion complete. Padding handled correctly.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 convert_alist.py input.alist output.txt")
        sys.exit(1)
        
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    # 默认 GF(64)
    convert_alist_to_kn(input_path, output_path, 64)
