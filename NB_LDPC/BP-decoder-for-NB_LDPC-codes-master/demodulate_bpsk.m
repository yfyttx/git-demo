% --- demodulate_bsc.m ---
% 版权归属：此文件基于您提供的 demodulate_bpsk.m 修改而来
% 
% 功能：将接收到的比特流（来自误翻概率为 fm 的 BSC 信道）
% 转换为非二进制解码器所需的 LLR 矩阵 (N x q)。
%
% 输入:
%   q   - 域大小 (例如 64)
%   received_bits - 接收到的比特列向量 (N*p x 1)
%   fm  - 信道误翻概率 (例如 0.01)
%
% 输出:
%   llr_matrix - 符号 LLR 矩阵 (N x q)

function [llr_matrix] = demodulate_bsc(q, received_bits, fm)
    
    % --- 1. 参数计算 ---
    p = log2(q); % 每个符号的比特数 (例如 log2(64) = 6)
    
    % 确保 fm 在有效范围内，避免 log(0)
    if fm <= 0
        fm = 1e-9;
    end
    if fm >= 1
        fm = 1 - 1e-9;
    end
    
    received_bits = received_bits(:); % 强制转为列向量
    n_bits = length(received_bits);
    n_syms = n_bits / p; % 符号数 (N)
    
    % --- 2. 维度检查 (来自您的原始代码) ---
    if mod(n_syms, 1) ~= 0
        warning('接收信号长度 %d 不能被 p=%d 整除。', n_bits, p);
        % 自动截断/补零
        ideal_bits = floor(n_syms) * p;
        if n_bits > ideal_bits
            received_bits = received_bits(1:ideal_bits);
        else
            received_bits = [received_bits; zeros(ideal_bits - n_bits, 1)];
        end
        n_bits = ideal_bits;
        n_syms = n_bits / p;
    end
    
    % --- 3. 计算比特LLR (核心修改) ---
    % 这是标准 BSC 信道的 LLR 计算公式
    % LLR = log( P(bit=0|rx) / P(bit=1|rx) )
    % 如果 rx=0, LLR = log((1-fm)/fm)
    % 如果 rx=1, LLR = log(fm/(1-fm))
    % 这可以合并为: LLR = (1 - 2*rx) * log((1-fm)/fm)
    
    base_llr = log((1-fm) / fm);
    bit_llr = base_llr * (1 - 2*received_bits); % (N*p x 1)
    
    % 重塑为“比特数 x 符号数”矩阵
    bit_llr_matrix = reshape(bit_llr, p, n_syms); % (p x N)
    
    % --- 4. 比特LLR 转 符号LLR (使用您文件中的逻辑) ---
    % 此逻辑将 (p x N) 的比特LLR 转换为 (N x q) 的符号LLR矩阵
    % 以匹配 decode_soft 解码器的输入要求
    
    llr_matrix = zeros(n_syms, q); % (N x q)
    
    % 预先计算所有 q 个符号对应的 p 个比特
    all_vals = 0:(q-1);
    all_bits = de2bi(all_vals, p, 'left-msb')'; % (p x q) 矩阵
                                             % 每一列是一个符号的比特
    
    % (2*all_bits - 1) 将 0/1 映射到 +1/-1
    bit_mapping = 2*all_bits - 1; 
    
    % 遍历每个符号
    for i = 1:n_syms
       % 提取当前符号的 p 个比特LLR
       current_bit_llrs = bit_llr_matrix(:, i); % (p x 1)
       
       % (1 x p) * (p x q) -> (1 x q)
       % 这实现了您在 demodulate_bpsk.m 中的 sum() 逻辑
       llr_matrix(i, :) = current_bit_llrs' * bit_mapping;
    end
end