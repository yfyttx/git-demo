% --- Simulate_CSS_Quantum.m ---
%
% == 全零码字测试版 ==
%
% 功能：通过只发送全零码字来隔离测试解码器 (decode_soft)。
% 这会绕过 MATLAB 的 ldpc_encode.m 和 h2g.m，
% 从而消除因 GF(q) 本原多项式不匹配而导致的编码错误。
%

    clear all; close all; clc;
    
    % --- 1. 初始化 ---
    s = RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(s);
    
    % --- 2. 核心仿真参数 ---
    fm_range = 0.01:0.01:0.08; % 误翻概率 f_m 范围
    FER_error_threshold = 100; % 累计100个块错误后停止
    q = 64; % 域大小 GF(2^6)
    p_bits = log2(q); % p = 6
    max_decode_iter = 100; % 解码器最大迭代次数
    
    % --- 3. 加载码本 (H_C 和 H_D) ---
    
    % 加载 H_C (H_gamma) 用于 Z 错误
    fprintf('正在加载 H_gamma (用于 Z 错误)...\n');
    alist_file_C = 'H_gamma.alist';
    h_C_skel = alist2sparse(alist_file_C);
    h_C = sparse(load('../NB-LDPC-toolbox-master/H_gamma.txt')); % 从 C++ 加载指数值
    if ~isequal(size(h_C_skel), size(h_C))
        error('H_gamma.alist 和 H_gamma.txt 维度不匹配!');
    end

    % 加载 H_D (H_delta) 用于 X 错误
    fprintf('正在加载 H_delta (用于 X 错误)...\n');
    alist_file_D = 'H_delta.alist';
    h_D_skel = alist2sparse(alist_file_D);
    h_D = sparse(load('../NB-LDPC-toolbox-master/H_delta.txt')); % 确保 C++ 保存了此文件
    if ~isequal(size(h_D_skel), size(h_D))
        error('H_delta.alist 和 H_delta.txt 维度不匹配!');
    end
    
    % --- 4. 生成编码/解码所需的矩阵 ---
    % 注意：我们不再使用 G_C 和 G_D，但我们保留 N_C 和 N_D 的计算
    [H_C, G_C] = h2g(h_C, q); 
    [K_C, N_C] = size(G_C);
    
    [H_D, G_D] = h2g(h_D, q);
    [K_D, N_D] = size(G_D);
    
    if N_C ~= N_D
        error('H_C 和 H_D 的码长 (N) 不相等!');
    end
    N = N_C; % 总符号码长
    
    K_logical = K_C + K_D - N; 
    Rate_Q = K_logical / N; 
    
    fprintf('量子CSS码参数: N_syms=%d, K_logical_syms=%d, R_Q = %.4f\n\n', ...
            N, K_logical, Rate_Q);
            
    % --- 5. 遍历所有 f_m 点进行仿真 ---
    num_fm_points = length(fm_range);
    FER_results = zeros(1, num_fm_points);
    
    for fm_idx = 1:num_fm_points
        current_fm = fm_range(fm_idx);
        
        FER_count = 0;
        Codewords_total = 0;
        
        fprintf('开始仿真: f_m = %.4f\n', current_fm);
        
        while FER_count < FER_error_threshold
            Codewords_total = Codewords_total + 1;
            
            % --- 6. 编码 (已修改为全零码测试) ---
            % 我们不再使用 ldpc_encode.m 来避免 GF(q) 不匹配问题
            
            % 线性码的全零信息 [0,...] 对应的码字必定是全零码字 [0,...]
            encoded_C = zeros(1, N_C); 
            encoded_D = zeros(1, N_D);
            % -------------------------------------
            
            % --- 7. 信道模型 (BSC) ---
            % Z 错误仿真
            bits_C = de2bi(encoded_C, p_bits, 'left-msb')'; 
            bits_C = bits_C(:); % 全零比特
            z_error_bits = (rand(size(bits_C)) < current_fm); % 引入 Z 错误
            z_received_bits = mod(bits_C + z_error_bits, 2); % 仅包含错误
            
            % X 错误仿真
            bits_D = de2bi(encoded_D, p_bits, 'left-msb')';
            bits_D = bits_D(:); % 全零比特
            x_error_bits = (rand(size(bits_D)) < current_fm); % 引入 X 错误
            x_received_bits = mod(bits_D + x_error_bits, 2); % 仅包含错误

            % --- 8. 解调 (BSC) ---
            % 根据 f_m 计算 LLR
            z_llr_matrix = demodulate_bpsk(q, z_received_bits, current_fm);
            x_llr_matrix = demodulate_bpsk(q, x_received_bits, current_fm);
            
            % --- 9. 解码 (双解码器) ---
            sparse2alist(q, h_C, 'simulate_C.alist');
            sparse2alist(q, h_D, 'simulate_D.alist');

            [ldpc_C, ~, ~] = decode_soft(0, 'simulate_C.alist');
            [~, ~, decoded_syms_C, ~] = decode_soft(2, ldpc_C, z_llr_matrix, max_decode_iter);
            
            [ldpc_D, ~, ~] = decode_soft(0, 'simulate_D.alist');
            [~, ~, decoded_syms_D, ~] = decode_soft(2, ldpc_D, x_llr_matrix, max_decode_iter);

            % --- 10. 统计错误 (检查是否恢复为全零) ---
            
            % 比较“发送的完整码字”(现在是全零)与“解码的完整码字”
            z_error_uncorrected = ~isequal(encoded_C, decoded_syms_C);
            x_error_uncorrected = ~isequal(encoded_D, decoded_syms_D);

            % 块错误 (FER): 任何一个解码器失败
            if z_error_uncorrected || x_error_uncorrected
                FER_count = FER_count + 1;
            end
            
            if mod(Codewords_total, 50) == 0
                fprintf('  -> Codewords: %d, Frame Errors: %d, FER: %.4f\n', ...
                        Codewords_total, FER_count, FER_count/Codewords_total);
            end
        end
        
        % --- 11. 保存结果 ---
        FER_results(fm_idx) = FER_count / Codewords_total;
        
        fprintf('完成: f_m = %.4f, FER = %.5f\n\n', ...
                current_fm, FER_results(fm_idx));
    end
    
    % --- 12. 绘制性能曲线 ---
    figure;
    semilogy(fm_range, FER_results, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('Flip Probability f_m'); % X轴是 f_m
    ylabel('Block Error Rate (FER)');
    title(sprintf('Quantum CSS Code Performance (GF(%d), R_Q=%.2f) - ALL ZERO TEST', q, Rate_Q));
    legend('FER (X or Z error)');
    ylim([1e-5 1]);