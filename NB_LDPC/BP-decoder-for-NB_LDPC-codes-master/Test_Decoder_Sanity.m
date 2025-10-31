% --- Test_Decoder_Sanity.m ---
%
% 最终理智检查 (Final Sanity Check)
%
% 目的: 确认 BP-decoder-for-NB_LDPC-codes-master 工具箱
%       (特别是 decode_soft.mexw64) 是否能用它 *自带* 的
%       矩阵文件正常工作。
%
% 步骤:
% 1. 加载一个 *工具箱自带* 的矩阵 (例如 '50x100x55.txt')
% 2. 运行 "全零码字测试" (发送全零码字)
% 3. 检查 FER 是否 < 1.0

    clear all; close all; clc;
    
    % --- 1. 初始化 ---
    s = RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(s);
    
    % --- 2. 核心仿真参数 ---
    fm_range = 0.01:0.01:0.08; % 误翻概率
    FER_error_threshold = 100;
    q = 64; % 假设工具箱的矩阵是为 GF(64) 设计的
    p_bits = log2(q);
    max_decode_iter = 100;
    
    % --- 3. 加载工具箱自带的码本 ---
    % *** 这是关键: 我们使用工具箱自带的文件 ***
    % (确保 '50x100x55.txt' 在您的MATLAB路径下)
    fprintf('正在加载工具箱自带的矩阵 (50x100x55.txt)...\n');
    try
        alist_file = '20x100x328.txt';
        h = alist2sparse(alist_file);
        h = sparse(h);
    catch ME
        fprintf('错误: 找不到 50x100x55.txt 文件。\n');
        fprintf('请确保该文件在 BP-decoder-for-NB_LDPC-codes-master 文件夹中。\n');
        rethrow(ME);
    end
    
    [H, G] = h2g(h, q); 
    [K, N] = size(G);
    
    fprintf('测试码参数: N_syms=%d, K_syms=%d, R=%.4f\n\n', N, K, K/N);
            
    % --- 5. 遍历 f_m 点 ---
    num_fm_points = length(fm_range);
    FER_results = zeros(1, num_fm_points);
    
    for fm_idx = 1:num_fm_points
        current_fm = fm_range(fm_idx);
        FER_count = 0;
        Codewords_total = 0;
        
        fprintf('开始仿真: f_m = %.4f\n', current_fm);
        
        % 我们需要运行足够多的码字来确认
        while (FER_count < FER_error_threshold) && (Codewords_total < 5000)
            Codewords_total = Codewords_total + 1;
            
            % --- 6. 全零码字测试 ---
            encoded_codeword = zeros(1, N); 
            
            % --- 7. 信道模型 (BSC) ---
            bits = de2bi(encoded_codeword, p_bits, 'left-msb')'; 
            bits = bits(:); % 全零比特
            error_bits = (rand(size(bits)) < current_fm);
            received_bits = mod(bits + error_bits, 2);

            % --- 8. 解调 (BSC LLR) ---
            llr_matrix = demodulate_bsc(q, received_bits, current_fm);
            
            % --- 9. 解码 ---
            % (我们只需要一个解码器)
            sparse2alist(q, h, 'simulate_temp.alist');
            [ldpc, ~, ~] = decode_soft(0, 'simulate_temp.alist');
            [~, ~, decoded_syms, ~] = decode_soft(2, ldpc, llr_matrix, max_decode_iter);

            % --- 10. 统计错误 (检查是否恢复为全零) ---
            if ~isequal(encoded_codeword, decoded_syms)
                FER_count = FER_count + 1;
            end
            
            % (减少打印，加快速度)
            if mod(Codewords_total, 500) == 0
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
    semilogx(fm_range, FER_results, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('Flip Probability f_m');
    ylabel('Block Error Rate (FER)');
    title('理智检查: 工具箱自带矩阵 (50x100x55.txt) - 全零测试');
    legend('FER');
    ylim([1e-5 1]);