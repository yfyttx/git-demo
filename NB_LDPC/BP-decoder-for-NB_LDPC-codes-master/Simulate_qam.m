%Copyright(c) 2013, Vasiliy Usatyuk
%All rights reserved.

%(...版权信息...)

% (与之前版本相同的版权声明)

    clear all; close all; clc;
    
    % --- 1. 初始化 ---
    s = RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(s);
    
    % --- 2. 核心仿真参数 ---
    ebn0_range = 4:1:8; 
    FER_error_threshold = 100;
    q = 64; 
    M = 64; 
    nbits = log2(q);
    max_decode_iter = 100;
    
    num_ebn0_points = length(ebn0_range);
    FER_results = zeros(1, num_ebn0_points);
    BER_results = zeros(1, num_ebn0_points);

    % --- 3. 加载码本并填充正确的值 ---
    fprintf('正在从 H_gamma.alist 加载码本结构...\n');
    h = alist2sparse('H_gamma.alist');
    
    fprintf('正在从 H_gamma.txt 加载正确的非二进制矩阵值...\n');
    try
        H_correct_values = load('H_gamma.txt');
        if ~isequal(size(h), size(H_correct_values))
            error('H_gamma.alist 和 H_gamma.txt 的维度不匹配！');
        end
        h = H_correct_values;
        fprintf('成功将正确的非二进制值填充到校验矩阵中。\n\n');
    catch ME
        fprintf('错误：无法加载 H_gamma.txt！请确保该文件在当前目录下。\n');
        rethrow(ME);
    end
    
    % --- 4. 生成编码/解码所需的矩阵 ---
    h = sparse(h); 
    [H, G] = h2g(h, q);
    [K, N] = size(G);
    Rate = K / N;
    fprintf('矩阵参数计算完成: H(%d x %d), G(%d x %d), 码率 R = %.4f\n\n', ...
            size(H,1), size(H,2), K, N, Rate);
    
    % --- 5. 遍历所有信噪比进行仿真 ---
    for ebn0_idx = 1:num_ebn0_points
        current_ebn0 = ebn0_range(ebn0_idx);
        EsN0_dB = current_ebn0 + 10*log10(nbits);
        EsN0_linear = 10^(EsN0_dB / 10);
        noise_var = 1 / (2 * Rate * EsN0_linear);
        sigma = sqrt(noise_var); 
        
        BER_total = 0;
        FER_count = 0;
        Codewords_total = 0;
        
        fprintf('开始仿真: Eb/N0 = %.2f dB (Es/N0 = %.2f dB)\n', current_ebn0, EsN0_dB);
        
        % =======================================================================
        % ========================== 兼容性修复区域开始 ==========================
        % =======================================================================
        
        % --- 关键修复：预先生成自然二进制映射的星座图 ---
        % 我们通过将 0..63 的比特流送入比特模式的 qammod 来实现
        all_symbols_as_bits = de2bi(0:M-1, nbits, 'left-msb')';
        all_symbols_as_bits = all_symbols_as_bits(:);
        constellation_points = qammod(all_symbols_as_bits, M, 'UnitAveragePower', true, 'InputType', 'Bit');

        % =======================================================================

        while FER_count < FER_error_threshold
            Codewords_total = Codewords_total + 1;
            
            % --- 6. 编码 (手动系统编码) ---
            info_symbols = randi([0 q-1], 1, K);
            P_matrix = G(:, 1:N-K);
            gf_info = gf(info_symbols, nbits);
            gf_P = gf(P_matrix, nbits);
            gf_parity = gf_info * gf_P;
            parity_symbols = gf_parity.x;
            encoded_symbols = [parity_symbols, info_symbols];
            
            % --- 7. 调制与信道 ---
            % --- 关键修复：使用通用的比特输入方式实现自然二进制映射 ---
            bits_to_modulate = de2bi(encoded_symbols, nbits, 'left-msb')';
            bits_to_modulate = bits_to_modulate(:);
            tx_signal = qammod(bits_to_modulate, M, 'UnitAveragePower', true, 'InputType', 'Bit');
            
            noise = (randn(size(tx_signal)) + 1i * randn(size(tx_signal))) * sigma;
            rx_signal = tx_signal + noise;
            
            % =======================================================================
            % ========================== 兼容性修复区域结束 ==========================
            % =======================================================================

            % --- 8. 解调 (计算符号级LLR) ---
            LLR_matrix = zeros(q, N);
            for i = 1:N
                distance_sq = abs(rx_signal(i) - constellation_points).^2;
                LLR_matrix(:, i) = -distance_sq' / noise_var; 
            end

            % --- 9. 解码 ---
            sparse2alist(q, H, 'simulate.alist');
            [ldpc, ~, ~] = decode_soft(0, 'simulate.alist');
            [~, ~, decoded_symbols, ~] = decode_soft(2, ldpc, LLR_matrix, max_decode_iter);
            
            % --- 10. 统计错误 ---
            decoded_info = decoded_symbols(N-K+1 : N)';
            
            if ~isequal(double(info_symbols), double(decoded_info))
                FER_count = FER_count + 1;
                original_bits = de2bi(info_symbols, nbits)';
                decoded_bits = de2bi(decoded_info, nbits)';
                bit_errors = sum(original_bits(:) ~= decoded_bits(:));
                BER_total = BER_total + bit_errors;
            end
            
            if mod(Codewords_total, 50) == 0
                fprintf('  -> Codewords: %d, Frame Errors: %d, FER: %.4f\n', ...
                        Codewords_total, FER_count, FER_count/Codewords_total);
            end
        end
        
        % --- 11. 保存当前信噪比点的结果 ---
        FER_results(ebn0_idx) = FER_count / Codewords_total;
        BER_results(ebn0_idx) = BER_total / (Codewords_total * K * nbits);
        
        fprintf('完成: Eb/N0 = %.2f dB, FER = %.5f, BER = %.2e\n\n', ...
                current_ebn0, FER_results(ebn0_idx), BER_results(ebn0_idx));
    end
    
    % --- 12. 绘制性能曲线 ---
    figure;
    semilogy(ebn0_range, FER_results, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'FER');
    hold on;
    semilogy(ebn0_range, BER_results, 's--', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BER');
    grid on;
    xlabel('Eb/N0 (dB)');
    ylabel('Error Rate');
    title(sprintf('NB-LDPC Performance over AWGN with 64-QAM (GF(%d), R=%.2f)', q, Rate));
    legend('show');
    ylim([1e-5 1]);