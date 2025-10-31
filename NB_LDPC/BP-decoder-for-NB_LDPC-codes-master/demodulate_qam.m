function LLR_matrix = demodulate_qam(hMod, rx_signal, sigma)
% demodulate_qam: Calculates symbol LLRs for a QAM signal in AWGN.
% This is a MODERNIZED version to replace the obsolete original file.

    % 从调制器对象中获取参数
    M = hMod.M;
    nbits = log2(M);
    N = length(rx_signal); % 码长 (符号数)

    % 获取星座图点
    % 关键：确保使用与调制时相同的 'SymbolOrder'
    if strcmp(hMod.SymbolOrder, 'binary')
        order = 'bin';
    else
        order = 'gray'; % 默认为格雷码
    end
    constellation_points = qammod(0:M-1, M, 'UnitAveragePower', true, 'InputType', 'Integer', 'SymbolOrder', order);
    
    % 计算符号级 LLR
    LLR_matrix = zeros(M, N);
    noise_var = sigma^2; % 噪声方差

    for i = 1:N % 遍历每个接收到的符号
        % 计算接收信号与每个星座图点的距离的平方
        distance_sq = abs(rx_signal(i) - constellation_points).^2;
        % LLR 正比于 -distance^2 / (2 * sigma_per_dimension^2)
        % 对于复高斯噪声，总方差 sigma^2 = 2 * sigma_per_dimension^2
        LLR_matrix(:, i) = -distance_sq' / noise_var; 
    end
end