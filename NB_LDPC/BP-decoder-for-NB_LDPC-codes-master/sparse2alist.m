function sparse2alist(Q, H, file_name)
% sparse2alist: Converts a sparse matrix H into the .alist format.
% This version includes a critical check for file write permissions.

    [M, N] = size(H);
    
    % --- 关键修改：打开文件并检查是否成功 ---
    file_id = fopen(file_name, 'wt'); % 使用 'wt' 以文本模式写入
    
    if file_id == -1
        % 如果 fopen 返回 -1, 说明文件打开失败。立即报错并停止。
        error('sparse2alist:FileOpenError', ...
              ['无法打开或创建文件 "%s" 进行写入。\n', ...
               '请检查以下几点：\n', ...
               '1. 您是否有在当前文件夹写入文件的权限？\n', ...
               '2. 文件夹路径是否真实存在？\n', ...
               '3. 文件是否正被另一个程序锁定或占用？'], ...
              file_name);
    end
    % --- 修改结束 ---

    % --- 从这里开始，后续所有逻辑保持不变，只是将 'file' 替换为 'file_id' ---
    
    % 写入文件头
    if Q == 2
        fprintf(file_id, '%d %d\n', N, M);
    else
        fprintf(file_id, '%d %d %d\n', N, M, Q);
    end
    
    % 计算并写入度数信息
    H_bin = (H > 0); % 创建一个逻辑矩阵来计算度数
    column_weights = sum(full(H_bin), 1); % 按列求和得到每个变量节点的度
    row_weights = sum(full(H_bin), 2)';   % 按行求和得到每个校验节点的度
    c_max = max(column_weights);
    r_max = max(row_weights); 
    
    fprintf(file_id, '%d %d\n', c_max, r_max);
    
    fprintf(file_id, '%d ', column_weights);
    fprintf(file_id, '\n');
    
    fprintf(file_id, '%d ', row_weights);
    fprintf(file_id, '\n');
    
    % 写入变量节点（列）视角列表
    for i = 1:N
        indexes = (find(H(:,i)))';
        values = (full(H(indexes, i)))';
        if Q == 2
            temp = indexes;
        else
            temp = [indexes; values];
        end
        fprintf(file_id, '%d ', temp);
        
        % 添加填充物
        if length(indexes) < c_max
            if Q == 2
                fillers = zeros(1, c_max - length(indexes));
            else
                fillers = zeros(1, 2*(c_max - length(indexes)));
            end
            fprintf(file_id, '%d ', fillers);
        end
        fprintf(file_id, '\n');
    end
    
    % 写入校验节点（行）视角列表
    for i = 1:M
        indexes = find(H(i,:));
        values = full(H(i, indexes));
        if Q == 2
            temp = indexes;
        else
            temp = [indexes; values];
        end
        fprintf(file_id, '%d ', temp);
        
        % 添加填充物
        if length(indexes) < r_max
            if Q == 2
                fillers = zeros(1, r_max - length(indexes));
            else
                fillers = zeros(1, 2*(r_max - length(indexes)));
            end
            fprintf(file_id, '%d ', fillers);
        end
        fprintf(file_id, '\n');
    end
    
    % 关闭文件
    fclose(file_id);
    
end