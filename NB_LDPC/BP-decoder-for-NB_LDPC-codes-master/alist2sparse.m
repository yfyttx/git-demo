function H = alist2sparse(file_name)
% alist2sparse: Converts an .alist file into a MATLAB sparse matrix.
% (已修正版本：增加了文件末尾检查)

    fid = fopen(file_name, 'r');
    if fid == -1
        error('alist2sparse:FileOpenError', '无法打开文件 "%s" 进行读取。', file_name);
    end
    
    % --- 读取文件头信息 ---
    line1 = fgetl(fid);
    if ~ischar(line1), fclose(fid); error('文件为空或已损坏: 无法读取 line 1'); end
    header1 = sscanf(line1, '%d');
    N = header1(1);
    M = header1(2);
    
    is_nonbinary = (length(header1) == 3);
    if is_nonbinary
        Q = header1(3);
    else
        Q = 2; % 默认为二进制
    end
    
    line2 = fgetl(fid);
    if ~ischar(line2), fclose(fid); error('文件已损坏: 无法读取 line 2'); end
    header2 = sscanf(line2, '%d');
    
    line3 = fgetl(fid);
    if ~ischar(line3), fclose(fid); error('文件已损坏: 无法读取 line 3'); end
    dv = sscanf(line3, '%d');
    
    line4 = fgetl(fid);
    if ~ischar(line4), fclose(fid); error('文件已损坏: 无法读取 line 4'); end
    dc = sscanf(line4, '%d');
    
    % --- 初始化稀疏矩阵向量 ---
    rows = [];
    cols = [];
    vals = [];
    
    % --- 读取变量节点（列）列表 ---
    for col_idx = 1:N
        line = fgetl(fid);
        
        % *********** 修正开始 ***********
        % 检查我们是否在 N 次循环完成前就到达了文件末尾
        if ~ischar(line)
            warning('alist2sparse:UnexpectedEOF', ...
                '文件头声明 N=%d, 但文件在第 %d 个变量节点后提前结束。', N, col_idx-1);
            break; % 提前跳出循环
        end
        % *********** 修正结束 ***********
        
        data = sscanf(line, '%d');
        
        if is_nonbinary
            connected_rows = data(1:2:end);
        else
            connected_rows = data;
        end
        
        num_connections = length(connected_rows);
        rows = [rows; connected_rows(:)];
        cols = [cols; repmat(col_idx, num_connections, 1)];
        vals = [vals; ones(num_connections, 1)]; 
    end
    
    fclose(fid);
    
    % --- 创建稀疏矩阵 ---
    % 即使文件提前结束 (break)，我们仍然使用头信息中的 M 和 N
    % 来创建矩阵，以保持维度一致性。
    H = sparse(rows, cols, vals, M, N);
end