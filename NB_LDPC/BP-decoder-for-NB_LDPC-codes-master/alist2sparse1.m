function [H] = alist2sparse(fname, q) 
% 非二进制版alist2sparse：读取非二进制LDPC码的alist文件，返回非二进制稀疏矩阵
% 输入：fname - alist文件名；q - 域大小（如GF(2^6)=64）
% 输出：H - 非二进制稀疏矩阵（元素为GF(q)中的非零值）

% 原版权声明保留
%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
%   2025 非二进制适配修改：保留非零元素数值，支持GF(q)

fid = fopen(fname); 
if fid == -1
    error('无法打开文件：%s，请检查路径是否正确', fname);
end

% 读取alist文件头部信息（非二进制alist格式与二进制一致，仅后续元素值不同）
n = fscanf(fid,'%d',1);       % 码长（符号数N）
m = fscanf(fid,'%d',1);       % 校验节点数（M）
maxinrow = fscanf(fid,'%d',1);% 每行最大非零元素数
junk = fscanf(fid,'%d',1);    % 无用参数（二进制alist中的maxincol，非二进制忽略）

% 读取每行的非零元素个数（num：校验节点行的非零数；num_col：符号节点列的非零数，无用）
num = fscanf(fid,'%d',[1 m]); % 第1~m行（校验节点）的非零元素个数
num_col = fscanf(fid,'%d',[1 n]); % 第1~n列（符号节点）的非零元素个数，无用

% 读取非零元素的【列索引】和【对应数值】（非二进制alist的核心：每个列索引后紧跟元素值）
position = zeros(m, maxinrow); % 存储列索引
values = zeros(m, maxinrow);    % 存储非二进制元素值（新增）
for i = 1:m  % 遍历每个校验节点行
    for j = 1:maxinrow     
        % 非二进制alist格式：列索引 + 元素值，成对出现
        pos = fscanf(fid,'%d',1);  % 列索引（1-based）
        val = fscanf(fid,'%d',1);  % 元素值（GF(q)中的非零值，1~q-1）
        position(i,j) = pos;
        values(i,j) = val;
    end
end 

% 构建稀疏矩阵的行、列、值索引
ii = zeros(1, sum(num));  % 行索引（校验节点）
jj = zeros(1, sum(num));  % 列索引（符号节点）
vv = zeros(1, sum(num));  % 非零元素值（新增）
k = 1; 
for i = 1:m  % 遍历每个校验节点行
    for j = 1:num(i)  % 仅遍历实际有非零元素的列
        ii(k) = i;                % 行索引（校验节点i）
        jj(k) = position(i,j);    % 列索引（符号节点）
        vv(k) = values(i,j);      % 非二进制元素值（保留原值，不设为1）
        k = k + 1;  
    end
end 

% 生成非二进制稀疏矩阵（m行校验节点，n列符号节点）
H = sparse(ii, jj, vv, m, n); 
fclose(fid);

% 验证：打印矩阵基本信息
fprintf('成功读取非二进制alist文件：%s\n', fname);
fprintf('矩阵维度：%d行（校验节点）×%d列（符号节点），域大小GF(%d)\n', m, n, q);
fprintf('非零元素个数：%d\n', nnz(H));
end