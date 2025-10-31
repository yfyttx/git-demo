function matrix2alist_nb(input_filename, output_filename, GF_size)
% matrix2alist_nb: Converts a non-binary parity matrix to the .alist format.
%
% FINAL CORRECTED VERSION - Fixes the logic for generating the row-perspective list.
% This version only writes positions to be compatible with binary readers.

    fprintf('Loading matrix from %s...\n', input_filename);
    H = load(input_filename);

    % Use logical matrix for structure (1 for non-zero, 0 for zero)
    H_bin = (H > 0); 

    [M, N] = size(H);
    dc = sum(H_bin, 2); % Row weights (degrees of check nodes)
    dv = sum(H_bin, 1); % Column weights (degrees of variable nodes)
    dcmax = max(dc);
    dvmax = max(dv);

    fprintf('Writing to %s...\n', output_filename);
    fileID = fopen(output_filename, 'w');

    % --- Header ---
    fprintf(fileID, '%d %d\n', N, M); % Standard .alist often omits GF_size here
    fprintf(fileID, '%d %d\n', dvmax, dcmax);

    % --- Degree Information ---
    fprintf(fileID, '%d ', dv);
    fprintf(fileID, '\n');
    fprintf(fileID, '%d ', dc);
    fprintf(fileID, '\n');

    % --- Column-Perspective List (Correctly generated) ---
    for i = 1:N % For each column (variable node)
        % Find the indices of non-zero elements in this column
        row_indices = find(H_bin(:, i));
        fprintf(fileID, '%d ', row_indices);
        fprintf(fileID, '\n');
    end

    % =======================================================================
    % --- Row-Perspective List (THIS IS THE CORRECTED LOGIC) ---
    % =======================================================================
    for i = 1:M % For each row (check node)
        % Find the indices of non-zero elements in this specific row
        col_indices = find(H_bin(i, :));
        fprintf(fileID, '%d ', col_indices);
        fprintf(fileID, '\n');
    end

    fclose(fileID);
    fprintf('Successfully created %s with correct logic.\n', output_filename);

end