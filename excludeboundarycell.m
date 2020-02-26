function [a, q] = excludeboundarycell(k, N, M, p)
% function to exclude boundary cells from domain
% Input:
    % k: control where to assign a
    % N: number of rows in model domain 
    % M: number of columns in model domain
    % p: specific cell indicies (boundary cells?)
% Output:
    % a: translated boundary cell grid
    % q: translated boundary cell index?
    
    [row, col] = ind2sub([N, M], p);

    if k == N
        a = find(col+1 <= M);
    end

    if k == -N
        a = find(col-1 > 0);
    end

    if k == -1
        a = find(row-1 > 0);
    end

    if k == 1
        a = find(row+1 <= N);
    end

    q = p + k; % the translated cell

end