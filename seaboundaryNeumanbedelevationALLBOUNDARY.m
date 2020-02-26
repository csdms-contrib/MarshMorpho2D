function [z] = seaboundaryNeumanbedelevationALLBOUNDARY(A, z)
% get the bed elevation at the sea boundary cells (Neumann BC)
% Input:
    % A: Cell types (0 = not in domain, 1 = normal cell, 2 = boundary cell)
    % z: bed elevation
% Output:
    % z: bed elevation

    [N, M] = size(A);
    p = find(A==2); % exclude the NOLAND CELLS (A==0)

    [row, col] = ind2sub(size(A), p);
    for k = [N -1 1 -N] 

        [a, q] = excludeboundarycell(k, N, M, p);
        a = a( A(q(a))==1 ); % only included the cells in whcih you can creep to

        z(p(a)) = z(q(a));

    end
    
end