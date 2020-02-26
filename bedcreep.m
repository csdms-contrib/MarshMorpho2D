function z = bedcreep(z, A, crMUD, crMARSH, dx, dt, VEG)
% function to compute change in bed elevation due to ... ?
% Inputs:
    % z - bed elevation
    % A - defines cell types in domain
        % 0: not in domain
        % 1: normal cell
        % 2: open sea boundary cell
    % crMUD - unvegetated creep coefficient [m2/day]
    % cdMARSH - vegetated creep coefficient [m2/day]
    % dx - cell size [m]
    % dt - length of one time step [days]
    % VEG - water depth above the lower limit for vegetation growth
    
% Outputs:
    % z - updated bed elevation
    
    creep = A*0;
    creep(VEG==0) = crMUD;           
    creep(VEG==1) = crMARSH;

    D = (creep)/(dx^2)*dt;

    G = 0*z;
    p = find(A==1); % exclude the NOLAND CELLS (A==0)
    NN = length(p);
    G(p) = [1:NN];
    rhs = z(p);
    [N,M] = size(G);
    i = [];
    j = [];
    s = [];

    S = 0*A;
    [row, col] = ind2sub(size(A), p);
    
    for k = [N -1 1 -N] 
        [a,q] = excludeboundarycell(k,N,M,p);
        a = a(A(q(a))==1); % only include the cells in which you can creep to

        gradF = abs(z(p(a))-z(q(a)))/dx;
        facNL = gradF>=tan(1/180*pi);

        value = (D(p(a))+D(q(a)))/2.*facNL;

        S(p(a)) = S(p(a))+value; % exit from that cell
        i = [i;G(q(a))];
        j = [j;G(p(a))];
        s = [s;-value]; % gain from the neighboring cell
    end

    % summary of the material that exits the cell
    i = [i;G(p)];
    j = [j;G(p)];
    s = [s;1+S(p)];
    ds2 = sparse(i,j,s); % solve the matrix inversion
    P = ds2\rhs;
    z(G>0) = full(P(G(G>0)));

end
