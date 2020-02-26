function [N,M,dx,A,z,x,y,msl] = initializegeometrySWENEY()
% initialize geometry from SWENEY site

    dx=2;
    load SWENEY; % topography for site
    z=-zSW3(1:end-10,:); % topo variable from SWENEY
    z=double(z); % conversion from single to double type

    %attach
    z(end-25:end,132:141)=-1.35;

    [N,M]=size(z);
    A=ones(N,M);
    A(z==-999 | z<-1.9)=0;
    A(end,141:155)=2;
    A(isnan(z))=0;
    z(A==0)=-2;

    x=[0:N-1]*dx/1000;
    y=[0:M-1]*dx/1000;

    msl=0;
    
end