function [N,M,dx,A,z,x,y,msl]=initializegeometry(P);

load WEST;
z=-PIE2mBELLO;
z(247:249,361:363)=-1.6;%tah house
z=z(:,1:end-24);


dx=2;
[N,M]=size(z);

%%%%%%%%%%%cell types
A=ones(N,M);

A(z==-999)=0;
A(:,end)=0;
A(95*4/dx:114*4/dx,end)=2;

A(isnan(z))=0;

z(A==0)=-2;


x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


msl=0;