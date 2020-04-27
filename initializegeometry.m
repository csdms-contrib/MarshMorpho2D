function [N,M,dx,A,z,Active,x,y,msl]=initializegeometry(P);

% %geometric parameters

dx=5*2;
N=1000/2*1;%x
M=600/2*1;%y

%%%%%%%%%%%cell types
A=ones(N,M);

%sea b.c.
A(end,:)=2;

%bathymetry
%initial profile. ELEVATION AT MSL
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


slope=1/1*0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1;%sloping

%channel in the middle
%z(end-100:end,M/2-5:M/2+5)=50;


%random
z=z+2*(rand(N,M)-0.5)*0.2;% *0.5   *0.2;

z(end-1:end,:)=10;

z=-z;

z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%




msl=0;
zbedo=-z-msl;
Active=zbedo<P.Trange/2;


