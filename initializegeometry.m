function [N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell]=initializegeometry(P);

% %geometric parameters
% dx=50*2;
% N=200/2; %x
% M=200/2; %y

% 
% dx=50/4;
% N=(200-70)*4; %x
% M=100*4  *0.6; %y

% dx=50/5*2;
% N=(200-86)*2*5/2; %x
% M=100*0.5*2*5/2; %y
dx=5*2;
N=400;%/2; %x
M=200;%/4; %y

%%%%%%%%%%%cell types
A=ones(N,M);
%sea b.c.
A(end,:)=2;
%A(end-1:end,:)=2;
%river b.c.
rivermouthfront=[];

% mouthW=1; %2
% A(1:2,1:M/2-1-mouthW)=0;%create a concreate wall on the sides
% A(1:2,M/2+1+mouthW:end)=0;%create a concreate wall on the sides
% A(1,M/2-mouthW:M/2+mouthW)=10;
% %these are the cells in front of the river mouth
% S=A*0;S(2,M/2-mouthW:M/2+mouthW)=1;rivermouthfront=find(S==1);clear S;

SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;

%bathymetry
%initial profile. ELEVATION AT MSL
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


%slope=0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-3;%sloping   0.5
%slope=0.2/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-1+1;%sloping   0.5
slope=1/1*0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M);%-0.5;%sloping   0.5
%z=z*0+5;
%slope=0.5/1000;z=[0:N-1]*slope*dx;z=z'*ones(1,M)-3;%sloping
%zin=5;z=zin*ones(N,M);%flat

% slopebreak=30;
% z=[1:N];
% slope=1/1000;z(1:N-slopebreak)=[1:N-slopebreak]*slope*dx;
% slope=1/500;z(slopebreak+1:N)=+z(slopebreak)+[slopebreak+1:N]*slope*dx;
% z=z'*ones(1,M);%sloping

% %barriers
%pos=1;bwidth=1;%250/dx;
%A(end-pos:end-pos+bwidth,1:M/2-15)=0;
%A(end-pos:end-pos+bwidth,M/2+15:end)=0;
%z(end-pos:end-pos+bwidth,1:M/2-1)=-5.1;
%z(end-pos:end-pos+bwidth,M/2+1:end)=-5.1;



%z(z>1)=1;
%channel
z(end-100:end,M/2-2/2:M/2+2/2)=5;
%z(:,M/2-2:M/2+2)=(ones(1,5)'*[0:N-1]/N*10)';
%z(1:150,1:2)=-2;

%THIS IS THE GOOD ONE
%z(end-pos-bwidth:end-pos,:)=4;%3;
%z(end-pos+1:end,:)=z(end-pos+1:end,:)+10;
%z(end-pos:end-pos+bwidth,1:M/2-4)=-2.1;
%z(end-pos:end-pos+bwidth,M/2+4:end)=-2.1;

%z(end-20:end,M/2-2:M/2+2)=3;

%random
z=z+2*(rand(N,M)-0.5)*0.2;% *0.5   *0.2;

%river
%z(A==10)=P.hmouth;
%z(1:40,M/2-mouthW:M/2+mouthW)=P.hmouth;
%load STORE
%z=-squeeze(STORE(:,:,end-60));

% [N,M]=size(z);
% dx=20;
% z=interp2([0:N-1],[0:M-1]',z',[0:2:N-1],[0:2:M-1]','nearest')';
% A=z*0+1;
% A(end,:)=2;
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


% [N,M]=size(z);
% dx=5;
% z=interp2([0:N-1],[0:M-1]',z',[0:0.5:N-1],[0:0.5:M-1]','nearest')';
% A=z*0+1;
% A(end,:)=2;
% [N,M]=size(z);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;

% % 
%load STORE
%z=-STORE(:,:,end);

z(A==0)=NaN;
%%%%%%%%%%%%%%%%%%


msl=0;
zbedo=-z-msl;
Active=zbedo<P.Trange/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy_3sediments(z,N,M,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%wave boundary conditions
%the lateral boundaries for waves
%postivie is left boundary; negative is right boudndary; 
%1 is no-gradient in wave. THIS IS ALSO no-gradient in wave-indcued lateral
%sediment transport
%2 is zero wave height . Also implies no sand transport at inlet
AW=A*0;%the internal points
%right boundary
AW(:,end)=-1;
%left boundary
AW(:,1)=1;
%AW(1:50,1)=1;
%%%%%%%%%%%%%%%%%%%%%
