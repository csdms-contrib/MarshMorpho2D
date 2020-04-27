function [S,deltaY2,pondlost]=isolatedpondexpansion(z,S,A,N,M,dx,dt,zpondcr,maxdpond,aPEXP,pondlost);

zoriginal=z;

p=find(S==1);%find the existing ponds

[row col]=ind2sub([N M],p);

for k = [N -1 1 -N]

%avoid to the the cells out of the domain (risk to make it periodic...)
if k==N;a=find(col+1<=M);end;
if k==-N;a=find(col-1>0);end;
if k==-1;a=find(row-1>0);end;
if k==1;a=find(row+1<=N);end;

q=p+k;%the translated cell
er=a(A(q(a))==1 & S(q(a))==0); %find the neightr cell that is a non-pond cell


%rng('shuffle');
r=rand(length(er),1);% you only need rng at the beginnign of the loop
aa=find(r<aPEXP*dt/365/dx);

%S(q(er(aa)))=1;


%you go over all pond cells (p) and you allow them to expand along all (up to
%four) edges with a marsh cell (q). So the cell that eventally erodes is q (not p)
for i=1:length(aa)
%dz=max(-z(q(er(aa(i))))-msl-zpondcr,0);%erode and set elevation fixed at zpondcr
%dz=max(-z(q(er(aa(i))))+z(p(er(aa(i)))),0);%erode and set elevation equal to next pond cell
dz=max(z(q(er(aa(i))))-z(p(er(aa(i)))),0);%erode and set elevation equal to next pond cell
pondbasemin=max(0,z(q(er(aa(i))))-zpondcr);
dz=min(dz,pondbasemin);
dz=min(maxdpond,dz);%maximum scour depth

% dz=max(-z(p)-msl-zpondcr,0);%base level of the pond
%dzvicino=max(-z(q(er(aa(i))))+z(p(er(aa(i)))),0);
%dz=min(dz,dzvicino);
% dz=min(maxdpond,dz);%maximum scour depth

    if dz>0
    z(q(er(aa(i))))=z(q(er(aa(i))))-dz;
    S(q(er(aa(i))))=1;%also update the pond ID
    pondlost=pondlost+dz;
    end
end

end

deltaY=zoriginal-z;
deltaY2=deltaY;

% figure
% imagesc(deltaY2)
% pause



% function [A,d]=isolatedpondexpansion(A,d,m,n,dx,dt,MF);
% 
% p = find(A==1);[row col]=ind2sub(size(A),p);
% 
% for k = [m -1 1 -m]
% 
% %avoid to the the cells out of the domain (risk to make it periodic...)
% if k==m;aa=find(col+1<=n);end;
% if k==-m;aa=find(col-1>0);end;
% if k==-1;aa=find(row-1>0);end;
% if k==1;aa=find(row+1<=m);end;
% 
% q=p+k;%the translated cell
% er=aa(A(q(aa))==3);
% 
% rng('shuffle');r=rand(length(er),1);% you only need rng at the beginnign of the loop
% a=find(r<0.1/365*dt*MF/dx);
% 
% A(p(er(a)))=3;
% end