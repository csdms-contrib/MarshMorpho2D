function [deltaY2,pondloss]=isolatedponddeepening(z,S,ponddeeprate,dt,pondloss,dBlo);

zoriginal=z;

%p=find(S==1);%find the existing ponds
p=find(S==1 & z>dBlo);%find the existing ponds

for i=1:length(p)
  dz=ponddeeprate*dt/365;
  z(p(i))=z(p(i))-dz;
  pondloss=pondloss+dz;  
% %dz=max(-z(q(er(aa(i))))-msl-zpondcr,0);
% dz=max(-z(p(i))+z(p(er(aa(i)))),0);
% 
% % dz=max(-z(p)-msl-zpondcr,0);%base level of the pond
% %dzvicino=max(-z(q(er(aa(i))))+z(p(er(aa(i)))),0);
% %dz=min(dz,dzvicino);
% % dz=min(maxdpond,dz);%maximum scour depth
% 
%     if dz>0
%     z(q(er(aa(i))))=z(q(er(aa(i))))+dz;
%     S(q(er(aa(i))))=1;%also update the pond ID
%     pondlost=pondlost-dz;
%     end
% end

end

deltaY=zoriginal-z;
deltaY2=deltaY;





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