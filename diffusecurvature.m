function F=diffusecurvature(A,F,alpha,dx,Ux,Uy);

%Uo=U;

U=sqrt(Ux.^2+Uy.^2);
%Ux=max(0.1,abs(Ux)./U);
%Uy=max(0.1,abs(Uy)./U);
Ux=abs(Ux)./U;
Uy=abs(Uy)./U;

%%%only perdendicular to flow
Dx=A*(alpha*3600*24)/(dx^2).*Ux;
Dy=A*(alpha*3600*24)/(dx^2).*Uy;

% %all directions
% %D=A*(alpha*3600*24)/(dx^2);
% Dx=A*(alpha*3600*24)/(dx^2);
% Dy=A*(alpha*3600*24)/(dx^2);



% Dx=A*(alpha*3600*24)/(dx^2).*0.5;
% Dy=A*(alpha*3600*24)/(dx^2).*0.5;

G=0*F;
p=find(A==1);%exclude the NOLAND CELLS (A==0)
NN=length(p);G(p)=[1:NN];rhs=F(p);[m,n]=size(G);i=[];j=[];s=[];

S=0*G;
[row col]=ind2sub(size(A),p);
for k = [m -1 1 -m]
%avoid to the the cells out of the domain (risk to make it periodic...)
if k==m;aa=find(col+1<=n);end;if k==-m;aa=find(col-1>0);end;if k==-1;aa=find(row-1>0);end;if k==1;aa=find(row+1<=m);end;

q=p+k;%the translated cell
a=aa(A(q(aa))==1 );%only inclued the cells in whcih you can creep to

if abs(k)==1;value=(Dx(p(a))+Dx(q(a)))/2;end
if abs(k)==m;value=(Dy(p(a))+Dy(q(a)))/2;end
%value=(D(p(a))+D(q(a)))/2;

S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;F(G>0)=full(P(G(G>0)));


