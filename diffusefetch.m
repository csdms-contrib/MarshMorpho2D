function F=diffusefetch(A,F,alpha,dx);


D=A*(alpha*3600*24)/(dx^2);


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

value=min(D(p(a)),D(q(a)));%value=(D(p(a))+D(q(a)))/2.*facNL;

S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;F(G>0)=full(P(G(G>0)));


