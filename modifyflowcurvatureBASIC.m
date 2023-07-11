function F=diffusefetch(A,F,alpha,dx,URx,URy,MASK);

Fo=F;
[N,M]=size(A);

G=0*F;
p=find(A==1);%exclude the NOLAND CELLS (A==0)
NN=length(p);G(p)=[1:NN];rhs=F(p);[m,n]=size(G);i=[];j=[];s=[];

store=0*G;

S=0*G;
[row col]=ind2sub(size(A),p);
for k = [m -1 1 -m]
%avoid to the the cells out of the domain (risk to make it periodic...)


% 
% if k==m;aa=find(col+1<=n);end;if k==-m;aa=find(col-1>0);end;if k==-1;aa=find(row-1>0);end;if k==1;aa=find(row+1<=m);end;
% q=p+k;%the translated cell
% a=aa(A(q(aa))==1 );%only inclued the cells in whcih you can creep to

[a,q,q2]=excludeboundarycell_opposite(k,N,M,p);
a=a(A(q(a))==1 );%only inclued the cells in whcih you can creep to



value=0*A(p(a));

if (k==N);UR=URy(p(a));up=find(UR>0);Y=UR(up);end %East-west
if (k==-N);UR=URy(p(a));up=find(UR<0);Y=-UR(up);end
if (k==1);UR=URx(p(a));up=find(UR>0);Y=UR(up);end  %North-south
if (k==-1);UR=URx(p(a));up=find(UR<0);Y=-UR(up);end

% if (k==N);UR=(URy(p(a))+URy(q(a)))/2;up=find(UR>0);Y=UR(up);end %East-west
% if (k==-N);UR=(URy(p(a))+URy(q(a)))/2;up=find(UR<0);Y=-UR(up);end
% if (k==1);UR=(URx(p(a))+URx(q(a)))/2;up=find(UR>0);Y=UR(up);end  %North-south
% if (k==-1);UR=(URx(p(a))+URx(q(a)))/2;up=find(UR<0);Y=-UR(up);end

value(up)=value(up)+alpha*Y/dx;
%value(up)=value(up).*(1+ double(MASK(q(a(up)))==0));
%value( up (MASK(q(a(up)))==0) )=value( up (MASK(q(a(up)))==0) )*10/dx;

%store(p(a( up (MASK(q(a(up)))==0) )))=1;



% q2(q2<1)=p(q2<1);
% q2(q2>N*M)=p(q2>N*M);
% leak=(MASK(p(a))==1 & (MASK(q2(a))==1 | MASK(q(a))==1));
% %leak=(MASK(p(a))==1 & (MASK(q2(a))==1));

    
S(p(a))=S(p(a))+value;%.*leak; %exit from that cell .*leak
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell - maybe this is what it gives to the neiburing cell


end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;F(G>0)=full(P(G(G>0)));


F(store==1)=max(F(store==1),Fo(store==1));
