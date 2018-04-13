function [U,Ux,Uy,q,P]=flowBasin(A,MANN,h,d,dx,DH,T,kro);

Uo=1;
MANN(isnan(MANN))=0.1;%
csi=h.^(1/3)./MANN.^2./Uo*24*3600;
D=csi.*h.^2/(dx^2);

G=0*d;a=find(A~=0);NN=length(a);G(a)=[1:NN];
rhs=ones(NN,1).*DH(a)/(T/2*3600*24); %in m/s!!!

[N,M] = size(G);i=[];j=[];s=[];

%boundary conditions imposed water level
a=find(A==2);
i=[i;G(a)]; j=[j;G(a)]; s=[s;ones(size(a))];rhs(G(a))=0;%water level zero

S=0*G;
%exclude the NOLAND CELLS (A==0)
p = find(A==1 | A==10);[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]

%the translated cells
[a,q]=excludeboundarycell(k,N,M,p);

a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells

DD=(D(p(a))+D(q(a)))/2;%.*(fM(p(a))+fM(q(a)))/2; %THA BEST!!!! BESTA! WITH THIS MORE STABLE

S(p(a))=S(p(a))+DD; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-DD]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)];

ds2 = sparse(i,j,s);%solve the matrix inversion
p=ds2\rhs;
P=G;P(G>0)=full(p(G(G>0)));
P(A==2)=0;  %need when swtinching q and p


D=D./h*dx;
Ux=0*A;Uy=0*A;
U1=0*A;Um1=0*A;UN=0*A;UmN=0*A;
p = find(A==1 | A==10 | A==2);[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]
%the translated cell
[a,q]=excludeboundarycell(k,N,M,p);


a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells
DD=min(D(p(a)),D(q(a)));%.*(fM(p(a))+fM(q(a)))/2; MEGLIO

if (k==1); U1(p(a))=U1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;
elseif (k==-1); Um1(p(a))=Um1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;
elseif (k==N); UN(p(a))=UN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;
elseif (k==-N); UmN(p(a))=UmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;
end

end

Uy=max(abs(U1),abs(Um1)).*sign(U1+Um1);
Ux=max(abs(UN),abs(UmN)).*sign(UN+UmN);

U=sqrt(Ux.^2+Uy.^2);
q=U.*h;








