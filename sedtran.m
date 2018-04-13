function [EmD,SSM,FLX]=sedtran(d,A,DiffS,h,ho,E,WS,dx,dt,rbulk,co,Ux,Uy,FLX,fTide,Ttide,kro);

A(ho<=0)=0; %eliminate the cells in which the water depth is too small

%rivermouthfront
p = find(A>0);%exclude the NOLAND CELLS (A==0)
G=0*d;NN=length(p);G(p)=[1:NN];
rhs=E(p); %in the rhs there is already the additon of the erosion input
[N,M]=size(G);
i=[];j=[];s=[];S=0*G;

%boundary conditions imposed SSC
a=find(A==2);rhs(G(a))=co.*h(a).*fTide(a);

Dxx=(DiffS*Ttide/2*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, it the h in the mass balance eq.
Dyy=(DiffS*Ttide/2*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*(ho>kro);%.*(hgross>0.01);

%the factor 24*3600 is used to convert the Ux and Uy from m/s to m/day

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]  %go over the 4 directions for the gradients

%avoid to the the cells out of the domain (risk to make it periodic...)
[a,q]=excludeboundarycell(k,N,M,p);

a=a(A(q(a))==1 | A(q(a))==2 | A(q(a))==10);%exlcude the translated cell that are NOLAND cells

%go over the extradiagonal direction for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;

DD=(D(p(a))+D(q(a)))/2;

Fin=0*DD;Fin(A(p(a))==1)=1; % (A(q(a))==1 | A(q(a))==10)=1; to conserve the mass a the river mouth = no input
Fout=0*DD;Fout(A(q(a))==1)=1; %needed not to affect the b.c. -> Do not choose 2 and 1p

%tidal dispersion component
value=DD./h(p(a))./fTide(p(a));

S(p(a))=S(p(a))+value.*Fin; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value.*Fout]; %gain from the neigborh cell

end

%summary of the material that exits the cell
settling=24*3600*WS(p)./h(p).*fTide(p);%.*(h(p)>3);

%sea boundary
a=find(A(p)==2);%find the co b.c.
settling(a)=0;%do not settle in the b.c. (to not change the SSC)
S(p(a))=1;%to impose the b.c.


i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)+settling];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;%P=mldivide(ds2,rhs);

SSM=0*d;SSM(G>0)=full(P(G(G>0)));%rescale the matrix

%update the bed
EmD=0*A;EmD(p)=(E(p)-SSM(p).*settling)/rbulk;









%%%%%%%%OUTPUT FOR SSM BALANCE. DOES NOT AFFECT COMPUTATION!!!
p = find(A==1);[row col]=ind2sub(size(A),p);

%sea boundary
Q=0;
for k = [N -1 1 -N]
   
[a,q]=excludeboundarycell(k,N,M,p);

%only get the cell that discharge into a b.c. A==2
a=a(A(q(a))==2);

%for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;
DD=(D(p(a))+D(q(a)))/2;

Q=Q+sum(DD.*(SSM(p(a))./h(p(a))./fTide(p(a))-SSM(q(a))./h(q(a))./fTide(q(a)))); %exit from that cell

end

FLX=FLX+dx*dt*Q/rbulk; %(note the the divided dx is for the gradient, not for the cell width!)
%%%%%%%%%%%%%%%%%%%



