function z=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,crbank,dx,dt,VEG,S,Qs,rbulk,alphaMUD,facQsbank,U);%,deltaUC,a_bankcreep);


A(Active==0)=0;

%%%%%%%!!!!!!!$$$$$$
%%%%%%%!!!!!!!$$$$$$
%A(A==2)=1;%%trucco per fare creep also at the boundary


%%downslope dependnets on sediment transport
Qs=Qs/rbulk;

creep=A*0;
creep(VEG==0)=crMUD;%+(alphaMUD*3600*24*Qs(VEG==0));          
creep(VEG==1)=crMARSH;    

creepQs=A*0;
creepQs(VEG==0)=(alphaMUD*3600*24*Qs(VEG==0));

%creepUc=A*0;
%creepUc(VEG==1)=(0.0001*3600*24*deltaUC(VEG==1));

% %%%%%%%!!!!!!!$$$$$$
% %trucco per tenere il seaward basso
% creep(end-100:end,:)=creep(end-100:end,:)+100;%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=(creep)/(dx^2)*dt;
DQs=creepQs/(dx^2)*dt;
Dcrbank=crbank/(dx^2)*dt;

Dcrbank=Dcrbank*dx;%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%$%##^Z$$
facQsbank=facQsbank/10*dx;

%DUc=a_bankcreep*0.1*deltaUC/(dx^2)*dt;

%consider the pond cells as A==0
%A(S==1)=0; %DO NOT CREEP INTO PONDS!!!! NOT USED ANYMOREEEE!!!

G=0*z;
p=find(A==1);%exclude the NOLAND CELLS(A==0)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NN=length(p);G(p)=[1:NN];rhs=z(p);[N,M]=size(G);i=[];j=[];s=[];

Spond=S;%%%%ATTENZIONE
S=0*G; %This S ia a different emanign that ponds. it is just to store values

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N] 
[a,q]=excludeboundarycell(k,N,M,p);
a=a(A(q(a))==1);%only inclued the cells in whcih you can creep!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

value=(D(p(a))+D(q(a)))/2 +(DQs(p(a))+DQs(q(a)))/2;%%  +max(DUc(p(a)),DUc(q(a)));%  +(DUc(p(a))+DUc(q(a)))/2.*max(VEG(p(a))==1,VEG(q(a))==1); 
%value=max(D(p(a)),D(q(a))) +(DQs(p(a))+DQs(q(a)))/2.*min(VEG(p(a))==0,VEG(q(a))==0);%  +(DUc(p(a))+DUc(q(a)))/2.*max(VEG(p(a))==1,VEG(q(a))==1); 
% 
% valQS=(DQs(p(a))+DQs(q(a)))/2;
% value( (VEG(p(a))==1 & VEG(q(a))==0) | (VEG(p(a))==0 & VEG(q(a))==1) )= Dcrbank+...
%     valQS( (VEG(p(a))==1 & VEG(q(a))==0) | (VEG(p(a))==0 & VEG(q(a))==1) )*facQsbank;


valQS=facQsbank*max(U(p(a)),U(q(a)))/(dx^2)*dt;
value( (VEG(p(a))==1 & VEG(q(a))==0) | (VEG(p(a))==0 & VEG(q(a))==1) )= Dcrbank+...
    valQS( (VEG(p(a))==1 & VEG(q(a))==0) | (VEG(p(a))==0 & VEG(q(a))==1) );

%value=(D(p(a))+D(q(a)))/2 +(DQs(p(a))+DQs(q(a)))/2;%  +(DUc(p(a))+DUc(q(a)))/2.*max(VEG(p(a))==1,VEG(q(a))==1); 
value=value.*min(Yreduction(p(a)),Yreduction(q(a)));


% gradF=abs(z(p(a))-z(q(a)))/dx;
% facNL=gradF>=tan(10/180*pi);
% value=value.*facNL;


%The standrd
%DO NOT CREEP AT THE POND EDGE
value(Spond(p(a))==1 & Spond(q(a))==0) = 0;
value(Spond(p(a))==0 & Spond(q(a))==1) = 0;

% %trucco per evitare the i canali si chiudana
% %DO NOT CREEP AT THE POND EDGE
% %"it is a pond but it would be a channel, i.e., it is not a deep channel" %cit. Giulio M
% value((Spond(p(a))==1 & VEG(p(a))==1) & Spond(q(a))==0) = 0;
% value(Spond(p(a))==0 & (Spond(q(a))==1 & VEG(q(a))==1)) = 0;


S(p(a))=S(p(a))+value; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;1+S(p)];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;z(G>0)=full(P(G(G>0)));


