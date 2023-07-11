function [EmD,SSM,FLX]=sedtran(numeric,d,A,SPCLcell,DoMUD,DiffS,h,ho,E,WS,dx,dt,rbulk,co,Ux,Uy,FLX,fTide,Ttide,URx,URy,UresX,UresY,periodic,computeriver,computetide,residualcurrents,kro,MUD,coMUD,Qsmouthi);

rivermouthfront=SPCLcell.rivermouthfront;

QmouthRiver=FLX(1);
QseaTide=FLX(2);
QseaRiver=FLX(3);
QmouthTide=FLX(4);
%consider the pond cells as A==0
%A(A==3)=1; %but do not update it!
%A(A==22)=2;

A(ho<=0)=0; %ATTENTION: eliminate the cells in which the water depth is too small. They cannot evolve at all!!!! They will not erode nor accrete!!!!!!!!!!!!!!!!!!!!!!
A(rivermouthfront)=1;%the cells in front of the river mouth: let them erode if needed

%h=max(0.1,h);

%rivermouthfront
p = find(A>0);%exclude the NOLAND CELLS (A==0)
G=0*d;NN=length(p);G(p)=[1:NN];
rhs=E(p); %in the rhs there is already the additon of the erosion input
[N,M]=size(G);
i=[];j=[];s=[];S=0*G;

%boundary conditions imposed SSC
a=find(A==2);rhs(G(a))=co.*h(a).*fTide(a);

%iver and mud (if mud, you impose the SSC at the inlet)
if MUD==1;
a=find(A==10);rhs(G(a))=coMUD.*h(a).*fTide(a);
end

%Dturb=0.1*24*3600;  %QUESTO PER ORA COME REFERENZA, POI METTI ALTRO!!!
%Dturb=5.9*0.01*h*24*3600;+Dturb
%Dxx=(Dwave*24*3600.*fTide+DiffS*Ttide/2.*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
%Dyy=(Dwave*24*3600.*fTide+DiffS*Ttide/2.*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);
%Dxx=(Dwave*24*3600+DiffS*Ttide/2.*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
%Dyy=(Dwave*24*3600+DiffS*Ttide/2.*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);

% Dxx=(Dwave*h*24*3600 +DiffS*Ttide/2.*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
% Dyy=(Dwave*h*24*3600 +DiffS*Ttide/2.*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);
%diffwave*24*3600 +
Dxx=(DoMUD*24*3600 +DiffS*Ttide/2.*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
Dyy=(DoMUD*24*3600 +DiffS*Ttide/2.*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);
%Dxx=(Dwave*10*24*3600*fTide +DiffS*Ttide/2.*(abs(Ux.*Ux))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
%Dyy=(Dwave*10*24*3600*fTide +DiffS*Ttide/2.*(abs(Uy.*Uy))*(24*3600).^2)/(dx^2).*h;%.*h;%.*(ho>kro);%.*(hgross>0.01);

%Dxx=(Dwave*24*3600+DiffS*Ttide/2.*(abs(Ux.*Ux*0.1))*(24*3600).^2)/(dx^2).*h;%.*h.^2;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
%Dyy=(Dwave*24*3600+DiffS*Ttide/2.*(abs(Uy.*Uy*0.1))*(24*3600).^2)/(dx^2).*h;%.*h.^2;%.*h;%.*(ho>kro);%.*(hgross>0.01);

%Dxx=(Dwave*24*3600+DiffS*Ttide/2.*(abs(Ux.^4*8/2))*(24*3600).^2)/(dx^2).*h;%.*h.^2;%.*h;%.*(ho>kro);%.*(hgross>0.01);%% the h is not the coefficient for diffusion, is the h in the mass balance eq.
%Dyy=(Dwave*24*3600+DiffS*Ttide/2.*(abs(Uy.^4*8/2))*(24*3600).^2)/(dx^2).*h;%.*h.^2;%.*h;%.*(ho>kro);%.*(hgross>0.01);

Dxx(A==10)=0;Dyy(A==10)=0; %no tidal flux at the river mouth
Dxx(rivermouthfront)=0;Dyy(rivermouthfront)=0; %no tidal flux at the river mouth
%the factor 24*3600 is used to convert the Ux and Uy from m/s to m/day

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]  %go over the 4 directions for the gradients

%avoid to the the cells out of the domain (risk to make it periodic...)
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end

a=a(A(q(a))==1 | A(q(a))==2 | A(q(a))==10);%exlcude the translated cell that are NOLAND cells

%go over the extradiagonal direction for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;

if numeric==1
DD=(D(p(a))+D(q(a)))/2;%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
elseif numeric==0
DD=min(D(p(a)),D(q(a)));%fa meno levees on the sides. sopratutto con la sand
end

Fin=0*DD;Fin(A(p(a))==1)=1; % (A(q(a))==1 | A(q(a))==10)=1; to conserve the mass a the river mouth = no input
Fout=0*DD;Fout(A(q(a))==1)=1; %needed not to affect the b.c. -> Do not choose 2 and 1p

%tidal dispersion component
value=DD./h(p(a))./fTide(p(a));

%river flow component
if computeriver==1
if (k==N);UR=URy(p(a));up=find(UR>0);F=UR(up);end %East-west
if (k==-N);UR=URy(p(a));up=find(UR<0);F=-UR(up);end
if (k==1);UR=URx(p(a));up=find(UR>0);F=UR(up);end  %North-south
if (k==-1);UR=URx(p(a));up=find(UR<0);F=-UR(up);end
value(up)=value(up)+F*3600*24/dx;
end

if residualcurrents==1;%tidal residual currents ansd transport. 
   % (I imposed no residual currents are the opne boundary to avoid
   % calcuitating the fluxes to get the mass balance at 100%)
if (k==N);UR=UresY(p(a));up=find(UR>0);F=UR(up);end
if (k==-N);UR=UresY(p(a));up=find(UR<0);F=-UR(up);end
if (k==1);UR=UresX(p(a));up=find(UR>0);F=UR(up);end
if (k==-1);UR=UresX(p(a));up=find(UR<0);F=-UR(up);end
value(up)=value(up)+F*3600*24/dx;
end



S(p(a))=S(p(a))+value.*Fin; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-value.*Fout]; %gain from the neigborh cell

end

%summary of the material that exits the cell
settling=24*3600*WS(p)./h(p).*fTide(p);%.*(h(p)>3);

%sea boundary
a=find(A(p)==2);%find the co b.c.
settling(a)=0;%do not settle in the b.c. (to not change the SSC)
S(p(a))=1;%to impose the b.c.

%river boundary
if MUD==1; %if mud, handlew this as an imposed SSC
a=find(A(p)==10);%find the co b.c.
settling(a)=0;%do not settle in the b.c. (to not change the SSC)
S(p(a))=1;%to impose the b.c.
end

i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)+settling];
ds2 = sparse(i,j,s);%solve the matrix inversion
P=ds2\rhs;%P=mldivide(ds2,rhs);

SSM=0*d;SSM(G>0)=full(P(G(G>0)));%rescale the matrix

%update the bed
EmD=0*A;EmD(p)=(E(p)-SSM(p).*settling)/rbulk;



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%OUTPUT FOR SSM BALANCE. DOES NOT AFFECT COMPUTATION!!!
p = find(A==1);[row col]=ind2sub(size(A),p);

%Sea boundary%%%%%%%%%%%%%%%
Q=0;QR=0;
for k = [N -1 1 -N]
   
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end

%only get the cell that discharge into a b.c. A==2
a=a(A(q(a))==2);

%for each gradient direction
if (k==N | k==-N);D=Dyy;else;D=Dxx;end;
if numeric==1
DD=(D(p(a))+D(q(a)))/2;%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
elseif numeric==0
DD=min(D(p(a)),D(q(a)));%fa meno levees on the sides. sopratutto con la sand
end

%Tide
Q=Q+sum(DD.*(SSM(p(a))./h(p(a))./fTide(p(a))-SSM(q(a))./h(q(a))./fTide(q(a)))); %exit from that cell
%River
if (k==N | k==-N);QR=QR+sum((SSM(p(a)).*URy(p(a))))*sign(k);end
if (k==1 | k==-1);QR=QR+sum((SSM(p(a)).*URx(p(a))))*sign(k);end
end

QseaTide=QseaTide+dx*dt*Q/rbulk; %(note the the divided dx is for the gradient, not for the cell width!)
QseaRiver=QseaRiver+dt*3600*24*QR/rbulk;
%%%%%%%%%%%%%%%%%%%



%River mouth%%%%%%%%%%%%%%%%
Q=0;
for k = [N -1 1 -N]
if periodic==0
[a,q]=excludeboundarycell(k,N,M,p);
elseif periodic==1;
[a,q]=periodicY(k,N,M,p); %for the long-shore
end
%exlcude the translated cell that are the b.c.
a=a(A(q(a))==10);
if (k==N | k==-N);D=Dyy;else;D=Dxx;end
if numeric==1%mud
DD=(D(p(a))+D(q(a)))/2;%THE ORIGINAL FOR MARSH fa piu' accumulo sui lati dei canali, cioe canali piou stretti
elseif numeric==0%sand
DD=min(D(p(a)),D(q(a)));%fa meno levees on the sides. sopratutto con la sand
end

%Tide
Q=Q+sum(DD.*(SSM(p(a))./h(p(a))./fTide(p(a))-SSM(q(a))./h(q(a))./fTide(q(a)))); %exit from that cell
%River
if (k==N | k==-N);QR=QR+sum((SSM(p(a)).*URy(p(a))))*sign(k);end
if (k==1 | k==-1);QR=QR+sum((SSM(p(a)).*URx(p(a))))*sign(k);end
end

QmouthTide=QmouthTide+dx*dt*Q/rbulk;
QmouthRiver=QmouthRiver+dt*3600*24*QR/rbulk;
%%%%%%%%%%%%%%%%%%%



%SUM OF THE RIVER INPUT (IMPOSED FROM THE OUTISE USING Qsmouthi, which is cl
QmouthRiver=FLX(1)+Qsmouthi*dt*24*3600/rbulk;
FLX=[QmouthRiver;QseaTide;QseaRiver;QmouthTide];





