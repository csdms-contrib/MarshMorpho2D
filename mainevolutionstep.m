function [IO,fIO,maxdeltaz,maxup,PLT,MMM]=mainevolutionstep(A,AW,SPCLcell,PARAMS,dx,dt,zb,IO,fIO,Hsurge,angleSWELL,angleWIND,t)


%Load the parameters and variables
names = fieldnames(IO);
for i=1:length(names);eval([names{i} '=IO.' names{i} ';' ]);end

names = fieldnames(fIO);
for i=1:length(names);eval([names{i} '=fIO.' names{i} ';' ]);end

names = fieldnames(PARAMS);
for i=1:length(names);eval([names{i} '=PARAMS.' names{i} ';' ]);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,M]=size(A);
optionBC=2; %impose water depth (for sea boundary)

Y=Y1+Y2+Y3;%total thickness of active later
zs=-zb+(Yb+plyr.*(dlyr)); %the heigth of the ground just below Y
z=zs+Y;%the absoulte bed elevation, positive is the bed is high ground, negative if it is low ground

zoriginal=z;

%Volumetric fraction
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
if reducefractionsediment==1;
fracY1=max(0,Y1)./Ytot;fracY1(Ytot==0)=0;%volumetric fraction of sand in active layer
fracY2=max(0,Y2)./Ytot;fracY2(Ytot==0)=0;%volumetric fraction of mud in active layer
fracY3=max(0,Y3)./Ytot;fracY3(Ytot==0)=0;%volumetric fraction of organic in active layer
else
fracY1=A*0+1;fracY2=A*0+1;fracY3=A*0+1;
end
    



%Hydrodynamic%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sea level
msl=msl+RSLR*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the upland cells active
Active((z-msl)<Trange/2)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Pond dynamics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Need to run it first because it affects the vegetation through S !!!!!
if calculateponddynamics==1     %%CHECK THE SIGN OF z!!!
zsill=NaN;%TrangeVEG/2*NaN;  
%use lev instead of z to be rleative to MSL!!!
[deltaY2,pondloss]=pondformation(A,dx,dt,Epondform,z-msl,zpondcr,maxdpond,zsill,pondloss,Active);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
[S,AC,DIF]=findisolatedponds(z-msl,A,N,M,dx,zntwrk,zsill,distdr,minponddepth,Active);%AC is only used for plotting, not in the computation   
%DIF is the amount of water impunded at low water. It is the remainng water depth in the pond!
%S(lev<dBlo)=0;% the ponds in the mudflats are not really a pond! You need to put it otherwise probelm with bedcreeppond. YOU CANNOT CREEP IN PONDS
[S,deltaY2,pondloss]=isolatedpondexpansion(z-msl,S,A,N,M,dx,dt,zpondcr,maxdpond,aPEXP,pondloss,Active);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
[deltaY2,pondloss]=isolatedponddeepening(z-msl,S,ponddeeprate,dt,pondloss,dBlo);Y2=Y2-deltaY2;z=zs+(Y1+Y2+Y3);
else;AC=A*0;S=A*0;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Accrete in the river mouth to maintain same water depth
Y1(A==10)=Y1(A==10)+RSLR*dt; %sediment on the bed of the river mouth!!
%z(A==10)=zs(A==10)-(Y1(A==10)+Y2(A==10));  %accrete the mouth ourlet with sand
z(A==10)=zs(A==10)+(Y1(A==10)-Y2(A==10));  %accrete the mouth ourlet with sand

%Water depth
if computeriver==1 & riverwaterlevel==1
hpRIV=max(0,z+5);%+5; %at the beginning, just a small water level from the river to avoid extra slopes
else;hpRIV=0;end
[h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV,tempdeltaMSL);
dtide=min(Trange,dtide);



%Vegetation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if VEGETATION==1;
Zlev=z-msl;
VEG=Zlev>dBlo;
B=4*(Zlev-dBup).*(dBlo-Zlev)/(dBlo-dBup)^2;B(Zlev>dBup)=0;B(Zlev<dBlo)=0;  %B=(dBup-lev)/(dBup-dBlo);B(lev>dBup)=0;B(lev<dBlo)=0;
else;VEG=A*0;B=A*0;end
B=B.*(S==0);%no biomass where there are ponds
VEG=VEG.*(S==0);%no vegeation where there are ponds %%%OCIOO!!! In some cases I did NOT put theVEG=VEG.*(S==0), so that the
%%%ponds count as vegeated area below (i.e., the sediment transport, and %%%maybe the creep). Don't remember why. Maybe for stability? Check lower
%%%tidal range???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Manning%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MANN=0*A+Cb;
%dragadapatationY=0.3;
%MANN(VEG==1)=Cb+min(1,max(0,Zlev(VEG==1)-dBlo)/dragadapatationY)*(Cv-Cb);
MANN(VEG==1)=Cv;
%MANN(VEG==1 & S==0)=Cv;??????????????????????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%River flow
if computeriver==1;
%     Umouth=NaN; %not used, diocane
%     Uo=A*0+1;%first attempt of velocity
%     [UR,URy,URx,q,hriver]=riverFLOWiter(A,Cb,max(0,h),dx,Qmouth,Umouth,Uo);%hriver=max(hriver, max(0,-z));
%     %[UR,URy,URx,q,hriver]=riverFLOWiter(A,Cb,max(0.5,h),dx,Qmouth,Umouth,Uo);%hriver=max(hriver, max(0,-z));
%         if riverwaterlevel==1
%             niter=5-floor(rand(1)+0.5);
%             for i=1:niter
%             hpRIV=(hpRIV+hriver)/2;
%             [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV); 
%             Uo=(Uo+UR)/2;
%             [UR,URy,URx,q,hriver]=riverFLOWiter(A,Cb,h,dx,Qmouth,Umouth,Uo);%hriver=max(hriver, max(0,-z));
%             %[UR,URy,URx,q,hriver]=riverFLOWiter(A,Cb,max(0.5,h),dx,Qmouth,Umouth,Uo);%hriver=max(hriver, max(0,-z));
%             end
%         else;hpRIV=A*0;end
%         hpRIV=hriver;
else;UR=0*A;URy=0*A;URx=0*A;end


%%%%%%%%%%%%%%%%%%%%%%%%%%
if depthlimiterflow_withVEG==1
h(VEG==1 & h<depthlimiterflow_withVEGVALUE)=depthlimiterflow_withVEGVALUE;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%

%CORRECTION RIVER FOR MOMENTUM ADVECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeriver==1 & rivermomemntumcorrection==1
[Vx,Vy]=correctvelocity4MOMENTUMebbflood(N,M,periodic,h,A,MANN,dx,URx,URy);
    %kind of EBB???
    URX=URx;a=find(URX.*Vx>0);URX(a)=URX(a)+Vx(a);
    URY=URy;a=find(URY.*Vy>0);URY(a)=URY(a)+Vy(a);
    URx=URX;URy=URY;
    UR=sqrt(URX.^2+URY.^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if riverwaterlevel==0
h(A==10)=msl-z(A==10);% the river is not influenced by tide or Hsurge (these two!)
end

%Pre-calculate sediment input at river mouth
if computeriver==1;
%[E,Ceq]=totalsedimenterosionSANDsine(h(A==10),0,ss,d50_1,ws1,1,0,1,0,1,URx(A==10),fMFswell,fMFsea,fMFriver,kro,MANN);
[E,Ceq]=totalsedimenterosionSANDsine(h(A==10),hlimC,0,ss,d50_1,ws1,1,URx(A==10),fMFriver,kro,MANN,VEG,U*0);          
[QsmouthSAND]=Ceq.*h(A==10).*URx(A==10);
QsmouthMUD=URx(A==10)*co2mouth*fMFriver.*h(A==10);
QsmouthSAND=sum(QsmouthSAND);
QsmouthMUD=sum(QsmouthMUD);
else;QsmouthSAND=0;QsmouthMUD=0;end



%Tide&Surge flow%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computetide==1;
%dtide(S==1)=min(dtide(S==1),(Trange/2)*0.2);%to limit the effect of ponds on the tidal prism. Speth 2019    %%%%dtide(S==1)=0;%tolgi tutta l aqaua del tidal prism dalle pozze%to limit the effect of ponds on the tidal prism. Speth 2019


if calculateponddynamics==1;
%DIF is the impounded water depth, need to subtract to the imput discharge from the ponded area (THE WATER REMAINS THERE!!!)
pondleakage=0.2;
DIF=max(DIF,0); %you cannot impound a negative water depth!!! This happens because of the trick used to swap the cell during the pond floodin
%dtideI=Trange/2-(z-msl)-DIF;%THSI IS JUST TO PLOT
dtide(S==1)=max(0, Trange/2-(z(S==1)-msl)-max(0,DIF(S==1)-pondleakage));
%reduce the hydperperio din the connecte dpond, becuase they do not exahcneg water as much as their wwater depth. Some of that depth is as if it was made of soil. only the top layer count as moving water!
[~,~,fTidePOND]=getwaterdepth(Trange,Hsurge,msl,z+DIF,kro,hpRIV,tempdeltaMSL);fTide(S==1)=fTidePOND(S==1);
%fTide(S==1)=0.01;%
end

%CALCULATE TIDAL PRSIM
DHeff=dtide+alpha_surge*min(Hsurge,dsurge);%the water level excusrion for the total water prismmin

%multiplier to make fake tidal prism...
%DHeff=min(3*Trange,dtide)+alpha_surge*min(Hsurge,dsurge);%the water level excusrion for the total water prismmin


[U,Uy,Ux]=tidalFLOW(A,MANN,h,ho,dHW,dx,DHeff,Ttide,periodic,kro);
U(A==10)=0;Ux(A==10)=0;Uy(A==10)=0;
else;U=0*A;Uy=0*A;Ux=0*A;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %CORRECTION FOR CURVATUE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if curvaturecorrection==1;
MASKCURVATURE=(VEG==0);

%normal velocity vector
Unx=-Uy;Uny=Ux;
Undiffx=-Uy;Undiffy=Ux;

Uo=sqrt(Ux.^2+Uy.^2);
Uxo=Ux./Uo;Uyo=Uy./Uo;
Uxo(MASKCURVATURE==0)=NaN;Uyo(MASKCURVATURE==0)=NaN;

cur=curlNAN(Uyo,Uxo)/dx;%.*U.*h;%.*U;%.*h.*U;
cur=cur.*U.*h;%.*U;%.*h.*U;
cur=sqrt(abs(cur)).*sign(cur);
     cur=diffusecurvature(isfinite(cur),cur,a_diffusecur,dx,Unx,Uny);
     %cur=cur.*U;%
     PLT.cur=cur;
     cur(isnan(cur))=0;
     %cur=sign(cur).*sqrt(abs(cur))*0.1;

Uno=sqrt(Unx.^2+Uny.^2);Unx=Unx./Uno.*cur;Uny=Uny./Uno.*cur;
   % Unx=diffusecurvature(isfinite(cur),Unx,a_diffusecur,dx,Unx,Uny);
   % Uny=diffusecurvature(isfinite(cur),Uny,a_diffusecur,dx,Unx,Uny);
   % cur(isnan(cur))=0;
%     figure;imagesc(Unx'.*(S'==0));caxis([-0.001 0.001]);
%     figure;imagesc(Uny'.*(S'==0));caxis([-0.001 0.001]);
%     pause
Unx(isnan(Unx))=0;Uny(isnan(Uny))=0;


% UCx=modifyflowcurvatureBASIC(MASKCURVATURE==1,Ux.*MASKCURVATURE,advectflow,dx,Unx.*U/0.2,Uny.*U/0.2,MASKCURVATURE);
% UCy=modifyflowcurvatureBASIC(MASKCURVATURE==1,Uy.*MASKCURVATURE,advectflow,dx,Unx.*U/0.2,Uny.*U/0.2,MASKCURVATURE); 
% UC=sqrt(UCx.^2+UCy.^2);


%UC=modifyflowcurvatureBASIC(MASKCURVATURE==1,U.*MASKCURVATURE,advectflow,dx,Unx.*U/0.2,Uny.*U/0.2,MASKCURVATURE); 
%UC=modifyflowcurvatureBASIC(MASKCURVATURE==1,U.*MASKCURVATURE,advectflow,dx,Unx.*abs(Uy)/0.2,Uny.*abs(Ux)/0.2,MASKCURVATURE); 
%UC=modifyflowcurvatureBASIC(MASKCURVATURE==1,U.*MASKCURVATURE,advectflow,dx,Unx.*U/0.2,Uny.*U/0.2,MASKCURVATURE); 
%factorU=modifyflowcurvatureBASICdiffuse(MASKCURVATURE==1,double(MASKCURVATURE),advectflow,dx,Unx,Uny,MASKCURVATURE,Undiffx,Undiffy,a_diffuseadvection); 

factorU=modifyflowcurvatureBASICdiffuse(MASKCURVATURE==1,double(MASKCURVATURE),advectflow,dx,Unx,Uny,MASKCURVATURE,Undiffx,Undiffy,a_diffuseadvection); 
%factorU=modifyflowcurvatureBASICdiffuse(MASKCURVATURE==1,double(MASKCURVATURE),advectflow,dx,Unx.*U,Uny.*U,MASKCURVATURE,Undiffx,Undiffy,a_diffuseadvection); 
%factorU=modifyflowcurvatureBASICdiffuse(MASKCURVATURE==1,double(MASKCURVATURE),advectflow,dx,Unx,Uny,MASKCURVATURE,Undiffx,Undiffy,a_diffuseadvection); 
%factorU=modifyflowcurvatureBASICdiffuse(MASKCURVATURE==1,double(MASKCURVATURE),advectflow,dx,Unx./max(1,h),Uny./max(1,h),MASKCURVATURE,Undiffx,Undiffy,a_diffuseadvection); 
%factorU=modifyflowcurvatureBASICdiffuse(MASKCURVATURE==1,double(MASKCURVATURE),advectflow,dx,Unx./max(1,h),Uny./max(1,h),MASKCURVATURE,Undiffx,Undiffy,a_diffuseadvection); 

%factorU=min(1,factorU);

fS=(factorU-1).*(MASKCURVATURE==1);

%pause
%UC

%facWS= (1-factorU).*(MASKCURVATURE==1);

% Upazzo=diffuseDeltaU(MASKCURVATURE==1,U,-0.001,dx,Undiffx,Undiffy);
% PLT.pazzo=Upazzo;


%flow into the bank- just storing it for now
%deltaUC=max(0,(UC).*(MASKCURVATURE==0));
%deltaUC=max(0,(UC-U).*(MASKCURVATURE==1));

Ubase=U;
PLT.Ubase=Ubase;
% %Actual flow modification
% Ux(MASKCURVATURE==1)=UCx(MASKCURVATURE==1);
% Uy(MASKCURVATURE==1)=UCy(MASKCURVATURE==1);
% U=sqrt(Ux.^2+Uy.^2);

%Actual flow modification
%U(MASKCURVATURE==1)=U(MASKCURVATURE==1).*factorU(MASKCURVATURE==1);
%fS=fS./max(1,h);
%U(MASKCURVATURE==1)=max(0,U(MASKCURVATURE==1)+fS(MASKCURVATURE==1));

%PLT.deltaUC=deltaUC;

%factorU=-factorU+2;
factorU(MASKCURVATURE==1)=min(1,factorU(MASKCURVATURE==1).^0.5);
factorU(MASKCURVATURE==0)=1;
[U,Uy,Ux]=tidalFLOW(A,MANN./factorU,h,ho,dHW,dx,DHeff,Ttide,periodic,kro);

PLT.fS=fS;
PLT.factorU=factorU;
%STOCASTIC BANK EROSION
if flowbankerosion==1
MASK=MASKCURVATURE;
%[deltaY1,deltaY2,deltaY3,Pedge,Y2OX,EdgeERY1,EdgeERY2,EdgeERY3]=EdgeerosionCURVEDFLOW(deltaUC*dx,z,a_bankerosion,999,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
[deltaY1,deltaY2,deltaY3,PedgeBANK,Y2OX,EdgeERY1,EdgeERY2,EdgeERY3]=EdgeerosionBANKPUSH(abs(cur),Unx,Uny,z,a_bankerosion,999,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
Y1=Y1-deltaY1;Y2=Y2-deltaY2;Y3=Y3-deltaY3;%erode the mardsh edge fully
EDGESED=diffuseedgesediments((A==1),EdgeERY2,1000*h,dx); %put eroded sediment only in chnannels, not on marsh
Y2=Y2+EDGESED;
end

%PLT.PedgeBANK=PedgeBANK;
%PLT.MASKCURVATURE=MASKCURVATURE;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
MASKflow=double(VEG==0);
MASKflow1=bwmorph(MASKflow,'bridge');
MASKflow2=bwmorph(MASKflow1,'diag');
fU=0.2;
Ufmax=max(U/fU,[U(:,1) U(:,1:end-1)]);
Ufmax=max(Ufmax,[U(:,2:end) U(:,end)]);
Ufmax=max(Ufmax,[U(1,:); U(1:end-1,:)]);
Ufmax=max(Ufmax,[U(2:end,:); U(end,:)]);
U(MASKflow2==1 & MASKflow==0)=fU*Ufmax(MASKflow2==1 & MASKflow==0);
%MANN(VEG==1)=Cb;%!!!!!!!!!!!!!!!!!!!ATTENTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     


%CORRECTION FOR MOMENTUM ADVECTION$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ebbfloodcorrection==1
[Vx,Vy]=correctvelocity4MOMENTUMebbflood(N,M,periodic,h,A,MANN,dx,Ux,Uy);
    %EBB
    VebbX=Ux;a=find(VebbX.*Vx>0);VebbX(a)=VebbX(a)+Vx(a);
    VebbY=Uy;a=find(VebbY.*Vy>0);VebbY(a)=VebbY(a)+Vy(a);
    Vebb=sqrt(VebbX.^2+VebbY.^2);
    %FLOOD
    VfloodX=-Ux;a=find(VfloodX.*Vx>0);VfloodX(a)=VfloodX(a)+Vx(a);
    VfloodY=-Uy;a=find(VfloodY.*Vy>0);VfloodY(a)=VfloodY(a)+Vy(a);
    Vflood=sqrt(VfloodX.^2+VfloodY.^2);
    
    %residual currents
    if residualcurrents==1;
    UresX=-(VebbX+VfloodX);UresY=-(VebbY+VfloodY);
    UresX(A==2 | A==0)=0;UresY(A==2 | A==0)=0;
    UresX(1:2,:)=0;UresX(end-1:end,:)=0;%to conserve mass at the open boundary. Need to get velocity zero next to A==2 (upwind!)
    else;UresX=0;UresY=0;
    end
else;Vebb=U;Vflood=U;UresX=0;UresY=0;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Swell waves%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeSwellwave==1;   
       
    hwave=ho;%just to isolate this water depth and not mess up    
%     if Hsurge>0
%     hwaverunup=(max(1,h)-h);
%     hwave=hwave+hwaverunup;
%     wlo=wlo+hwaverunup;
%     end   
   
        %THIS IS THE OLD SINGLE FREQUENCY AND SINGLE DIRECTION
         % kwave=wavek(1/Tp_swell,hwave);
         % [Hs]=SwellWaves(A,AW,Ho,N,M,hwave,Tp_swell,kwave,dx,periodic,angleSWELL,gridDIR);  
    if multifrequency==1; %Case multiple frequency
        [Tperiodi Ejonswap]=getJONSWAPspectrum(Tp_swell,Ho,[1 1.5 2 2.5]);
        [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,hwave,ho,hwSwell_lim,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap,nrefrac,wavediffraction);
    elseif multifrequency==0 %Case single frequency  
         %%%%%%%%%%Tperiodi=Tp_swell;Ejonswap=1;[Hs,waveANGLE,wavePERIOD,PWswell]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,hwave,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap);
        [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave,deltaPW]=SwellWavesMultiDirection(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,hwave,ho,hwSwell_lim,Tp_swell,dx,periodic,angleSWELL,gridDIR,nrefrac,wavediffraction);
    end
   waveANGLE(isnan(waveANGLE))=0;
   %Hs=0*h;Hs(h>0.5)=Ho;%to reproduce ortiz results
   %figure;
   %subplot(1,2,1);imagesc(PWswell);set(gca,'YDir','normal');colormap('jet');
   %subplot(1,2,2);imagesc(deltaPW);set(gca,'YDir','normal');colormap('jet');caxis([0 20]);
   %pause
          if computesand==1;
           [QsWslope,QsWon]=WaveSedimentTransport(Hs,hwave,kwave,rhos,N,M,wavePERIOD,dx,ss,ws1,hwSwell_lim,fTide);
%            QsWon(hwave<=hwSwelltransport_lim)=0;
%            QsWslope(hwave<=hwSwelltransport_lim)=0;
%            deltaPW(hwave<=hwSwelltransport_lim)=0;   
           QsWon(hwave<=hwSwell_lim)=0;
           QsWslope(hwave<=hwSwell_lim)=0;
           deltaPW(hwave<=hwSwell_lim)=0;
           
           %PWswell(hwave<=hwSwell_lim)=0;
          end
else;QsWslope=A*0;QsWon=A*0;Hs=A*0;wavePERIOD=A*0;Uwave=A*0;waveANGLE=A*0;PWswell=0*A;deltaPW=0*A;end
%if periodic==0
%mA=mean(waveANGLE(:));waveANGLE(waveANGLE*mA<1)=0.001*sign(mA);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%swell wave effect on destorying marshes
B(Hs>0.1)=0;
VEG(Hs>0.1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%SeaWaves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if computeSeaWaves==1
%hwave=ho;   %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
MASK=0*A+1;
MASK(ho<=hwSea_lim | A==0 | VEG==1 | Zlev>=dBlo)=0;
[Uwave_sea,Tp_sea,Hsea,Fetch,kwave,PWsea]=SeaWaves(h,angleWIND,hwSea_lim,Trange,wind,MASK,64,dx); %72,dx
Uwave_sea=Uwave_sea.*(VEG==0 & S==0); Hsea=Hsea.*(VEG==0 & S==0); %vegetation effect and no waves in isolated pond 9because we also redcued ws!!1)%Uwave_sea=Uwave_sea.*(VEG==0); Hsea=Hsea.*(VEG==0); %vegetation effect
    if computesand==1;[QsWslope_sea]=WaveSedimentTransport(Hsea,h,kwave,rhos,N,M,Tp_sea,dx,ss,ws1,hwSea_lim,fTide);QsWslope_sea(Hsea==0)=0;end
else;Uwave_sea=0*A;Tp_sea=0*A;Hsea=0*A;Fetch=0*A;QsWslope_sea=0*A;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%Wave-induced edge erosion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (computeEdgeErosionSea==1 | computeEdgeErosionSwell==1)%%%MASK=0*A+1;MASK(h<hwSea_lim | A==0)=0;
PW=A*0;
    if computeEdgeErosionSea==1
    PW=PW+PWsea*fMFsea.*fTide; %Wave power reduction for hydroperiod
    end
    if computeEdgeErosionSwell==1;
    PW=PW+PWswell.*fMFswell.*fTide;
    end       
[deltaY1,deltaY2,deltaY3,Pedge,Y2OX,EdgeERY1,EdgeERY2,EdgeERY3]=EdgeerosionSTRAT_3sedimentsXXX(PW,z,aw,maxedgeheight,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
Y1=Y1-deltaY1;Y2=Y2-deltaY2;Y3=Y3-deltaY3;%erode the mardsh edge fully

%Redistribute the eroded sediment 
%EDGESED=diffuseedgesediments((A==1),EdgeERY2,1*h,dx); %Original ADR 2020 article 
%EDGESED=diffuseedgesediments((A==1),EdgeERY2,1*ho,dx); 
EDGESED=diffuseedgesediments((A==1),EdgeERY2,1*ho,dx);% SIII!
Y2=Y2+EDGESED;
else;Pedge=A*0;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%MORPHODYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%CURRENT-DRIVEN TRANSPORT (Tide and River)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (computetide==1 | computeriver==1)
    
    U(A==10)=0;Uwave(A==10)=0;Uwave_sea(A==10)=0; %in the river mouth only resuspension from river flow
    
    %(1)Total sediment resupension SAND
    if computesand==1;
            if ebbfloodcorrection==1
                 [E1E]=totalsedimenterosionSANDsine(h,hlimC,Vebb,ss,d50_1,ws1,fTide,UR,fMFriver,kro,MANN,VEG,Uwave,wavePERIOD);
                 [E1F]=totalsedimenterosionSANDsine(h,hlimC,Vflood,ss,d50_1,ws1,fTide,UR,fMFriver,kro,MANN,VEG,Uwave,wavePERIOD);
                 E1=(E1E+E1F)/2;
            else
                 [E1]=totalsedimenterosionSANDsine(h,hlimC,U,ss,d50_1,ws1,fTide,UR,fMFriver,kro,MANN,VEG,Uwave,wavePERIOD);
            end            
            E1(A==2)=0; %needed for b.c.
    end;
    
    %(2)Total sediment resupension MUD
    if computemud==1
                    if computeSwellwave==1;%Swell waves for MUD only!
                    Uwave4mud=Uwave.*(VEG==0); %vegetation effect. Plants put to zero wave erosion
                    else;QsWslope=zeros(N,M);QsWon=zeros(N,M);BRK=zeros(N,M);Uwave4mud=zeros(N,M);Hs=zeros(N,M);Tp_swellMUD=1;HswellMUD=0;end
            if ebbfloodcorrection==1
                [E2E,E2tideE]=totalsedimenterosionMUDsine(Vebb,MANN,VEG,fTide,UR,Uwave_sea,Uwave4mud,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,tcrgradeint,leveltauincrease,taucrVEG,me,h,Zlev,NaN,computeSeaWaves,computeSwellwave,computeriver,limitertauC);
                [E2F,E2tideF]=totalsedimenterosionMUDsine(Vflood,MANN,VEG,fTide,UR,Uwave_sea,Uwave4mud,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,tcrgradeint,leveltauincrease,taucrVEG,me,h,Zlev,NaN,computeSeaWaves,computeSwellwave,computeriver,limitertauC);
                E2=(E2E+E2F)/2;E2tide=(E2tideE+E2tideF)/2;
            else
                [E2,E2tide]=totalsedimenterosionMUDsine(U,MANN,VEG,fTide,UR,Uwave_sea,Uwave4mud,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,tcrgradeint,leveltauincrease,taucrVEG,me,h,Zlev,NaN,computeSeaWaves,computeSwellwave,computeriver,limitertauC);
            end      
            E2(A==2)=0; %needed for b.c.
          
    %(3)Total sediment resupension ORGANIC
    E3=E2;
    end

    
    %Erosion limiters
    if computesand==1
    E1=E1.*fracY1;%Reduced for fraction of sediemnt
    Elimit=max(0,Y1*conces)/dt*rbulk1;a=find(E1>Elimit);E1(a)=Elimit(a);%this is the limit to avoid E to scour more than Y1 or Y2. 
    end
    if computemud==1;    
       % conces=1;%how much to extra erode, a parameter
    E2=E2.*fracY2;    
    Elimit=max(0,Y2*conces)/dt*rbulk2;a=find(E2>Elimit);E2(a)=Elimit(a);
    E3=E3.*fracY3;   
    Elimit=max(0,Y3*conces)/dt*rbulk2;a=find(E3>Elimit);E3(a)=Elimit(a);
    end
         
    
    %Advection-Diffusion Sediment transport
    if computesand==1;
        WS=A*0+ws1;
         [EmD1,SSM,FLX1]=sedtran(0,h,A,SPCLcell,0,DiffSsand,h,ho,E1,WS,dx,dt,rbulk1,co1,Ux,Uy,FLX1,fTide,Ttide,URx,URy,UresX,UresY,periodic,computeriver,computetide,residualcurrents,kro,0,0,QsmouthSAND);  
     else;EmD1=0*A;SSM=0*A;end
    SSC1=SSM./h;
    if computemud==1;
        WS=A*0+ws2;
        WS(VEG==1)=wsB;           
        WS(S==1)=ws2;%SHOULD NOT BE NECEEARY BECUASE VEG alreeady set equal to zero where S=1 (see above).  ->Do not add the vegetation settling velocity in the ponds! %WS(S==1)=0.000000000001;%E2(S==1)=0;
        %WS=WS+facWS*5/1000;
        [EmD2,SSM,FLX2]=sedtran(1,h,A,SPCLcell,DoMUD,DiffSmud,h,ho,E2,WS,dx,dt,rbulk2,co2,Ux,Uy,FLX2,fTide,Ttide,URx,URy,UresX,UresY,periodic,computeriver,computetide,residualcurrents,kro,1,co2mouth*fMFriver,QsmouthMUD);   
    else;EmD2=0*A;SSM=0*A;end
    SSC2=SSM./h;%./fTide;  %devi metter i fTide per farti plottare la b.c quando il fondo e' sopra il MLW (ftide<1)
    if VEGstratigraphy==1;
        WS=A*0+ws2;WS(VEG==1)=wsB;
        [EmD3,SSM,FLX3]=sedtran(1,h,A,SPCLcell,DoMUD,DiffSmud,h,ho,E3,WS,dx,dt,rbulk2,co3,Ux,Uy,FLX3,fTide,Ttide,URx,URy,periodic,computeriver,computetide,kro,1,0,0);
    else;EmD3=0*A;SSM=0*A;end
    SSC3=SSM./h;%./fTide;  %devi metter il fTide per farti plottare la b.c quando il fondo e' sopra il MLW (ftide<1)
    
else;EmD1=0*A;Qs1=0*A;tideE1=A*0;E1=A*0;EmD2=0*A;SSC=0*A;SSC1=A*0;SSC2=A*0;SSC3=A*0;end 
%end of compute tide and river%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Bed evolution erosion/depositon from tidal and river transport
z=imposeboundarydepth(A,z,optionBC,NaN);
if computesand==1;    Y1=Y1-dt*EmD1;end
if computemud==1;     Y2=Y2-dt*EmD2;end
if VEGstratigraphy==1;Y3=Y3-dt*EmD3;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%Current drive bank erosion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%z=zs-(Y1+Y2+Y3);
% if compute_currentbankerosion==1
% MASK=0*A+1;MASK((z+msl)<Trange/2 | A==0)=0;
% PW=20*U;%100
% [deltaY1,deltaY2,deltaY3,Pedge,Y2OX]=Bankerosion(PW,z,aw,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
% Y1=Y1-deltaY1;Y2=Y2-deltaY2;Y3=Y3-deltaY3;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Organic accretion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=B.*(S==0);
if AccreteOrganic==1
if VEGstratigraphy==1;
Y3=Y3+B*Korg*dt; %accrete the organic
else; %put it with mud or sand
    if VEGonsand==0
    Y2=Y2+B*Korg*dt; % putorganic on mud!!!
    else
    Y1=Y1+B*Korg*dt; % putorganic on sand!!!
    end
end
KBTOT=KBTOT+sum(B(A==1))*Korg*dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %Water table
% zbed=z+msl;
%   facdrain=0.05;
%   xdrain1=-0.1+bwdist(zbed>0)*facdrain;
%   xdrain2=-0.1+0.2+bwdist(zbed>-0.2)*facdrain;
%   xdrain3=-0.1+0.4+bwdist(zbed>-0.4)*facdrain;
%   xdrain4=-0.1+0.6+bwdist(zbed>-0.6)*facdrain;
%   xdrain5=-0.1+0.8+bwdist(zbed>-0.8)*facdrain;
%   xdrain6=-0.1+1+bwdist(zbed>-1)*facdrain;
%   xdrain7=-0.1+1.2+bwdist(zbed>-1.2)*facdrain;
%   xdrain=min(min(min(min(min(min(min(-zbed,xdrain1),xdrain2),xdrain3),xdrain4),xdrain5),xdrain6),xdrain7);
%   hdrain=max(0,min(1,-zbed-xdrain));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% %Organic carbon oxydation
% Koxidation=0.01/365;
% KBTOT=KBTOT-sum(Y3(A==1 & Y3>0).*hdrain(A==1 & Y3>0)*Koxidation*dt);
% Y3(A==1 & Y3>0)=Y3(A==1 & Y3>0)-Y3(A==1 & Y3>0).*hdrain(A==1 & Y3>0)*Koxidation*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%Compaction ->simulated as lowering of the basal level
%zPOT=2;  This is the potential for compaction (in meters). If greater than
%zero can compact more by lowering the water table. If 0 then it has
%compacted 100%
%Kcompaction=0.01/365;
%zb(A==1)=zb(A==1)+hdrain(A==1)*Kcompaction*dt;
%zb(A==1)=zb(A==1)+hdrain(A==1).*zPOT(A==1)*Kcompaction*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% %%%%%%%%%%%%%%%%%%%%%%Wave-induced edge erosion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (computeEdgeErosionSea==1 | computeEdgeErosionSwell==1)
% %MASK=0*A+1;MASK(h<hwSea_lim | A==0)=0;
% PW=A*0;
%     if computeEdgeErosionSea==1
%     PW=PW+PWsea*fMFsea.*fTide; %Wave power reduction for hydroperiod
%     end
%     if computeEdgeErosionSwell==1;
%     PW=PW+PWswell.*fMFswell.*fTide;
%     end       
% [deltaY1,deltaY2,deltaY3,Pedge,Y2OX]=EdgeerosionSTRAT_3sediments(PW,z,aw,maxedgeheight,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
% %[deltaY1,deltaY2,deltaY3,Pedge,Y2OX,EdgeERY1,EdgeERY2,EdgeERY3]=EdgeerosionSTRAT_3sedimentsXXX(PW,z,aw,maxedgeheight,fox,dt,dx,MASK,A,fracY1,fracY2,fracY3,Y2OX);
% Y1=Y1-deltaY1;Y2=Y2-deltaY2;Y3=Y3-deltaY3;
% else;Pedge=A*0;end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%BED EVOLUTION VERTICAL FLUXES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update the  bed using: 1)Current transport 2)Edge Erosion 3)Organic growth
%this will allow to next compute the evolution by divergence!!!
z=zs+(Y1+Y2+Y3);
znew=z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%BED EVOLUTION DIVERGENCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EVOLUTION OF Y1
if computesand==1;
Yreduction=(1-min(1,exp(-10*(Y1)))).*fracY1;  %ONLY FOR THE WAVE TRANSPORT, DIOCANE
Qs1=E1./(ws1*3600*24).*(U+UR*0.3).*max(hlimC,h); %kg/m/s
%wlo is the water level in which the dpeht h can actually go to zer (h is
%the water level in which the depth is at minimum kro, do avoid explidign
%with flow and disperive transport. wlo is used in wave transpor
%wlo=wlo+hwaverunup;%add extra water level due to wave run up. Already
%added where you calculate the swell waves around line 189
[znew,FQsW_L,FQsW_R,longshore,angleshore]=bedevolutionDIVlongshore(fMFswell*deltaPW,U,fTide,A,AW,z,Zlev,wlo,ho,Yreduction,N,M,dt,dx,Trange,Qs1,rbulk1,alphaSAND,hwSwelltransport_lim,PWfactorlongshore,computeSwellwave,fMFswell*QsWslope+fMFsea*QsWslope_sea*downslopeSANDseawaves,fMFswell*QsWon,angleSWELL,waveANGLE,Active,periodic,optionBC,gridDIR,FQsW_L,FQsW_R);
deltaY1=z-znew;Y1=Y1-deltaY1;
else;longshore=0;angleshore=0;end

%EVOLUTION OF Y2 
if computemud==1;
z=znew;  %NEED TO RE-UDPATED FROM THE Y1
Yreduction=(1-min(1,exp(-10*(Y2)))).*fracY2;%fracY2; %this is just for the marsh stratigraphy%FORSE PRIMA ERA PER UNA NO PER 10    %Yreduction=fracY2;  %this is for most of the simulations big basin % Yreduction(Yreduction<0.5)=0;
%Qs2=E2./(ws2*3600*24).*(U+UR*0.3).*max(hlimC,h); %kg/m/s

%Qs2=E2tide./(ws2*3600*24).*(U+UR*0.3).*max(hlimC,h); %kg/m/s %[hcorrected]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV);  %VEG=(z-msl)>dBlo;
Qs2=E2tide./(ws2*3600*24).*(U+UR*0.3).*max(0,h); %kg/m/s %[hcorrected]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV);  %VEG=(z-msl)>dBlo;

%%%ALOW TO CREEP IN PONDS.!!!
%VEG=VEG.*(S==0);%Qs2=Qs2.*(S==0); 
%znew=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,dx,dt,VEG,A*0,Qs2,rbulk1,alphaMUD);  %MUD CREEP  MARSH
znew=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,crbank,dx,dt,VEG,S,Qs2,rbulk2,alphaMUD,facQsbank,Ubase);%,deltaUC,a_bankcreep);  %MUD CREEP  MARSH
deltaY2=z-znew;
deltaY2(A==2)=0;  %DO NOT UPDATE THE BOUNDARY
Y2=Y2-deltaY2;
end
 
%EVOLUTION OF Y3 
if computemud==1;
    if VEGstratigraphy==1;
    z=znew;  %NEED TO RE-UDPATED FROM THE Y1
    Yreduction=fracY3;%fracY3; (1-min(1,exp(-1*(Y3)))).*%Yreduction(Yreduction<0.1)=0;%to reduce some starneg probelms with creep of two sediment swiht stratigraphy. negative creep!
    znew=bedcreepponds(z,A,Active,Yreduction,crMUD,crMARSH,crbank,dx,dt,VEG,S);  %MUD CREEP
    deltaY3=z-znew;Y3=Y3-deltaY3;
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if evolvestratigraphy==1;
[Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr]=stratigraphy2D_3sediments(A,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr,nlyr,dlyr,tlyrU,tlyrD);
end

%%%%IMPOSE "NEUMAN" boundary condition for morphodynamics%%%%%%%%%%%%%%%%%%
%Traslate the first boundary cell bed elevetion
if imposeseaboundarydepthmorphoALL==1;
[plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=seaboundaryNeumanbedelevationALLBOUNDARY(A,plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3);
else
    %remeber that zb has the oppisite sign!!!
    if computesand==1
        Y1(A==2)=Y1(A==2)+RSLR*dt;%%%ADDED SEPT 2019
    else
        Y2(A==2)=Y2(A==2)+RSLR*dt;%%%ADDED SEPT 2019
    end
end


%%%%%%%%%%%%%%%%%END OF MORPHODYNAMICS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % %%%%%%%%%%%%%%%%
% znew=z;
% znew(end-50:end,:)=-5+msl;
% deltaY2=z-znew;
% Y2=Y2-deltaY2;
% %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% NUMERICAL
%%%%%%%%%% CHECKS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxdeltaz=prctile(abs(znew(:)-zoriginal(:)),99.9);
%ORGINAL VERSION muchup=max(0,max(0,znew-zoriginal)-dHW);
MMM=(znew-zoriginal).*((znew)>(msl+Trange/2+Hsurge));
muchup=max(0,max(0,znew-zoriginal)).*((znew)>(msl+Trange/2+Hsurge));%modieif on Oct 2019
maxup=max(muchup(:));%maxup=prctile(muchup(:),99.9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%OUTPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IO.Y1=Y1;IO.Y2=Y2;IO.Y3=Y3;
% IO.zb=zb;
% IO.flyr1=flyr1;IO.flyr2=flyr2;IO.flyr3=flyr3;
% IO.flyrb1=flyrb1;IO.flyrb2=flyrb2;IO.flyrb3=flyrb3;
% IO.plyr=plyr;IO.Yb=Yb;IO.msl=msl;
% IO.Active=Active;
% fIO.FLX1=FLX1;fIO.FLX2=FLX2;fIO.FLX3=FLX3;
% fIO.KBTOT=KBTOT;fIO.Y2OX=Y2OX;
% fIO.pondloss=pondloss;
% fIO.FQsW_R=FQsW_R;
% fIO.FQsW_L=FQsW_L;

names = fieldnames(IO);
for i=1:length(names);eval(['IO.' names{i} '=' names{i} ';' ]);end

names = fieldnames(fIO);
for i=1:length(names);eval(['fIO.' names{i} '=' names{i} ';' ]);end

PLT.PW=PW;
PLT.angleshore=angleshore;
PLT.wl=wl;
PLT.wlo=wlo;
%PLT.hwaverunup=hwaverunup;
if calculateponddynamics==1
PLT.DIF=DIF;
PLT.dtide=dtide;
%PLT.dtideI=dtideI;
end
%PLT.Fetch=Fetch;
%PLOT outputs
PLT.U=U;
PLT.Vebb=Vebb;
PLT.Vflood=Vflood;
PLT.UR=UR;
%if computetide==1 | computeriver==1;
PLT.SSC1=SSC1;
PLT.SSC2=SSC2;
PLT.SSC3=SSC3;
%end
PLT.Tp_sea=Tp_sea;
PLT.Uwave_sea=Uwave_sea;
PLT.Hsea=Hsea;
PLT.Fetch=Fetch;
PLT.EmD2=EmD2;
PLT.Pedge=Pedge;
%PLT.SALTC=SALTC;
%PLT.QsWon=QsWon;
PLT.URx=URx;
PLT.URy=URy;
PLT.Ux=Ux;
PLT.Uy=Uy;
    PLT.Hs=Hs;
if computeSwellwave==1;
    PLT.waveANGLE=waveANGLE;
    PLT.wavePERIOD=wavePERIOD;
    %PLT.kwave=kwave;
    PLT.Uwave=Uwave;
end
if ebbfloodcorrection==1
%PLT.Vx=Vx;PLT.Vy=Vy;
    PLT.Vebb=Vebb;PLT.Vflood=Vflood;
    PLT.VebbX=VebbX;PLT.VfloodX=VfloodX;
    PLT.VebbY=VebbY;PLT.VfloodY=VfloodY;
end
PLT.VEG=VEG;
%PLT.MARSH=MARSH;
%PLT.Cmin=Cmin;
%PLT.Cmax=Cmax;
PLT.deltaPW=deltaPW;
%PLT.HCURV=HCURV;
%PLT.CD=CD;
PLT.wsB=wsB;
%PLT.hriver=hriver;
PLT.h=h;
PLT.AC=AC;
PLT.S=S;
PLT.B=B;
PLT.fTide=fTide;
%PLT.fTidePOND=fTidePOND;
%[angleSWELL ]
PLT.EmD1=EmD1;
%PLT.hdrain=hdrain;
%PLT.qswell=qswell;
%PLT.E1=E1;
PLT.hpRIV=hpRIV;
PLT.ho=ho;
PLT.longshore=longshore;

PLT.DHeff=DHeff;
%%%OLD STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%SALT transport
% [SALT]=salttran(h,A,SPCLcell,1,DiffS,h,ho,E2*0,ws2,dx,dt,rbulk2,30,Ux,Uy,FLX2,fTide,Ttide,fMFriver*URx,fMFriver*URy,periodic,computeriver,computetide,kro,1,0,QsmouthMUD);
% SALTC=SALT./h./fTide;
%SALTC=0;


% if imposeseaboundarydepthmorphoEAST==1;%east
% [plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=seaboundaryNeumanbedelevationEAST(plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3);
% end
% if imposeseaboundarydepthmorphoNORTH==1;%north
% [plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=seaboundaryNeumanbedelevation(plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3);
% end
% sumY1=sumSedcolum(Yb,flyrb1,flyr1,dlyr,Y1);sumY1=sum(sumY1(A==1))









% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% z=zs-(Y1+Y2+Y3);%update the bed elevation
% [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(Trange,Hsurge,msl,z,kro,hpRIV);
% 
% % %Swell waves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if computeSwellwave==1;   
%     hwave=ho;%just to isolate this water depth and not mess up
%     
%         %THIS IS THE OLD SINGLE FREQUENCY AND SINGLE DIRECTION
%          % kwave=wavek(1/Tp_swell,hwave);
%          % [Hs]=SwellWaves(A,AW,Ho,N,M,hwave,Tp_swell,kwave,dx,periodic,angleSWELL,gridDIR);   
%     if multifrequency==1; %Case multiple frequency
%         [Tperiodi Ejonswap]=getJONSWAPspectrum(Tp_swell,Ho,[1 1.5 2 2.5]);
%         [Hs,waveANGLE,wavePERIOD,PWswell,kwave]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,hwave,hwSwell_lim,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap,nrefrac,wavediffraction);
%     elseif multifrequency==0 %Case single frequency  
%          %%%%%%%%%%Tperiodi=Tp_swell;Ejonswap=1;[Hs,waveANGLE,wavePERIOD,PWswell]=SwellWavesMultiDirectionMultiFrequency(A,AW,Ho,N,M,hwave,Tperiodi,dx,periodic,angleSWELL,gridDIR,Ejonswap);
%         [Hs,waveANGLE,wavePERIOD,PWswell,kwave,Uwave]=SwellWavesMultiDirection(A,AW,Ho,N,M,Cbr,Cbed,wavefrictionCollins,hwave,hwSwell_lim,Tp_swell,dx,periodic,angleSWELL,gridDIR,nrefrac,wavediffraction);
%     end
%    waveANGLE(isnan(waveANGLE))=0;
%    %Hs=0*h;Hs(h>0.5)=Ho;%to reproduce ortiz results
%    
%           if computesand==1;
%            [QsWslope,QsWon]=WaveSedimentTransport(Hs,hwave,kwave,rhos,N,M,wavePERIOD,dx,ss,ws1,hwSwell_lim,fTide);
%            QsWon(hwave<=hwSwell_lim)=0;QsWslope(hwave<=hwSwell_lim)=0;
%           end
% else;QsWslope=A*0;QsWon=A*0;Hs=A*0;wavePERIOD=A*0;Uwave=A*0;waveANGLE=A*0;end
% 
% %SeaWaves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if computeSeaWaves==1
% %hwave=ho;   %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% MASK=0*A+1;
% MASK(ho<=hwSea_lim | A==0 | VEG==1 | lev>=0)=0;
% [Uwave_sea,Tp_sea,Hsea,Fetch,kwave,PWsea]=SeaWaves(h,angleWIND,hwSea_lim,Trange,wind,MASK,64,dx);
% Uwave_sea=Uwave_sea.*(VEG==0); Hsea=Hsea.*(VEG==0); %vegetation effect
% [QsWslope_sea]=WaveSedimentTransport(Hsea,h,kwave,rhos,N,M,Tp_sea,dx,ss,ws1,hwSea_lim,fTide);
% QsWslope_sea(Hsea==0)=0;
% else;Uwave_sea=0*A;Tp_sea=0*A;Hsea=0*A;Fetch=0*A;QsWslope_sea=0*A;end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



