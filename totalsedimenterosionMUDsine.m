function [E,Etide]=totalsedimenterosionMUDsine(U,MANN,VEG,fTide,UR,Uwave_sea,Uwave,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,tcrgradeint,leveltauincrease,taucrVEG,me,h,lev,TrangeVEG,computeSeaWaves,computeSwellwave,computeRiver,limitertauC);
fUpeak=pi/2;

taucro=U*0+taucr;
taucro(VEG==1)=taucrVEG;


%increase tcr with depth (asusming an existing vertical distribution. 
%USE with cauton, only ok for simulation of small marsh domain
xi=-lev-leveltauincrease;xi(xi<0)=0;%
taucro=taucro+xi*tcrgradeint;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%figure;imagesc(taucro);pause
%%%%%%%%%%%%%%%%%%%%%tidal current erosion
ncyc=10;
E=0;
for i=0:ncyc
Utide=U*fUpeak*sin(i/ncyc*pi/2);
tauC=1030*9.81*MANN.^2.*h.^(-1/3).*Utide.^2; 
tauC(h<0.5)=min(limitertauC,tauC(h<0.5));%limiter for tay, since the smallest water depth kro = 5cm
E=E+1/(ncyc+1)*me.*(sqrt(1+(tauC./taucro).^2)-1);%E=E+1/(ncyc+1)*me.*taucr./taucro.*(sqrt(1+(tauC./taucr).^2)-1);
end
Etide=E;

%%Swell erosion
if computeSwellwave==1
%ko=1/1000;%0.1/1000;%0.1*3
%aw=Tp_swell.*Uwave/(2*pi);
%fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;
fw=0.015;
tauWswell=0.5*1030*fw.*Uwave.^2;
%E=E+me*max(0,tauWswell-taucr)./taucro*fMFswell;
E=E+me.*(sqrt(1+(tauWswell./taucro).^2)-1)*fMFswell.*fTide; %OCIO AI TOLTO QUIESTI fTide IN APRIL 18th 2018
end

%Sea waves erosion
if computeSeaWaves==1
Uwave_sea=max(0,Uwave_sea);

%ko=0.5/1000;%%%
ko=1/1000;
aw=Tp_sea.*Uwave_sea/(2*pi);
fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;
%fw=0.015;
tauWsea=0.5*1030*fw.*Uwave_sea.^2;
%E=E+me*max(0,tauWsea-taucr)./taucro*fMFsea;
E=E+me.*(sqrt(1+(tauWsea./taucro).^2)-1)*fMFsea.*fTide;  %OCIO AI TOLTO QUIESTI IN APRIL 18th 2018
end

%River erosion
if computeRiver==1
%tauCRiver=1030*0.04^2*9.81*h.^(-1/3).*UR.^2;
tauC=1030*9.81*MANN.^2.*h.^(-1/3).*UR.^2; 
E=E+me.*(sqrt(1+(tauC./taucro).^2)-1)*fMFriver.*fTide;  
end


% tauCRiver=1030*0.04^2*9.81*h.^(-1/3).*UR.^2;
 

%E=E+1/(ncyc+1)*me./taucro.*(sqrt(taucr.^2+tauC.^2)-taucr);
%E=E+1/(ncyc+1)*me*max(0,tauC-taucr)./taucro;

%E=me./taucro.*(sqrt(taucr.^2+tauC.^2)-taucr)*fMFtide;


%NONCONTA!!!
%E=me*max(0,tauC-taucr)./taucro*fMFtide;% +me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea +me*max(0,tauCRiver-taucr)./taucro*fMFriver;
%E=me*(sqrt(1+(tauC/taucr).^2)-1)*fMFtide;



%+me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea +me*max(0,tauCRiver-taucr)./taucro*fMFriver;
%%%
%E=h*NaN;
%Ec=me*max(0,tauC-taucr)./taucro*fMFtide;
%Ew=+me*max(0,tauWswell-taucr)./taucro*fMFswell +me*max(0,tauWsea-taucr)./taucro*fMFsea;

%E=E*0;
%E(h<=0.1)=0;
% Ec(h<=0.1)=0;
% Ew(h<=0.1)=0;

%Ew=Ew*0;
%Ec=Ec*0;