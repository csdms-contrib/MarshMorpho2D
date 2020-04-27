function [E,Etide]=totalsedimenterosionMUDsine(U,MANN,VEG,fTide,UR,Uwave_sea,Uwave,Tp_sea,Tp_swell,fMFswell,fMFsea,fMFriver,taucr,taucrVEG,me,h,lev,TrangeVEG,computeSeaWaves,computeSwellwave,computeRiver);
fUpeak=pi/2;

taucro=U*0+taucr;
taucro(VEG==1)=taucrVEG;


%%%%%%%%%%%%%%%%%%%%%tidal current erosion
ncyc=10;
E=0;
for i=0:ncyc
Utide=U*fUpeak*sin(i/ncyc*pi/2);
tauC=1030*9.81*MANN.^2.*h.^(-1/3).*Utide.^2; 
E=E+1/(ncyc+1)*me.*(sqrt(1+(tauC./taucro).^2)-1);
end
Etide=E;


%Sea waves erosion
if computeSeaWaves==1
Uwave_sea=max(0,Uwave_sea);

ko=0.1/1000*3;
aw=Tp_sea.*Uwave_sea/(2*pi);
fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;
tauWsea=0.5*1030*fw.*Uwave_sea.^2;
E=E+me.*(sqrt(1+(tauWsea./taucro).^2)-1)*fMFsea.*fTide;  
end

%River erosion
if computeRiver==1
%tauCRiver=1030*0.04^2*9.81*h.^(-1/3).*UR.^2;
tauC=1030*9.81*MANN.^2.*h.^(-1/3).*UR.^2; 
E=E+me.*(sqrt(1+(tauC./taucro).^2)-1)*fMFriver.*fTide;  
end


