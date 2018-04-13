function [E]=totalsedimenterosionMUDsine(U,MANN,VEG,fTide,taucr,tcrgradeint,taucrVEG,me,h,lev,TrangeVEG);
fUpeak=pi/2;

taucro=U*0+taucr;
taucro(VEG==1)=taucrVEG;

%increase tcr with depth (asusming an existing vertical distribution. 
%USE with cauton, only ok for simulation of small marsh domain
xi=-lev-TrangeVEG/2;xi(xi<0)=0;
taucro=taucro+xi*tcrgradeint;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%tidal current erosion
ncyc=10;
E=0;
for i=0:ncyc
Utide=U*fUpeak*sin(i/ncyc*pi/2);
tauC=1030*9.81*MANN.^2.*h.^(-1/3).*Utide.^2; 
E=E+1/(ncyc+1)*me.*(sqrt(1+(tauC./taucro).^2)-1);
end

