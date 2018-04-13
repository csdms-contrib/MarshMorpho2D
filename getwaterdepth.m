function [h,ho,fTide,dtide,dHW,wl]=getwaterdepth(range,msl,z,kro);
DH=range; %total water displacement
dHW=max(0,z+msl+range/2);%water depth at MHW
dtide=max(0,z+msl+range/2);%water depth at MHW
fTide=min(1,max(10^-3,dHW/DH));%hydroperiod
h=0.5*(max(0,dHW)+max(0,dHW-DH)); %for the flow
%ho=h; %the original water depth without the limiter
h(h<kro)=kro;%reference water depth for flow calculation
ho=h;

%the depth that gives the tidal prism
ho(-z>(msl+range/2))=0;

wl=h-z;

% figure;
% imagesc(ho)
% pause