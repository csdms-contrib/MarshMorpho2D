function [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(range,msl,z,kro,hpRIV);
DH=range; %total water displacement
dHW=max(0,-z+hpRIV+msl+range/2);%water depth at MHW
dtide=max(0,-z+hpRIV+msl+range/2);%water depth at MHW
dsurge=max(0,-z+hpRIV+msl);%water depth at MHW
fTide=min(1,max(10^-3,dHW/DH));%hydroperiod
h=0.5*(max(0,dHW)+max(0,dHW-range)); %for the flow
ho=h; %the original water depth without the limiter %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!

relax=10;
hxxx=h(h<kro);
h(h<kro)=max(0.02,kro*(1-exp(-hxxx*relax))/(1-exp(-kro*relax)));

wl=z+h;
wlo=z+ho;




