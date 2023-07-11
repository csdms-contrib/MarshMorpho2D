function [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(range,Hsurge,msl,z,kro,hpRIV,tempdeltaMSL);


z=z-tempdeltaMSL;%temporary shift to change MSL at ever tide, trick!!!

DH=range+Hsurge; %total water displacement
dHW=max(0,-z+hpRIV+msl+range/2+Hsurge);%water depth at MHW
dtide=max(0,-z+hpRIV+msl+range/2);%water depth at MHW
dsurge=max(0,-z+hpRIV+msl+Hsurge);%water depth at MHW
fTide=min(1,max(10^-5,dHW/DH));%hydroperiod
%h=0.5*(max(0,dHW)+max(0,dHW-DH)); %for the flow  THE ORGINIAL
h=0.5*(max(0,dHW)+max(0,dHW-range)); %for the flow
ho=h; %the original water depth without the limiter %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!


% relax=10;
% hxxx=h(h<kro);
% h(h<kro)=max(0.05,kro*(1-exp(-hxxx*relax))/(1-exp(-kro*relax)));

h(h<kro)=kro;%reference water depth for flow calculation


% kl1=0.;
% kl2=0.5;
% a=find(h<kl2 & h>kl1);
% h(a)=kl2/2+(h(a)-kl1)/(kl2-kl1)*kl2/2;
% h(h<=kl1)=kl2/2;
% 

%hi is the depth that gives the tidal prism %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
%ho=h;
%ho(-z>(hpRIV+msl+range/2+Hsurge))=0; %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!

wl=z+h;
wlo=z+ho;

% figure;
% imagesc(ho)
% pause











%Original November 2019!!!!!!!!!!!!!
% function [h,ho,fTide,dtide,dsurge,dHW,wl,wlo]=getwaterdepth(range,Hsurge,msl,z,kro,hpRIV);
% DH=range+Hsurge; %total water displacement
% dHW=max(0,-z+hpRIV+msl+range/2+Hsurge);%water depth at MHW
% dtide=max(0,-z+hpRIV+msl+range/2);%water depth at MHW
% dsurge=max(0,-z+hpRIV+msl+Hsurge);%water depth at MHW
% fTide=min(1,max(10^-3,dHW/DH));%hydroperiod
% %h=0.5*(max(0,dHW)+max(0,dHW-DH)); %for the flow  THE ORGINIAL
% h=0.5*(max(0,dHW)+max(0,dHW-range)); %for the flow
% ho=h; %the original water depth without the limiter %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% h(h<kro)=kro;%reference water depth for flow calculation
% 
% % kl1=0.;
% % kl2=0.5;
% % a=find(h<kl2 & h>kl1);
% % h(a)=kl2/2+(h(a)-kl1)/(kl2-kl1)*kl2/2;
% % h(h<=kl1)=kl2/2;
% % 
% 
% %hi is the depth that gives the tidal prism %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% %ho=h;
% %ho(-z>(hpRIV+msl+range/2+Hsurge))=0; %CAMBIATO MAGGIO 14 2018!!!!!!!!!!!!!!!!!
% 
% wl=z+h;
% wlo=z+ho;
% 
% % figure;
% % imagesc(ho)
% % pause