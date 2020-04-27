clear; close all;clc


%NOTES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qs is the total lateral transport, E is the vertical erosion flux

%lateral boundary options
% 0 is closed (no flux if nothign specified. no-gradient if AW equal to 1 or -1)
% 1 is periodic

%A
%0: not in the domain
%1: a normal cell
%2: the open sea boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('units','normalized','outerposition',[0 0 1 1]) %set(gca,'Color','k')


rng(3)%to get always the same random numbers. you can pick any number in there

P.g=9.81; %gravity [m/s2]
P.rho=1030; %water density [kg/m3] 
P.rhos=2650; %sediment density (quartz-mica) [kg/m3]
P.ss=(P.rhos-P.rho)/P.rho; %relative density
P.kro=0.1;%[m] minimum water depth for hydrodynamics. NEEDS TO BE SMALLER THAN hwSea_lim
P.DiffSmud=1; %[-]coefficient for tidal dispersion [-]. 1 DO NOT CHANGE
P.DoMUD=1;%base diffusivity of suspedned mud [m2/s]. Process not related to tides (e.g. wind and waves, other ocean circulation)

%Sea level rise
P.RSLR=2.5/1000/365;  %from mm/yr to m/day (the time unit is the day!)

%Tide
P.Ttide=12.5/24; %tidal period [day]
P.Trange=0.7; %Tidal Trange [m]
P.TrangeVEG=P.Trange;%tidal range for vegetation [m]. Generally same of tidal range

%Wind for sea waves
P.wind=7;%reference wind speed [m/s]

%Edge erosion
P.aw=0.3/365; %wave edge erodability m/yr/W/m2
P.maxedgeheight=2;
P.fox=0;%fraction of edge eroded material that is oxidized.

%Wind waves and swell numerics
P.hwSwelltransport_lim=1;
P.hwSea_lim=0.2;%0.5; %limiter water deth of sea waves %THIS IS USED TO FIND THE "EDGE"%NEEDS TO BE LARGER THAN KO!!!!

%SSC at the sea boundary
P.co2=60/1000;%40/1000; %Sea boundary SSC for mud [g/l]

%Manning coeffinent unvegeated (same for sand and mud)
P.Cb=0.02;

%Mud
P.d50_2=0.02/1000/1000;%mud grain size [m]
P.ws2=0.2/1000;%0.2/1000;% m/s
P.por2=0.7;P.rbulk2=P.rhos*(1-P.por2);%P.por2=0.7

%Mud parameters
P.me=0.1*10^-4*24*3600;  %per day!!!
P.taucr=0.2;
P.crMARSH=0.1/365;%creep coefficient vegetated
P.crMUD=3.65/365;%creep coeffcinet
P.alphaMUD=0.25; %coefficient for bedload downslope of mud. added April 2019. Similar to P.alphaSAND

%Vegetation parameters
P.dBlo=0;%-0.2;%-0.2;%0;%-(0.237*P.TrangeVEG-0.092);
P.dBup=P.TrangeVEG/2;%-0.2;
P.Cv=0.1;%Manning for vegetated ara
P.wsB=1/1000;%Mud Settling velocity for vegetated ara
P.taucrVEG=0.5;%Critical sheak stress for vegetated areas

%Organic accretion by vegetation
P.AccreteOrganic=1;
P.Korg=6/1000/365;% [mm/yr]

%ON/OFF processes
P.VEGETATION=1;%vegeation resistance and settling (DOES NOT controll organic accretion)
P.computeSeaWaves=1;
P.computeEdgeErosionSea=1;
P.computetide=1;

%Various boundary condtions
P.periodic=0;
P.imposeseaboundarydepthmorphoALL=1; %to use when a channel mouth is at a boundary

%Pond dynamics
P.calculateponddynamics=1;
    P.Epondform=0.4*10^-5;%4*10^-4/10;%probabiliy of new pond formation per year (area/area)
    P.zpondcr=-0.2;%P.Trange/4;%base of new pond formation with respect to MSL
    P.minponddepth=0.1;%1; %minimum depth to define a pond after you identified the "lakes"
    P.maxdpond=max(0.2,max(P.minponddepth*1.1,0.15*P.Trange));%0.5;%0.5;%maximum depth scour of a new pond
    P.zntwrk=(P.Trange/2)*0.5;%(P.Trange/2)*0.2;%P.Trange/2*0.9;%P.Trange/2-0.3;%0.3; %depth above msl that defines the channel network.  the smaller the harder to drain!
    P.aPEXP=0.015*10;%isolated pond expansion rate m/yr
    P.ponddeeprate=0.003;%m/yr

%Global numerical limiters
limitdeltaz=2;%[m]
limitmaxup=1;%[m]

%Time parameters
tmax1000;% number of time steps
tINT=1;%how many time steps you want to do the plot (if 1 you plot every time step). Does not affect the computation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time series
numberserie=10000;
dtOserie=ones(numberserie,1)*365;

time=cumsum(dtOserie)/365; %converted to years, just to plot. Does not affect the computation
time=[time(2:end);time(end)];

%SeaWave direction
angleWINDserie=rand(numberserie,1)*360; %every time step a random direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%Geometry Initilization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M,dx,A,z,Active,x,y,msl]=initializegeometry(P);%use this to create a new marsh
%load REFERENCEMARSHr07; %use this to load an existign marsh

%savegeometry(N,M,dx,A,z,Active,x,y,msl,'name_SAVEYOURCURRENTCONFIGURATION');%execute this once to save you current configuration

%Store value for mass balance check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumzIN=sum(z(A==1));
FLXz=zeros(4,1);
KBTOT=0;
zOX=0;
pondloss=0;

IO.z=z;
IO.msl=msl;
IO.Active=Active;
fIO.FLXz=FLXz;
fIO.pondloss=pondloss;
fIO.KBTOT=KBTOT;
fIO.zOX=zOX;

zbedo=z-msl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



makevideo=0;%put one if you want to create a video
%v=VideoWriter('nameofvideo','Motion JPEG AVI');

%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if makevideo==1;open(v);end 
s=0;step=0;tic;
for t=1:tmax;  %iteration over the tmax time stpes
    
%SeaWave direction
angleWIND=angleWINDserie(t);%rand(1)*360; %every time step a random direction
%Length of time step
dtO=dtOserie(t);

if t==1;dto=0.00001;else;dto=dtO;end   

dti=0;dt=dto;
while dti<dto;
    firstattemp=1;maxdeltaz=limitdeltaz+1;maxup=limitmaxup+1;
        while maxdeltaz>limitdeltaz | maxup>limitmaxup
        if firstattemp==1;else;dt=dt/2*min(limitdeltaz/maxdeltaz,limitmaxup/maxup);end;firstattemp=0;
        if t<=2;dt=min(0.5*365,dt);end
        [IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstep(A,P,dx,dt,IO,fIO,angleWIND,t);
        step=step+1; %this is how many time you called the function mainevolution step
        end

    %the partial updating step was succefull! Keep going
    IO=IOtemp;
    fIO=fIOtemp;
    dti=dti+dt;%how much you moved forward
    dt=min(dt*2,max(0,dto-dti));%the remaining time in the time step
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%PLOT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(t,tINT)==0;s=s+1;  
%read the variables
names = fieldnames(IO);
for i=1:length(names);eval([names{i} '=IO.' names{i} ';' ]);end

%read the fluxes
names = fieldnames(fIO);
for i=1:length(names);eval([names{i} '=fIO.' names{i} ';' ]);end
    
%read the plot
names = fieldnames(PLT);
for i=1:length(names);eval([names{i} '=PLT.' names{i} ';' ]);end

zbed=z-msl;

IM=imagesc(x,y,zbed');axis equal;set(IM,'alphadata',~(A'==0));set(gca,'YDir','normal');
cmp=demcmap([-3 P.Trange/2],256); %-3 1 P.Trange/2
colormap(cmp)
caxis([-3 P.Trange/2]);
title(strcat(num2str(time(t)),' years ',num2str(step)))




%THIS IS JUST TO MAKE THE SEDIMENT BUDGET AND CHECK THAT THE SEDIMENT
%CONSERVE SEDIMENT!!! NOT NECESSARY FOR THE COMPUTATION
QseaTide=FLXz(2);
QmouthTide=FLXz(4);
sumFLUX2=-QseaTide*dx-QmouthTide*dx;
sumz=sum(z(A==1));
%NOTE: Thsi is the equivalent volumetric flux, not the mass flux
checksum= [(sumzIN-sumz)+sumFLUX2/dx^2]-pondloss+KBTOT-zOX  ;
if abs(checksum(1))>0.1  ;checksum,pause;end% if this is not zero, the code has an error somewhere


%this is to make the video
if makevideo==1;V=getframe(figure(1));writeVideo(v,V);end


pause(0.1)
end%end of plot
end; %end of panel plotting



%this is to make the video
if makevideo==1;close(v);end
