clear; close all;clc
%NOTES
%Qs is the total lateral transport, E is the vertical erosion flux

%lateral boundary options
% 0 is closed (no flux if nothign specified. no-gradient if AW equal to 1 or -1)
% 1 is periodic

%options for lateral wave condtions


%A
%0: not in the domain
%1: a normal cell
%2: the open sea boundary
%10: the river boundary
%%%FALSE%3 and -3: no-gradient boundaries (for the lateral). 3 is left (1)  -3 is right (end)

%AP= 0 is not a pond  1 if an isolate dpond


%Various initiliaztion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first line is bottom %second line is top
% %this is the good one
 cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmap.cpt', cdata, ' ');

 cdata = [-1    50  255 255     0     255    255    255
          0    255 255 255     1     0    0    0];
 cdata = [-1    205  165 0     0     0    0    0
          0    255 255 255     1     0    0    0];
%this is the original one. boring
% cdata = [-1   0   0  255   0   255   0   0
%          0   255 255 255   1   0     0   0];
dlmwrite('mycmapCLR.cpt', cdata, ' ');

cdata = [-1    205  165 0     0     0    0    0];
dlmwrite('mycmap1.cpt', cdata, ' ');

%to get always the same random numbers
rng(2)


%%%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%
P=struct;
P.g=9.81; %gravity [m/s2]
P.rho=1030; %water density [kg/m3] 
P.rhos=2650; %sediment density (quatz-mica) [kg/m3]
P.ss=(P.rhos-P.rho)/P.rho; %relative density
P.kro=0.02; % minimum water depth [m]
P.DiffS=0.1; %coefficient for tidal dispersion [-]

P.RSLR=2.6/1000/365;  %from mm/yr to m/day (the time unit is the day!)

%tide
P.Ttide=12.5/24; %[tidal period [day]
P.Trange=3.1;%1.5; %tidal Trange [m]
P.TrangeVEG=3.1;%2.7;%1.5; %tidal Trange [m]

%storm surge
P.alpha_surge=0.25;%how much the surge contributes to tidal prsim

%swell waves
P.gridDIR=1; %1: waves propoagation is top to bottom;   -1: waves propoagation is bottom to top
P.Pswelldir=0.5;  %if 0.5, then is symmetric  (1 is left or right) (0 is rigth to left)
P.Pswellhighangle=0.5; %if zero, only low angle waves
P.Ho=4; %boundary swell height Hs [m]
P.Tp_swell=8;% %boundary swell period Tp [m]

%wind for sea waves
P.wind=10;
P.aw=0.3/365; %wave edge erodability m/yr/W/m2

P.hwSwell_lim=0; %limiter water depth for swell
P.hwSea_lim=0.5; %limiter water deth of sea waves

%these will impose the sediment discharge input at the river mouth
P.Qmouth=5*1; %river discharge per unit of cell [m2/s]
P.hmouth=15; %water depth [m]
%P.Umouth=0;%P.Qmouth/P.hmouth;  % -->  velocity -->> Qs
%P.co2mouth=500/1000; %SSC of mud at the river [g/l]

%SSC at the sea boundary
P.co1=0/1000; % Sea boundary SSC for sand [g/l]
P.co2=5/1000; %Sea boundary SSC for mud [g/l]
P.co3=0/1000; %Sea boundary SSC for mud [g/l]

P.d50_1=0.25/1000;% %sand grain size [m]
P.d50_2=0.02/1000/1000;%mud grain size [m]
P.ws1=0.02;%sand with D50=500um  0.05 %m/s
P.ws2=0.2/1000;%
P.por1=0.4;P.rbulk1=P.rhos*(1-P.por1);
P.por2=0.7;P.rbulk2=P.rhos*(1-P.por2);

%mud parameters
P.me=0.1*10^-4*24*3600;  %per day!!!
P.taucr=0.2;
P.crMUD=2/365;
P.crMARSH=0.2/365;
P.DoMUD=0*100;%base diffusivity os suspedned mud [m2/s]

P.fUpeak=pi/2;
P.fMFtide=NaN;%1/(pi/2);%0.2
P.fMFswell=0.04;%0.1; %to scale waves and tidal transport. Waves do not occur all the time
P.fMFswellMUD=0.04;%0.1;
P.fMFsea=0.3;
P.fMFriver=30/365;

%Vegetation parameters
P.Bmax=2000;
P.dBlo=-(0.237*P.TrangeVEG-0.092);
P.dBup=P.TrangeVEG/2;%-0.2;
%dlo=0.737*range-0.092;
P.alphaSAND=4; %coefficient for bedload downslope of sand
P.Cb=0.02;%9.8/45^2;

Bpeak=2.0;
nuGp=0.0138;%1/day%rate at which organic matter is stored. to convert to AMC, total mass of organic accumalted
chiref=0.158;%refractory fraction
rorg=1200;%density of organic matter * 1 - water content
Worg=0.1;
org=(Bpeak*365/2*nuGp)*chiref/(rorg*Worg); %in m/hour
P.Korg=8/1000/365;%5/1000/365;
%P.Korg=org/365;%5/1000/365;

P.fox=0;

%ON/OFF
P.VEGETATION=1;%vegeation resistance and settling
P.computemud=1;
P.computesand=0;
P.computeEdgeErosion=0;
P.computeSwellwave=0;
P.computeSeaWaves=0;
P.computetide=1;
P.computeriver=0;
P.riverwaterlevel=0;

P.periodic=0;

P.imposeseaboundarydepthmorphoNORTH=0;
P.imposeseaboundarydepthmorphoEAST=1;
P.evolvestratigraphy=0;
P.VEGstratigraphy=0;

P.calculateponddynamics=0;

%stratigraphy
P.nlyr=20; %max number of layers
P.dlyr=0.3; %thickenss of layers
P.tlyrU=0.5; %max depth to add layer %must be larger than dlyr
P.tlyrD=0.1; %min depth merge layers %mus be larger than dlyr
P.tcko=10;%tickness of bed layer
P.levo=19;%intial level occupied
P.YUi=100;%initial thickess of active layer
P.initialfU=0;%initial composition of the active layer (sand content) 
P.initialf=0;%initial composition of all the layers including bottom (sand content)

%numerical limiters
limitdeltaz=5;
limitmaxup=5;

tmax=55;%501;%in days
tINT=1;

dtO=2*365;
time=[0:tmax-1]*dtO/365;

%time series
surge=0*wblrnd(0.1,0.5,tmax,1);  surge(surge>5)=5; %surge(1:10)=0;surge(surge<1)=0;

makevideo=0;





%%%%%%%%%%%%%%%Geometry Initilization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,z,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell]=initializegeometry_3sediments(P);

%set the initial mouth no-flux bathymetry. set same zb and have the same sediment thickness!
zb=setmouthbathymetrynoflux(zb,2); %1 is north, 2 is east 

S=0*A;  %the status of a cell; (for ponds dynamics)

%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'WESTR1co5_eq');
%load WESTR1co5_eq




%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'S');
%load S

%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'WESTcr03_5Korg66');
%load WESTcr03_5Korg66


%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'sineWESTSLR1co20dx2mr31_org66wsB05crMUD3Cvt01e05');
%load sineWESTSLR1co20dx2mr31_org66wsB05crMUD3



%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'PIESLR1co20dx2mr31_org66velpi_piwsB05');
%load PIESLR1co20dx2mr31_org66velpi_piwsB05

%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'PIESLR1co20dx2mr31_org66velpi_piwsB05');
%load PIEco10SLR1dx2m_org5vel02

%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'PIEco10SLR1dx2m_mud05');
%load PIEco10SLR1dx2m_mud05%PIEco10SLR1mud05_1938
% for i=1:length(A(1,:))
%     a=find(A(:,i)>0);
%     A(1:a(1)-1,i)=999;
%     A(a(end)+1:end,i)=999;
% end
% a=find(A==0);
% A(a)=1;zb(a)=100;
% Y2(a)=100-10-P.dlyr*P.levo+2+1;
% A(A==999)=0;


%savegeometry(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,flyr,flyrb,x,y,msl,SPCLcell,'PIEco10SLR1');
%savegeometry(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,flyr,flyrb,x,y,msl,SPCLcell,'PIEco10SLR1');

%S(20:52,30:75)=1;
%[N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,flyr,flyrb,x,y,msl,SPCLcell]=initializegeometry(P);
%savegeometry(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,flyr,flyrb,x,y,msl,SPCLcell,'BTrange15slr0');
%savegeometry_3sediments(N,M,dx,A,AW,Yb,plyr,zb,z,Y1,Y2,Y3,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,x,y,msl,SPCLcell,'BTrange15slr0');
%load BTrange15slr0

%river
% mouthW=6; %2
% A(1:10,1:M/2-1-mouthW)=0;%create a concreate wall on the sides
% A(1:10,M/2+1+mouthW:end)=0;%create a concreate wall on the sides
% A(1,M/2-mouthW:M/2+mouthW)=10;
% %these are the cells in front of the river mouth
% S=A*0;S(2,M/2-mouthW:M/2+mouthW)=1;rivermouthfront=find(S==1);clear S;
% SPCLcell=struct;
% SPCLcell.rivermouthfront=rivermouthfront;
% %river
% %z(A==10)=P.hmouth;
% z(1:52,M/2-mouthW:M/2+mouthW)=z(1:52,M/2-mouthW:M/2+mouthW)+P.hmouth;
% [Yb,Y1,Y2,zb,zs,plyr,flyr,flyrb]=initializestratigraphy(z,N,M,P);


% load E:\VCR\zr
% dx=200;
% %figure;imagesc(zr);
% z=-zr(end:-1:1,:);
% A=0*z+1;A(isnan(z))=0;
% A(1,:)=2;
% AW=A*0;
% AW(:,1)=1;
% AW(:,end)=-1;
% [N,M]=size(A);
% x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;
% SPCLcell.rivermouthfront=[];
% [Yb,Y1,Y2,zb,zs,plyr,flyr,flyrb]=initializestratigraphy(z,N,M,P);
% msl=0;


%  a=find(z<2);
%  Y2(a)=Y2(a)+Y1(a);Y1(a)=0;
% for i=1:P.levo;
%     tmp=squeeze(flyr(:,:,i));
%     tmp(a)=0;
%     flyr(:,:,i)=tmp;
% end

%Store value formass balance check
sumY1IN=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1IN=sum(sumY1IN(A==1));
sumY2IN=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2IN=sum(sumY2IN(A==1));
sumY3IN=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3IN=sum(sumY3IN(A==1));
FLX1=zeros(4,1);FLX2=zeros(4,1);FLX3=zeros(4,1);KBTOT=0;Y2OX=0;
FQsW_L=0;FQsW_R=0;
pondloss=0;

IO.S=S;
IO.Y1=Y1;IO.Y2=Y2;IO.Y3=Y3;
IO.flyr1=flyr1;IO.flyr2=flyr2;IO.flyr3=flyr3;
IO.flyrb1=flyrb1;IO.flyrb2=flyrb2;IO.flyrb3=flyrb3;
IO.plyr=plyr;IO.Yb=Yb;IO.msl=msl;
fIO.FLX1=FLX1;fIO.FLX2=FLX2;fIO.FLX3=FLX3;
fIO.pondloss=pondloss;
fIO.KBTOT=KBTOT;fIO.Y2OX=Y2OX;
fIO.FQsW_R=FQsW_R;fIO.FQsW_L=FQsW_L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


z=zb-(Yb+plyr*P.dlyr)-(Y1+Y2+Y3);
Y=Y1+Y2+Y3;
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
flyrU1=max(0,Y1)./Ytot;flyrU1(Ytot==0)=1;
flyrU2=max(0,Y2)./Ytot;flyrU2(Ytot==0)=1;
flyrU3=max(0,Y3)./Ytot;flyrU3(Ytot==0)=1;
zbedo=z+msl;




%%%MAIN LOOP
%figure('units','normalized','outerposition',[0 0 1 1])
%set(gca,'Color','k')
%Seac20RSLR5mm _fillSeamudco0RSL0Riverco200Q2x5h5_1months   _RIVERU1SAND&MUD
%_FillRiverW3Q5H5Mud200SLR0   
if makevideo==1;v=VideoWriter('BTrange15slr0_MUD_ORG','Motion JPEG AVI');open(v);end %_fillSeamudco20RSL0 _fillSeamudco20RSL2
s=0;step=0;tic;
for t=1:tmax;  
    
%swell wave direction
randdir=rand(1);if randdir>P.Pswelldir;dirsign=1;else;dirsign=-1;end %dirsign=sign((mod(t,2)-0.5));
rndhl=rand(1);if rndhl>P.Pswellhighangle;angleswell=dirsign*(rand(1)*45);else;angleswell=dirsign*min(80,(rand(1)*45+45));end
%angleswell=-70;%

%SeaWave direction
angleWIND=rand(1)*360; %every time step a random direction

if mod(t,2)==0;P.Qmouth=5;P.co2mouth=000/1000;else;P.Qmouth=5;P.co2mouth=0/1000;end

if t==1;dto=0.00001;else;dto=dtO;end   
%dto=10^-20;
dti=0;dt=dto;
while dti<dto;
    firstattemp=1;maxdeltaz=limitdeltaz+1;maxup=limitmaxup+1;
        while maxdeltaz>limitdeltaz | maxup>limitmaxup
        if firstattemp==1;else;dt=dt/2*min(limitdeltaz/maxdeltaz,limitmaxup/maxup);end;firstattemp=0;%
        if t<=1;dt=min(0.2*365,dt);end
        [IOtemp,fIOtemp,maxdeltaz,maxup,PLT]=mainevolutionstep(A,AW,SPCLcell,P,dx,dt,zb,IO,fIO,surge(t),t,angleswell,angleWIND);
        step=step+1;
        end

    %the partial updating step was succefull
    IO=IOtemp;
    fIO=fIOtemp;
    dti=dti+dt;%how much you moved forward
    dt=min(dt*2,max(0,dto-dti));%the remaining time in the time step
end

%if time(t)>100;P.co2=200/1000;end
%if time(t)>2000;P.RSLR=10/1000/365;end


%%%PLOT
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


z=zb-(Yb+plyr*P.dlyr)-(Y1+Y2+Y3);
Y=Y1+Y2+Y3;
Ytot=(max(0,Y1)+max(0,Y2)+max(0,Y3));
flyrU1=max(0,Y1)./Ytot;flyrU1(Ytot==0)=1;
flyrU2=max(0,Y2)./Ytot;flyrU2(Ytot==0)=0;
flyrU3=max(0,Y3)./Ytot;flyrU3(Ytot==0)=0;
zbed=z+msl;





ax1 = subplot(1,3,1); %set(IM,'alphadata',~(A==0));
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(y,x,-zbed);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-3 P.Trange/2],256); %-3 1
colormap(ax1,cmp)
caxis([-3 P.Trange/2]);
%colorbar('hori') 
 %colorbar
 
 ax1 = subplot(1,3,2); %set(IM,'alphadata',~(A==0));
%IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
IM=imagesc(y,x,-zbedo);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
cmp=demcmap([-3 P.Trange/2],256); %-3 1
colormap(ax1,cmp)
caxis([-3 P.Trange/2]);
%colorbar('hori') 
 %colorbar
 
 
%    ax2 = subplot(1,3,3);
%   IM=imagesc(y,x,U);set(IM,'alphadata',~(A==0));axis equal;set(gca,'YDir','normal');%colormap('jet');%colorbar('hori') 
%   caxis([0 0.3]);colormap('jet')
%  colormap(ax2,flipud(hot))
 

 max(1000*SSC(:))
  ax2 = subplot(1,3,3);
  IM=imagesc(y,x,1000*SSC);set(IM,'alphadata',~(A==0));axis equal;set(gca,'YDir','normal');%colormap('jet');%colorbar('hori') 
  caxis([0 10]);colormap('jet')
 colormap(ax2,flipud(hot))
 %colorbar
 

%   ax2 = subplot(2,3,2);
%  IM=imagesc(y,x,h);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
% % colormap(ax3,flipud(gray))
%  %first line is bottom %second line is top
% %cptcmap('mycmap1', 'mapping', 'direct'); 
% colormap('jet')
%  caxis([0 10])

% ax3 = subplot(3,1,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,h);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% cmp=demcmap([-20 3]);
% colormap(ax3,cmp)
% caxis([-20 3]);
% %colorbar('hori') 

%  ax3 = subplot(2,3,3);
%  IM=imagesc(y,x,hriver);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
% % colormap(ax3,flipud(gray))
%  %first line is bottom %second line is top
% %cptcmap('mycmap1', 'mapping', 'direct'); 
% colormap('jet')
%  caxis([0 3])

%  ax3 = subplot(2,3,3);
%  IM=imagesc(y,x,1-flyrU1);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
% % colormap(ax3,flipud(gray))
%  %first line is bottom %second line is top
% cptcmap('mycmap1', 'mapping', 'direct'); 
%  caxis([0 1])
% % 

%S(AC==1)=2;
% ax3 = subplot(2,3,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,U);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% caxis([0 0.5]);%colormap('jet')
%  %colorbar
 
%  ax4 = subplot(2,3,4); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,AC);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% caxis([0 1]);%colormap('jet')
%  %colorbar

% ax1 = subplot(2,2,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,PLT.Ux);axis equal;set(IM,'alphadata',~(A==0));%set(gca,'YDir','normal');%colormap('jet');
% 
% ax1 = subplot(2,2,4); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,PLT.Uy);axis equal;set(IM,'alphadata',~(A==0));%set(gca,'YDir','normal');%colormap('jet');



%  ax3 = subplot(2,3,3);
%  IM=imagesc(y,x,Y2);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');colormap('jet');%colorbar('hori') 
%  colormap(ax3,flipud(gray))
%  caxis([0 2])
 
%  ax1 = subplot(2,3,3); %set(IM,'alphadata',~(A==0));
% %IM=imagesc(x,y,-z'-msl);set(IM,'alphadata',~(A'==0));axis equal;set(gca,'YDir','normal');xlim([0 x(end)]);colormap('jet');caxis([-20 4]);%colorbar('hori') 
% IM=imagesc(y,x,Hs);axis equal;set(IM,'alphadata',~(A==0));set(gca,'YDir','normal');%colormap('jet');
% caxis([0 4])
% colorbar
 




% subplot(3,1,3);%MUD vs total
% %transx_y=1;Mi=20; %cross section along direction x(1) or y(2)
% transx_y=2;Mi=N; %cross section along direction x(1) or y(2)
% [si,zi,Ai]=getstrat2plot(-zb,flyrU2,flyrb2,flyr2,P.nlyr,P.dlyr,plyr,Y,N,M,Mi,Yb,dx/1000,transx_y);
% hold off;pcolorCENTER(si-dx/2/1000,zi-msl,-Ai,dx/1000);%axis equal;set(gca,'YDir','normal');
% cptcmap('mycmap', 'mapping', 'direct','ncol',256); 
% caxis([-1 1.001]);shading flat
% %colorbar  
% hold on;plot(si,si*0+msl*0+P.Trange/2+surge(t),'-c',si,si*0+msl*NaN,'-k',si,si*0+0*msl+P.Trange/2,'--k',si,si*0+0*msl-P.Trange/2,'--k')
% ylim([-5 2]);%caxis([0 1]) si,-zbed(Mi,:)*NaN,'.-b',
% 

% subplot(4,1,3);%MUD vs total
% transx_y=1;Mi=20; %cross section along direction x(1) or y(2)
% [si,zi,Ai]=getstrat2plot(-zb,flyrU1+flyrU3,flyrb1+flyrb3,flyr1+flyr3,P.nlyr,P.dlyr,plyr,Y,N,M,Mi,Yb,dx/1000,transx_y);
% hold off;pcolorCENTER(si-dx/2/1000,zi,-Ai,dx/1000);%axis equal;set(gca,'YDir','normal');
% cptcmap('mycmap', 'mapping', 'direct','ncol',256); 
% caxis([-1 1.01]);shading flat
% %colorbar
% hold on;plot(si,si*0+msl+P.Trange/2+surge(t),'-c',si,si*0+msl*NaN,'-k',si,si*0+msl+P.Trange/2,'--k',si,si*0+msl-P.Trange/2,'--k')
% ylim([-15 3]);%caxis([0 1])

% subplot(4,1,4);%ORGANIC vs total
% transx_y=1;Mi=20; %cross section along direction x(1) or y(2)
% [si,zi,Ai]=getstrat2plot(-zb,flyrU1+flyrU2,flyrb1+flyrb2,flyr1+flyr2,P.nlyr,P.dlyr,plyr,Y,N,M,Mi,Yb,dx/1000,transx_y);
% hold off;pcolorCENTER(si-dx/2/1000,zi,-Ai,dx/1000);%axis equal;set(gca,'YDir','normal');
% cptcmap('mycmapCLR', 'mapping', 'direct','ncol',256); 
% caxis([-1 1.01]);shading flat
% hold on;plot(si,si*0+msl+P.Trange/2+surge(t),'-c',si,si*0+msl*NaN,'-k',si,si*0+msl+P.Trange/2,'--k',si,si*0+msl-P.Trange/2,'--k')
% ylim([-15 3]);%caxis([0 1])


title(strcat(num2str(time(t)),' years ',num2str(step)))




%TERM1; Qouthriver: if postive it enters
%TERM2; Qseatide: if postive it exits
%TERM3; Qseariver: if postive it exits
%TERM4 Qmouth tide. THIS IS IMPOSED ZERO BY setting D=0 at the mouth in sedtran

%SAND
QmouthRiver=FLX1(1);QseaTide=FLX1(2);QseaRiver=FLX1(3);QmouthTide=FLX1(4);
sumFLUX1=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%MUD
QmouthRiver=FLX2(1);QseaTide=FLX2(2);QseaRiver=FLX2(3);QmouthTide=FLX2(4);
sumFLUX2=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;
%ORG
QmouthRiver=FLX3(1);QseaTide=FLX3(2);QseaRiver=FLX3(3);QmouthTide=FLX3(4);
sumFLUX3=dx*QmouthRiver-QseaTide*dx-QseaRiver*dx-QmouthTide*dx;

sumY1=sumSedcolum(Yb,flyrb1,flyr1,P.dlyr,Y1);sumY1=sum(sumY1(A==1));
sumY2=sumSedcolum(Yb,flyrb2,flyr2,P.dlyr,Y2);sumY2=sum(sumY2(A==1));
sumY3=sumSedcolum(Yb,flyrb3,flyr3,P.dlyr,Y3);sumY3=sum(sumY3(A==1));

%NOTE: Thsi is the equivalent volumetric flux, not the mass flux
checksum=[[(sumY1IN-sumY1)+sumFLUX1/dx^2]+fIO.FQsW_L+fIO.FQsW_R  [(sumY2IN-sumY2)+sumFLUX2/dx^2]+pondloss  [(sumY3IN-sumY3)-Y2OX+sumFLUX3/dx^2+KBTOT]]
%[fIO.FQsW_R fIO.FQsW_R]
%if abs(checksum(1))>1;pause;end
%RfIO.FQsW_L
%fIO.FQsW_R
%+fIO.FQsW_L+fIO.FQsW_R

if makevideo==1;V=getframe(figure(1));writeVideo(v,V);end
pause(0.1)
end

end

if makevideo==1;close(v);end %UNCOMMENT THIS TO CREATE A VIDEO



