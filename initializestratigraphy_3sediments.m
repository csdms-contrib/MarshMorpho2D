function [Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy(z,N,M,P)
Yb=ones(N,M)*P.tcko;%tickness of bed layer
plyr=ones(N,M)*P.levo; %initial level occupied 

zb=z+(Yb+plyr.*(P.dlyr))+P.YUi;%substrate (hard rock) elevation

%the heigth of the statigraphy column (without Y, the mobile layer)
zs=zb-(Yb+plyr.*(P.dlyr)); 

%the hight of the mobile layer
Y=zs-z;

%initial composition of the active layer
Y1=P.initialfU*Y;
Y2=(1-P.initialfU)*Y;
Y3=0*Y;
%composition of the stratigraphy column
flyr1=NaN*zeros(N,M,P.nlyr);
flyr2=NaN*zeros(N,M,P.nlyr);
flyr3=NaN*zeros(N,M,P.nlyr);

flyr1(:,:,1:plyr)=P.initialf;
flyr2(:,:,1:plyr)=1-P.initialf;
flyr3(:,:,1:plyr)=0;

flyrb1=ones(N,M)*P.initialf;%composition of bottom layer
flyrb2=ones(N,M)*(1-P.initialf);%composition of bottom layer
flyrb3=ones(N,M)*0;%composition of bottom layer