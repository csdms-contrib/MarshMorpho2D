function [Um,TP,HS,F,kwave,PW]=SeaWaves(h,angle,hwSea_lim,range,wind,MASK,ndir,dx);

%angle=0;
Um=0*h;TP=0*h;HS=0*h;

extrafetch=1000;%[m}
Lbasin=1000/dx;
Fetchlim=max(50,dx*2);%dx*2;%600;%dx*2*10;
dlo=hwSea_lim; %minimum water depth to calculate wave. below this you don't calculate it



MASK(:,1)=0;MASK(:,end)=0;


extra=1;%if (angle<90 | angle>270);extra=1;else;extra=1;end


%The standard way
if extra==0
F=calculatefetch(MASK,ndir,dx,angle);
end
% 
% %For Georgia
%extrafetch=0;%[m}
%F=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,angle,extrafetch);
% F1=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,angle,extrafetch);
% F2=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,mod(angle+90,360),extrafetch);
% F3=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,mod(angle+180,360),extrafetch);
% F4=calculatefetchWITHEXTRASonlyOCEAN(MASK,ndir,dx,mod(angle+270,360),extrafetch);
% F=(F1+F2+F3+F4)/4;

%For the idealize basin
%extrafetch=10000;%[m}
if extra==1;
F=calculatefetchWITHEXTRAS(MASK,ndir,dx,angle,extrafetch,Lbasin,MASK);
end
% F1=calculatefetchWITHEXTRAS(MASK,ndir,dx,angle,extrafetch);
% F2=calculatefetchWITHEXTRAS(MASK,ndir,dx,mod(angle+90,360),extrafetch);
% F3=calculatefetchWITHEXTRAS(MASK,ndir,dx,mod(angle+180,360),extrafetch);
% F4=calculatefetchWITHEXTRAS(MASK,ndir,dx,mod(angle+270,360),extrafetch);
% F=(F1+F2+F3+F4)/4;

Fo=F;
            
%             %For all the modified ways. Creates a buffer on the side
%             %boundaries. Just used as a mask, the actual value is not
%             %importnat, just need to be larger than fetchlim.
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*
%             [N,M]=size(h);
%             Fo(2+floor(N*0.5):end-1,1:20)=9999;
%             Fo(2+floor(N*0.5):end-1,end-20:end)=9999;
%             %%%%%%%%%%%%%%%%%%%
%             %%%%%%%%%%%%%%%%%%%%%$%$#$&^$#^$&#^$#^$#&^$$*&^%&*%$*^%$%&*$*&%%$&*%$&%*$%&*


F(Fo<=Fetchlim)=0;

%usa questo per isolared la mudflat
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% if extra==1;
% MASK(end-Lbasin:end,:)=1;
% F(end-Lbasin:end,:)=extrafetch;
% Fo(end-Lbasin:end,:)=extrafetch;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%diffuse the fetch field    
alphadiffusefetch=0.1;   %messo 10 for the VCR wave validation 10;%0;%%%QUESTO ERA 1 FINO AD APRILE 23 2018!!!!!
F=diffusefetch(MASK,F,alphadiffusefetch,dx); 
F(Fo<=Fetchlim | MASK==0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% imagesc(F)
% colormap('jet')
% caxis([0 max(F(:))])
% pause
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=find(Fo>Fetchlim & h>dlo & F>0 & MASK==1);%h>dlo & %a=find(Fo>dx*2);%h>dlo & %a=find(h>dlo);
D=h(a);Ff=F(a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%TRCUCCOZO TO AVOID depths too small
% hbedsheatresslim=0.5;
% h(h<hbedsheatresslim)=hbedsheatresslim;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[Hs,Tp]=YeV_correction(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
[Hs,Tp]=YeV(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
HS(a)=Hs;
TP(a)=Tp;TP(TP==0)=1;


%do not diffuse in cells outside the MASK
%HS=diffusefetch(MASK,HS,alpha,dx);
%TP=diffusefetch(MASK,TP,alpha,dx);


%hlimbedshearstress=0.5;
%h=max(hlimbedshearstress,h);% to reduce the bed shear stress for very small water depth


kwave=0*h;
kk=wavek(1./TP(a),h(a));%kk=wavekFAST(1./Tp,D);
kwave(a)=kk;
kwave(kwave==0)=1;

Um=pi*HS./(TP.*sinh(kwave.*h));

cg=(2*pi./kwave./TP)*0.5.*(1+2*kwave.*h./(sinh(2*kwave.*h)));
PW=cg*1030*9.8.*HS.^2/16;

Um(MASK==0)=0;
PW(MASK==0)=0;



% h=[0.1:0.1:3.5];
% [Hs,Tp]=YeV(3000,7,h);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
% kk=wavek(1./Tp,h);%kk=wavekFAST(1./Tp,D);
% Um=pi*Hs./(Tp.*sinh(kk.*h));
% ko=0.1/1000*3;
% aw=Tp.*Um/(2*pi);
% fw=0.00251*exp(5.21*(aw/ko).^-0.19);fw(aw/ko<pi/2)=0.3;
% figure;plot(h,0.5*1000*0.015*Um.^2,h,0.5*1000*fw.*Um.^2)