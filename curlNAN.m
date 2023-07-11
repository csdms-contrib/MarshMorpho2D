function c=curlNAN(Ux,Uy);
%c=(circshift(Ux,[1 0])-circshift(Ux,[-1 0]))/2-(circshift(Uy,[0 1])-circshift(Uy,[0 -1]))/2;
%c=([Ux(2:end) Ux(end)]-circshift(Ux,[-1 0]))/2-(circshift(Uy,[0 1])-circshift(Uy,[0 -1]))/2;

u1=circshift(Ux,[1 0])-Ux;
u2=Ux-circshift(Ux,[-1 0]);
v1=circshift(Uy,[0 1])-Uy;
v2=Uy-circshift(Uy,[0 -1]);

% u1=[Ux(1,:)*NaN; Ux(1:end-1,:)]-Ux;
% u2=Ux-[Ux(2:end,:); Ux(end,:)*NaN];
% v1=[Ux(:,1)*NaN Ux(:,1:end-1)]-Uy;
% v2=Uy-[Ux(:,2:end) Ux(:,end)*NaN];


%c=(u1+u2)/2-(v1+v2)/2;


%old version, has more NaNs
% %c=Ux*0;
% %c=sum(u1,u2,'omitnan')/2+sum(v1,v2,'omitnan')/2;
% tmp1 = cat(3,u1,u2);
% tmp2 = cat(3,v1,v2);
% c = mean(tmp1,3,'omitnan')-mean(tmp2,3,'omitnan');
% %c=(u1+u2)/2-(v1+v2)/2;
% 

%%new version
%c=Ux*0;
%c=sum(u1,u2,'omitnan')/2+sum(v1,v2,'omitnan')/2;
tmp1 = cat(3,u1,u2);
tmp2 = cat(3,v1,v2);
c1=mean(tmp1,3,'omitnan');c1(isnan(c1))=0;
c2=mean(tmp2,3,'omitnan');c2(isnan(c2))=0;
c=c1-c2;
c(isnan(Ux))=NaN;

c=c/2;
% c(:,end)=0;
% c(:,1)=0;
% c(end,:)=0;
% c(1,:)=0;

c(1,:)=NaN;
c(end,:)=NaN;
c(:,1)=NaN;
c(:,end)=NaN;