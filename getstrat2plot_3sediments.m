function [si,zi,Ai]=getstrat2plot(zb,flyrU1,flyrU2,flyrU3,flyrb1,flyrb2,flyrb3,flyr1,flyr2,flyr3,nlyr,dlyr,plyr,Y,N,M,Mi,Yb,dx,transx_y);

% if transx_y==1
% Q=squeeze(isfinite(flyr(:,Mi,:)));Q=cumsum(Q,2);
% zi=[zb(:,Mi)  zb(:,Mi)+Yb(:,Mi)   zb(:,Mi)*ones(1,nlyr)+Yb(:,Mi)*ones(1,nlyr)+Q*dlyr   zb(:,Mi)+Yb(:,Mi)+plyr(:,Mi)*dlyr+Y(:,Mi)];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
% si=(ones(nlyr+3,1)*[0:N-1]*dx)';
% Ai=[flyrb(:,Mi)  squeeze(flyr(:,Mi,:))   flyrU(:,Mi)   flyrU(:,Mi)*0];
% 
% elseif transx_y==2;
% Q=squeeze(isfinite(flyr(Mi,:,:)));Q=cumsum(Q,2);
% zi=[0*Yb(Mi,:)'   Yb(Mi,:)'   Yb(Mi,:)'*ones(1,nlyr)+Q*dlyr   Yb(Mi,:)'+plyr(Mi,:)'*dlyr+Y(Mi,:)'];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
% si=(ones(nlyr+3,1)*[0:M-1]*dx)';
% Ai=[flyrb(Mi,:)'  squeeze(flyr(Mi,:,:))   flyrU(Mi,:)'   flyrU(Mi,:)'*0];
% 
% end

%Y(:,Mi)
%pause
%with the white on top
if transx_y==1
Q=squeeze(isfinite(flyr1(:,Mi,:)));Q=cumsum(Q,2);
zi=[zb(:,Mi)  zb(:,Mi)+Yb(:,Mi)   zb(:,Mi)*ones(1,nlyr)+Yb(:,Mi)*ones(1,nlyr)+Q*dlyr   zb(:,Mi)+Yb(:,Mi)+plyr(:,Mi)*dlyr+Y(:,Mi) zb(:,Mi)+Yb(:,Mi)+plyr(:,Mi)*dlyr+max(0,Y(:,Mi))];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
si=(ones(nlyr+4,1)*[0:N-1]*dx)';
Ai=[flyrb1(:,Mi)  squeeze(flyr1(:,Mi,:))   flyrU1(:,Mi)   flyrU1(:,Mi)*0-0.1 flyrU1(:,Mi)*0-0.1];

elseif transx_y==2;
    %%%%% I THINK SOMETHING IS MISSING HERE
Q=squeeze(isfinite(flyr1(Mi,:,:)));Q=cumsum(Q,2);
zi=[0*Yb(Mi,:)'   Yb(Mi,:)'   Yb(Mi,:)'*ones(1,nlyr)+Q*dlyr   Yb(Mi,:)'+plyr(Mi,:)'*dlyr+Y(Mi,:)'];%*ones(1,5);zi(:,2:3)=zi(:,2:3)+2.5;
si=(ones(nlyr+3,1)*[0:M-1]*dx)';
Ai=[flyrb(Mi,:)'  squeeze(flyr1(Mi,:,:))   flyrU1(Mi,:)'   flyrU1(Mi,:)'*0];

end