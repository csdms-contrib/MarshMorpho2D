function [S,AC,DIF]=findisolatedponds(Z,A,N,M,dx,Zntw,Zsill,distdr,minponddepth);

% %Find everything that is isolated at least by the ponddepth threhsold
% Zold=Z;
%       
%         %This is to avoid pond driange on the sides!
%         ZB=A*0;ZB(:,1)=1;ZB(:,end)=1;ZB(1,:)=1;ZB(end,:)=1;
%         Z(ZB==1 & A==1)=10;
%         Z(A==0)=10;
%       
%       Z=max(Zntw,Z);%below Znet is just a mudflat
%       Z=min(Zsill,Z);%flooding sill height
%       ZZ=Z;ZZ(isnan(ZZ))=0;
%       %ZZ=imfill(ZZ,8); %THIS DOES ALL THE HEAVY LIFTING!!!       ZZ=imfill(ZZ); %THIS DOES ALL THE HEAVY LIFTING!!!
%       ZZ=imfill(ZZ,4); %THIS DOES ALL THE HEAVY LIFTING!!!       ZZ=imfill(ZZ); %THIS DOES ALL THE HEAVY LIFTING!!!
%       DIF = ZZ-Z;%the depth of pond filling
%       S=DIF>minponddepth;%what constitutes a pond
% 
% AC=(Zold<Zntw & DIF<=minponddepth);%%IF<10^-8
% 


%Find everything that is isolated at least by the ponddepth threhsold
Zold=Z;
      
        %This is to avoid pond driange on the sides!
        ZB=A*0;ZB(:,1)=1;ZB(:,end)=1;ZB(1,:)=1;ZB(end,:)=1;
        Z(ZB==1 & A==1)=10;
        Z(A==0)=10;
      
      Z=max(Zntw,Z);
      Z=min(Zsill,Z);%flooding sill height.  This is not very importnat. It is just needed for the upland.. becuase we set Znnte = Trnage/2...
      ZZ=Z;ZZ(isnan(ZZ))=0;
      
      
      
%lateral always drain. CODE: KREMLIN
%ZZ(:,1)=-10;ZZ(:,end)=-10;

%trucco per far drenare meglio
ZZ=min(ZZ,min([ZZ(:,1)*0 ZZ(:,1:end-1)],min([ZZ(1,:)*0; ZZ(1:end-1,:)],min([ZZ(:,2:end) ZZ(:,end)*0],[ZZ(2:end,:); ZZ(end,:)*0]))));%([ZZ(:,1)*0 ZZ(:,1:end-1)], [ZZ(1,:)*0; ZZ(1:end-1,:)]):%, min([ZZ(:,2:end) ZZ(:,end)*0], min([ZZ(2:end,:); ZZ(end,:)*0]))));
      


      %ZZ=imfill(ZZ,8); %THIS DOES ALL THE HEAVY LIFTING!!!       ZZ=imfill(ZZ); %THIS DOES ALL THE HEAVY LIFTING!!!
      ZZ=imfill(ZZ,4); %THIS DOES ALL THE HEAVY LIFTING!!!       ZZ=imfill(ZZ); %THIS DOES ALL THE HEAVY LIFTING!!!
      DIF = ZZ-Z;%the depth of pond filling.   oNLY USED TO DEFINE WHAT IS AN impounded area!!!
      DIF(DIF>0.001)=ZZ(DIF>0.001)-Zold(DIF>0.001);     %tHIS IS THE ACTUAL water depth in the pond
      S=DIF>minponddepth;%what constitutes a pond
      

AC=(Zold<Zntw & DIF<=minponddepth);%%IF<10^-8

%lateral never a pond. CODE: KREMLIN
%S(:,1)=0;S(:,end)=0;




% function [S,AC]=findisolatedponds(z,A,N,M,dx,zntw,distdr,minponddepth);
% 
% 
% %Find everything that is isolated at least by the ponddepth threhsold
% zold=z;
% Z=z;
%       
%         %This is to avoid pond driange on the sides!
%         ZB=A*0;ZB(:,1)=1;ZB(:,end)=1;ZB(1,:)=1;ZB(end,:)=1;
%         Z(ZB==1 & A==1)=10;
%         Z(A==0)=10;
%       
%       Z=min(zntw,Z);%flooding sill height
%       ZZ=Z;ZZ(isnan(ZZ))=0;
%       ZZ=imfill(ZZ); %THIS DOES ALL THE HEAVY LIFTING!!!
%       DIF = ZZ-Z;%the depth of pond filling
%       S=DIF>minponddepth;%what constitutes a pond
% 
% %AC=(-zold<zntw & DIF<=minponddepth);%%IF<10^-8
% AC=(zold<zntw & DIF>=minponddepth);%%IF<10^-8
% 

















% figure
% imagesc(S)
% pause
%AC=A*0;%the channel netwrk. 0 if not a channel netwrk  1 is a channel network
%Need to add to the driange systems those channel that are created as the model evolves%%%%%
% CC = bwconncomp((z>-zntw),8); %find everything that migth be a channel: deep and connected to the main channel network
% for i=1:CC.NumObjects
%     ind=CC.PixelIdxList{i};
%     a=find(A(ind)==2);%if that blob includes a cell connect to a channel (A==2)
%     if a>0;AC(ind)=1;end %all the elements become a channel
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%G=AC;
% %%%%%%%%%MAKE IT FATTER BY A CERTAIN LAG
% if distdr>0
% fat=floor(distdr/dx);  %how much thick to make the rim
% for lag=1:fat
%     p = find(AC==1);[row col]=ind2sub(size(A),p);
%     for k = [N -1 1 -N]
%     %avoid to the the cells out of the domain (risk to make it periodic...)
%     if k==N;a=find(col+1<=M);end;
%     if k==-N;a=find(col-1>0);end;
%     if k==-1;a=find(row-1>0);end;
%     if k==1;a=find(row+1<=N);end;   
%     q=p+k;%the translated cell
%     AC(q(a))=1; %make the translated cell a channel network
%     end
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CC = bwconncomp((S==1),8); %find the isolated ponds
% for i=1:CC.NumObjects
%     ind=CC.PixelIdxList{i};
%     a=find(AC(ind)==1);%if that pond includes a cell within the channel network
%     if a>0;S(ind)=0;end %all the elements in the pond become drained (S==0)
% end

