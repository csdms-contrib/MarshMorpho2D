function [deltaY2,Pedge,Y2OX,EdgeERY2]=EdgeerosionXXX(P,z,aw,maxedgeheight,fox,dt,dx,MASK,A,Y2OX);



%MASK 1 are the mudflat
%MASK 2 is where the are no waves
%fox=0;

%the non-periodic cell are walls
MASK(1,:)=0;
MASK(end,:)=0;

%you get the wave power from the 4 sourronding cells
Pedge=[P(:,1)*0 P(:,1:end-1)]+[P(1,:)*0; P(1:end-1,:)]+[P(:,2:end) P(:,end)*0]+[P(2:end,:); P(end,:)*0];
Pedge(MASK==1)=0;%the wave power in the mudflat becomes zero!
%find the marsh cells with some wave power around it
edg=find(A==1 & MASK==0 & Pedge>0);
%rng('shuffle');
r=rand(length(edg),1);% you only need rng at the beginnign of the loop
a=find(r<aw*Pedge(edg)*dt/dx);
%edgG=edg;save edgG edgG

%these cells will erode (the "high" cells)
edger=edg(a);

deltaY2=MASK*0;
EdgeERY2=MASK*0;



 %m  max heigth that is eroded, to avoid strange erosion of channels
[N,M]=size(z);
for i=1:length(edger)
    
    %these are the adjacent cells
    [I,J] = ind2sub(size(z),edger(i));
    p=[sub2ind(size(z),mod(I+1-1,N)+1,J) sub2ind(size(z),mod(I-1-1,N)+1,J) sub2ind(size(z),I,mod(J+1-1,M)+1) sub2ind(size(z),I,mod(J-1-1,M)+1)];
 
    %standard way
    %a=find(MASK(p)==1 & A(p)==1); %only choses the adjacent cells if they are mudlfat and if they are standard cells
    
    %to arode also at the opne boundary
    a=find(MASK(p)==1); 
 
    
 %pause(1)
    if a>0 %yes, there are adjacent cells   
            dz=([z(edger(i))-z(p(a))]);
            dz=mean(max(dz,0));%dz=min(max(dz,0));            %dz=max(max(dz,0));
            dz=min(dz,maxedgeheight);
            %dz=2;%maxedgeheight;
                        
            %this is how much the bed is lowered, includes everything!!!
            deltaY2(edger(i))=deltaY2(edger(i))+dz;

            %This goes into resuspension. %This is what is conserved!!!
            EdgeERY2(edger(i))=EdgeERY2(edger(i))+dz.*(1-fox);
                                   
            %Keep track of how much you oxydeized!!!!
            Y2OX=Y2OX+dz.*fox;
    

    end
   
end


    