function [S,AC]=drainpond(z,A,S,N,M,dx,dt,zntw,distdr);

AC=A*0;%the channel netwrk. 0 if not a channel netwrk  1 is a channel network

%need to add to the driange systems those channel that are created in a transitory way%%%%%
CC = bwconncomp((z>-zntw),8); %find everything that migth be a channel: deep and connected to the main channel network
for i=1:CC.NumObjects
    ind=CC.PixelIdxList{i};
    a=find(A(ind)==2);%if that blob includes a cell connect to a channel (A==2)
    if a>0;AC(ind)=1;end %all the elements become a channel
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%G=AC;
%%%%%%%%%MAKE IT FATTER BY A CERTAIN LAG
fat=floor(distdr/dx);  %how much thick to make the rim
for lag=1:fat
    p = find(AC==1);[row col]=ind2sub(size(A),p);
    for k = [N -1 1 -N]
    %avoid to the the cells out of the domain (risk to make it periodic...)
    if k==N;a=find(col+1<=M);end;
    if k==-N;a=find(col-1>0);end;
    if k==-1;a=find(row-1>0);end;
    if k==1;a=find(row+1<=N);end;   
    q=p+k;%the translated cell
    AC(q(a))=1; %make the translated cell a channel network
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CC = bwconncomp((S==1),8); %find the isolated ponds
for i=1:CC.NumObjects
    ind=CC.PixelIdxList{i};
    a=find(AC(ind)==1);%if that pond includes a cell within the channel network
    if a>0;S(ind)=0;end %all the elements in the pond become drained (S==0)
end

