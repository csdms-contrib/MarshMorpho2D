function zb=setmouthbathymetrynoflux(zb,direction); %1 is north, 
%set the initial mouth no-flux bathymetry. set same zb and have the same sediment thickness!

if direction==1
%NORTH
zb(end,:)=zb(end-1,:);%tho make boundary cell equal to first cell, diomaile

elseif direction==2
%EAST
zb(:,end)=zb(:,end-1);%tho make boundary cell equal to first cell, diomaile
end
