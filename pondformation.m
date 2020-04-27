function [deltaY2,pondlost]=pondformation(A,dx,dt,Epondform,z,zpondcr,maxdpond,zsill,pondlost);

zoriginal=z;

p=find(A==1); %how many points can become a pond, how much normal cells are present.

dz=max(z(p)-zpondcr,0);%base level of the pond
dz=min(maxdpond,dz);%maximum scour depth

%a=find(dz>0);%only create new ponds when the scour can be made!
a=find(dz>0 & z(p)<zsill);%only create new ponds when the scour can be made!

r=rand(length(a),1);

aa=a(r<Epondform*dt/365*dx^2);
for i=1:length(aa)
    if dz(aa(i))>0
    %S(p(aa(i)))=1;
    z(p(aa(i)))=z(p(aa(i)))-dz(aa(i));
    pondlost=pondlost+dz(aa(i));
    end
end


deltaY=zoriginal-z;
deltaY2=deltaY;