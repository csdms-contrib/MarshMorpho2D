clear;close all

load zSmall
z=zSmall;A=z*0;
figure('units','normalized','outerposition',[0 0 1 1])
for i=1:200
imagesc(z);caxis([-1 0]);xlim([0 1200]);colormap('jet')%ylim([1600 3700])
pause
[xca(i,1),yca(i,1)]=ginput(1);
[xcb(i,1),ycb(i,1)]=ginput(1);
a=[xca(i,1),yca(i,1)];b=[xcb(i,1),ycb(i,1)];
ra=abs(a(2)-a(1));rb=abs(b(2)-b(1));
x=[a(2) b(2)];y=[a(1) b(1)];
if ra<rb;
X=a(2):0.2*sign(b(2)-a(2)):b(2);Y=interp1(x,y,X);
else
Y=a(1):0.2*sign(b(1)-a(1)):b(1);X=interp1(y,x,Y);
end
X=round(X);Y=round(Y);
A(sub2ind(size(A),X,Y))=1;  
z(sub2ind(size(A),X,Y))=2;
end

xyditches=find(A==1);