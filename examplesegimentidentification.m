clear;close all
%I  = imread('circuit.tif');

load zBarnsmall
I=zBarnsmall;


%rotI = imrotate(I,33,'crop');

%BW = edge(rotI,'canny');
BW=I>-1.1 & I<-0.6;


figure;
subplot(1,2,1);imagesc(I)
subplot(1,2,2);imagesc(BW)

figure
[H,T,R] = hough(BW);
imshow(H,[],'XData',T,'YData',R,...
            'InitialMagnification','fit');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;

P  = houghpeaks(H,250,'threshold',ceil(0.3*max(H(:))));
a=find(P(:,2)>30 & P(:,2)<90); 
P=P(a,:);
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');

lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',20);
%figure, imshow(rotI), hold on
figure;imagesc(BW), hold on
max_len = 0;


for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

%    % Determine the endpoints of the longest line segment
%    len = norm(lines(k).point1 - lines(k).point2);
%    if ( len > max_len)
%       max_len = len;
%       xy_long = xy;
%    end
end


A=I*0;
%figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(lines)
   xy = [lines(i).point1; lines(i).point2];

xca(i,1)=xy(1,1)
yca(i,1)=xy(1,2)
xcb(i,1)=xy(2,1)
ycb(i,1)=xy(2,2)
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
%Z(sub2ind(size(A),X,Y))=-1;
end

figure;imagesc(A)

Barnditeches=A;