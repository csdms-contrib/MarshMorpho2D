function pcolorCENTER(x,y,z,dx);
[N,M]=size(x);

X=ones(N*2,M);
Y=ones(N*2,M);
Z=ones(N*2,M);

i=[1:N];
%for i=1:N
    X(1+2*(i-1),:)=x(i,:);
    X(2+2*(i-1),:)=x(i,:)+dx;
    Y(1+2*(i-1),:)=y(i,:);
    Y(2+2*(i-1),:)=y(i,:);
    Z(1+2*(i-1),:)=z(i,:);
    Z(2+2*(i-1),:)=z(i,:);
%end
pcolor(X,Y,Z)