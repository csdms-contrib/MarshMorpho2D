function [plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=seaboundaryNeumanbedelevationALLBOUNDARY(A,plyr,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3)

[N,M]=size(A);
p=find(A==2);%exclude the NOLAND CELLS (A==0)

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N] 

[a,q]=excludeboundarycell(k,N,M,p);
a=a(A(q(a))==1);%only inclued the cells in whcih you can creep to



plyr(p(a))=plyr(q(a));
Y1(p(a))=Y1(q(a));
Y2(p(a))=Y2(q(a));
Y3(p(a))=Y3(q(a));
Yb(p(a))=Yb(q(a));
flyrb1(p(a))=flyrb1(q(a));
flyrb2(p(a))=flyrb2(q(a));
flyrb3(p(a))=flyrb3(q(a));

[px,py]=ind2sub(size(A),p(a));
[qx,qy]=ind2sub(size(A),q(a));
for i=1:length(a)
flyr1(px(i),py(i),:)=flyr1(qx(i),qy(i),:);
flyr2(px(i),py(i),:)=flyr2(qx(i),qy(i),:);
flyr3(px(i),py(i),:)=flyr3(qx(i),qy(i),:);
% plyr(px(i),py(i))=plyr(qx(i),qy(i));
% Y1(px(i),py(i))=Y1(qx(i),qy(i));
% Y2(px(i),py(i))=Y2(qx(i),qy(i));
% Y3(px(i),py(i))=Y3(qx(i),qy(i));
% Yb(px(i),py(i))=Yb(qx(i),qy(i));
% flyrb1(px(i),py(i))=flyrb1(qx(i),qy(i));
% flyrb2(px(i),py(i))=flyrb2(qx(i),qy(i));
% flyrb3(px(i),py(i))=flyrb3(qx(i),qy(i));
end


end

%z=zb-(Yb+plyr*0.3)-(Y1+Y2+Y3);
%figure;imagesc(z);colormap('jet');pause

