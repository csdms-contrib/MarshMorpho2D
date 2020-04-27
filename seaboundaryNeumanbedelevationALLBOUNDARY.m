function [Y2]=seaboundaryNeumanbedelevationALLBOUNDARY(A,Y2)

[N,M]=size(A);
p=find(A==2);%exclude the NOLAND CELLS (A==0)

[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N] 

[a,q]=excludeboundarycell(k,N,M,p);
a=a(A(q(a))==1);%only inclued the cells in whcih you can creep to


Y2(p(a))=Y2(q(a));

end


