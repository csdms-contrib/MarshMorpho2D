function [a,q]=excludeboundarycell(k,N,M,p);
[row col]=ind2sub([N M],p);
if k==N;a=find(col+1<=M);end;
if k==-N;a=find(col-1>0);end;
if k==-1;a=find(row-1>0);end;
if k==1;a=find(row+1<=N);end;

q=p+k;%the translated cell