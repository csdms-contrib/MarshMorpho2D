function [a,q]=excludeboundarycellDIAG(k,N,M,p);
[row col]=ind2sub([N M],p);
if k==N+1;a=find(col+1<=M & row+1<=N);end;
if k==-N-1;a=find(col-1>0 & row-1>0);end;
if k==-1+N;a=find(row-1>0 &col+1<=M);end;
if k==1-N;a=find(row+1<=N & col-1>0);end;

q=p+k;%the translated cell