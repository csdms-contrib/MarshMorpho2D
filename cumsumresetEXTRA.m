function F=cumsumreset(A,extrafetch)

% a=find(A==0);
% A(a) = 1-diff([0; a]);
% F=cumsum(A,1);

S=A;
S(1,2:end-1)=extrafetch;
S(S==0)=NaN;
S(S==1)=0;
G=cumsum(S,1);
G(isnan(G))=0;


a=find(A==0);
A(a) = 1-diff([0; a]);
F=cumsum(A,1);

F=F+G;
% % S=S(2:end,2:end-1);
% % [x,y]=find(S==0);
% % 
% % figure;plot(x,y,'.')
% % figure;plot(x)
% figure;imagesc(S);
% figure;imagesc(G);
% pause