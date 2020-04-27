function F=cumsumreset(A)

a=find(A==0);
A(a) = 1-diff([0; a]);
F=cumsum(A,1);
