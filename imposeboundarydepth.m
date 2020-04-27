function d=imposeboundarydepth(A,d,opt,dimposed);

if opt ==1
    
[m,n] = size(A);
%boundary cell impose depth
p = find(A==1 | A==21);[row col]=ind2sub(size(A),p);
for k = [m -1 1 -m]
if k==m;aa=find(col+1<=n);end;if k==-m;aa=find(col-1>0);end;if k==-1;aa=find(row-1>0);end;if k==1;aa=find(row+1<=m);end;
q=p+k;%the translated cell
a=aa(A(q(aa))==2);%only select the cells where you imposed co (A==2)
d(q(a))=d(p(a));
end


elseif opt==2
d=d;%    
%or you can imposed the depth you want
%d(A==2 | A==22)=dimposed;

end
    
    