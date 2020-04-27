function [a,q]=periodicY(k,N,M,p);
%q is the translated cell
%CAN BE OPTIMIZED

[row col]=ind2sub([N M],[1:N*M]');

%row = rem([1:N*M]'-1,N)+1;
%row=[1:N]'*ones(1,M);row=reshape(row,1,N*M)';
%col = ([1:N*M]'-row)/N + 1;

if k==-1;
    q=p+k;
    a=find(row(p)>1);
end;

if k==1;
    q=p+k;
    a=find(row(p)<N);
end;

if k==N;a=find(p>0);
    G=[1:N*M]';
    Q=G+k;
    s1=find(col==M);
    s2=sub2ind([N M],row(s1),1+0*col(s1));
    %s2 = row(s1) + (1+0*col(s1)-1)*N;
    Q(s1)=G(s2);
    q=Q(p);
end;

if k==-N;a=find(p>0);
    G=[1:N*M]';
    Q=G+k;
    s1=find(col==1);
    s2=sub2ind([N M],row(s1),M+0*col(s1));
    %s2 = row(s1) + (M+0*col(s1)-1)*N;
    Q(s1)=G(s2);
    q=Q(p);
end;




%%%%%%%%%%%%%%%ORIGINAL
% function [a,q]=periodicY(k,N,M,p);
% [row col]=ind2sub([N M],p);
% 
% q=p+k;%the translated cell
% 
% 
% if k==-1;a=find(row>1);end;
% if k==1;a=find(row<N);end;
% 
% if k==N;a=find(p>0);
%     s1=find(col==M);
%     s2=find(col==1);
%     q(s1)=p(s2);
% end;
% if k==-N;a=find(p>0);s1=find(col==1);s2=find(col==M);q(s1)=p(s2);end;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% if k==N;a=find(p>0);
%     s1=find(col==M);
%     s2=find(col==1);
%     %q(s1)=p(s1-M*N+M+N-1);
% s2-(s1-length(p)+N-1)
%     q(s1)=p(s1-length(p)+N-1);
% end;
% pause
% if k==-N;a=find(p>0);s1=find(col==1);s2=find(col==M);q(s1)=p(s2);end;

%if k==N;a=find(p>0);s=find(col==M);col(s)=1;q=sub2ind(N,M,row,col);end;
%if k==-N;a=find(p>0);s1=find(col==1);s2=find(col==M);q(s1)=q(s2);end;
%if k==-N;a=find(col-1>0);end;
