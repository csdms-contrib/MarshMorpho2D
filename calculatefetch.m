function F=calculatefetch(A,ndir,dx,angle)
angle=angle-180;
[N,M]=size(A);
F=zeros(N,M);
A(isnan(A))=0;  %any NaN is a boundary, that is, a 0

%Every boundary is a wall, need to do this so the fetch does not warp around!
A(1,1:end)=0;A(end,1:end)=0;A(1:end,1)=0;A(1:end,end)=0;

di=1+mod(floor(angle/360*ndir+0.5),ndir);
alfa=(di-1)/ndir*2*pi;
m=max(abs(cos(alfa)),abs(sin(alfa)));
if (di<=(ndir*1/8) | (di>ndir*3/8 & di<=ndir*5/8) | di>ndir*7/8) 
    IND=[1+mod(round([1:N]'*cos(alfa)/m)*ones(1,M),N)]+[mod(ones(N,1)*[1:M]+round([1:N]'*sin(alfa)/m)*ones(1,M)-1,M)]*N;
else
    IND=([1+mod(ones(M,1)*[1:N]+round([1:M]'*cos(alfa)/m)*ones(1,N)-1,N)]+[mod(+round([1:M]'*sin(alfa)/m)*ones(1,N),M)]*N);
end

F(IND)=cumsumreset(A(IND))/m*dx;

%at the boundary, impose the fetch of the nearby cell
F(1,1:end)=F(2,1:end);F(end,1:end)=F(end-1,1:end);F(1:end,1)=F(1:end,2);F(1:end,end)=F(1:end,end-1);

% F1=[F(:,1) F(:,1:end-1)];
% F2=[F(1,:); F(1:end-1,:)];
% F3=[F(:,2:end) F(:,end)];
% F4=[F(2:end,:); F(end,:)];
% F=(F+F1+F2+F3+F4)/5;


%NM=max(M,N);
% P1=zeros(NM,ndir);P2=zeros(NM,ndir);
% for i=1:ndir;
%     dir=(i-1)/ndir*2*pi;
%     m=max(abs(cos(dir)),abs(sin(dir)));
%     mm(i)=m;
%     P1(:,i)=-round([1:NM]*cos(dir)/m);
%     P2(:,i)=round([1:NM]*sin(dir)/m);
% end
% 

% di=1+mod(floor(angle/360*ndir+0.5),ndir);
% alfa=(di-1)/ndir*2*pi;
% m=max(abs(cos(alfa)),abs(sin(alfa)));
% if (di<=(ndir*1/8) | (di>ndir*3/8 & di<=ndir*5/8) | di>ndir*7/8) 
%     IND=[1+mod(-P1(1:N,di)*ones(1,M),N)]+[mod(ones(N,1)*[1:M]-1-P2(1:N,di)*ones(1,M),M)]*N;
% else
%     IND=([1+mod(ones(M,1)*[1:N]-1-P1(1:M,di)*ones(1,N),N)]+[mod(-P2(1:M,di)*ones(1,N),M)]*N);
% end
% 
% F(IND)=cumsumreset(A(IND),1)/mm(di)*dx;





% di=1+mod(floor(angle/360*ndir+0.5),ndir);
% alfa=(i-1)/ndir*2*pi;
% m=max(abs(cos(alfa)),abs(sin(alfa)));
% 
% if (di<=(ndir*1/8) | (di>ndir*3/8 & di<=ndir*5/8) | di>ndir*7/8) 
%     IND=[1+mod(round([1:N]'*cos(alfa)/m)*ones(1,M),N)]+[mod(ones(N,1)*[1:M]-1-round([1:N]'*sin(alfa)/m)*ones(1,M),M)]*N;
% else
%     IND=([1+mod(ones(M,1)*[1:N]-1+round([1:M]'*cos(alfa)/m)*ones(1,N),N)]+[mod(-round([1:M]'*sin(alfa)/m)*ones(1,N),M)]*N);
% end
% 
% F(IND)=cumsumreset(A(IND),1)/mm(di)*dx;


% P1=zeros(NM,ndir);P2=zeros(NM,ndir);
% for i=1:ndir;
%     dir=(i-1)/ndir*2*pi;
%     m=max(abs(cos(dir)),abs(sin(dir)));
%     mm(i)=m;
%     P1(:,i)=-round([1:NM]*cos(dir)/m);
%     P2(:,i)=round([1:NM]*sin(dir)/m);
% end

