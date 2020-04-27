function F=calculatefetchWITHEXTRAS(A,ndir,dx,angle,extrafetch,Lbasin,h)
%angle=315-180-0.
angleo=angle;
%angle=mod(angle-180,360);
angle=angle-180;
[N,M]=size(A);
F=zeros(N,M);
A(isnan(A))=0;  %any NaN is a boundary, that is, a 0c

%Every boundary is a wall, need to do this so the fetch does not warp around!
A(1,1:end)=0;A(end,1:end)=0;A(1:end,1)=0;A(1:end,end)=0;


%di=1+mod(floor(angle/360*ndir+0.5),ndir);
di=1+mod(floor(angle/360*ndir),ndir);%modified Dec 2019
alfa=(di-1)/ndir*2*pi;
m=max(abs(cos(alfa)),abs(sin(alfa)));
if (di<=(ndir*1/8) | (di>ndir*3/8 & di<=ndir*5/8) | di>ndir*7/8) 
    IND=[1+mod(round([1:N]'*cos(alfa)/m)*ones(1,M),N)]+[mod(ones(N,1)*[1:M]+round([1:N]'*sin(alfa)/m)*ones(1,M)-1,M)]*N;
else
    IND=([1+mod(ones(M,1)*[1:N]+round([1:M]'*cos(alfa)/m)*ones(1,N)-1,N)]+[mod(+round([1:M]'*sin(alfa)/m)*ones(1,N),M)]*N);
end

%F(IND)=cumsumreset(A(IND),1)/m*dx;

%error
if angleo<0 | angleo>360;A='error' ;end

% %if (angleo<44 | angleo>316)
% if (angleo<134 | angleo>226)
% F(IND)=cumsumresetEXTRA(A(IND),extrafetch/dx)/m*dx;
% else
% F(IND)=cumsumreset(A(IND))/m*dx;
% end

padding=Lbasin*2;%must be larger than fetchlim!!!
if (angleo<45 | angleo>315)
F(IND)=cumsumresetEXTRA(A(IND),extrafetch/dx)/m*dx;
        if angleo<44;
        %ll=length(F(2:end-1-floor(N*0.5),end));    
        %F(2+floor(N*0.5):end-1,end-padding:end)=2*extrafetch*[1:ll]'/ll*ones(padding+1,1)';
        %F(end,end-padding:end)=extrafetch; 
        a=find(h(:,end)<0);if a>0;Lside=a(end)-1;else;Lside=1;end
        F(Lside:end,end-100:end)=max(extrafetch,F(Lside:end,end-100:end)); %TOGLI PER EVITRARE IL FETCH ALTO NELLE MUDFLAT SEAWARD
        %F(Lside:end,end-10:end)=max(extrafetch,F(Lside:end,end-10:end)); 
        %F(:,end)=extrafetch; 
        else
         %ll=length(F(2:end-1-floor(N*0.5),1));    
        %F(2+floor(N*0.5):end-1,1:1+padding)=2*extrafetch*[1:ll]'/ll*ones(padding+1,1)';
        %F(1,end-padding:end)=extrafetch; 
        %F(1,:)=extrafetch; 
        a=find(h(:,1)<0);if a>0;Lside=a(end)-1;else;Lside=1;end
        F(Lside:end,1:101)=max(extrafetch,F(Lside:end,1:101)); %TOGLI PER EVITRARE IL FETCH ALTO NELLE MUDFLAT SEAWARD
        %F(Lside:end,1:11)=max(extrafetch,F(Lside:end,1:11)); 
        end
    
     
elseif (angleo>=45 & angleo<90)%134)
a=find(h(:,end)<0);if a>0;Lside=N-a(end)-1;else;Lside=N-1;end
F(IND)=cumsumresetEXTRAlateral1(A(IND),extrafetch/dx,Lside)/m*dx;
%(IND)=cumsumreset(A(IND))/m*dx;
%a=find(MASK(
    F(end-Lside:end,:)=extrafetch;       %offshore boudnary
       %ll=length(F(2:end-1-N/2,end));    
       %F(2+N/2:end-1,end-padding:end)=extrafetch*[1:ll]'/ll*ones(padding+1,1)';   
       
elseif (angleo>275 & angleo<=315)
a=find(h(:,1)<0);if a>0;Lside=N-a(end)-1;else;Lside=N-1;end
F(IND)=cumsumresetEXTRAlateral1(A(IND),extrafetch/dx,Lside)/m*dx;
%F(IND)=cumsumreset(A(IND))/m*dx;
    F(end-Lside:end,:)=extrafetch;      %offshore boudnary  
     %ll=length(F(2:end-1-N/2,1));    
     %F(2+N/2:end-1,1:1+padding)=extrafetch*[1:ll]'/ll*ones(padding+1,1)';



else
F(IND)=cumsumreset(A(IND))/m*dx;
end

% padding=5;
% if (angleo<44 | angleo>316)
% F(IND)=cumsumresetEXTRA(A(IND),extrafetch/dx)/m*dx;
%      if angleo<44;
%      ll=length(F(2:end-1,end));    
%        F(2:end-1,end-padding:end)=extrafetch*[1:ll]'/ll*ones(padding+1,1)';
%      else
%      ll=length(F(2:end-1,1));    
%      F(2:end-1,1:1+padding)=extrafetch*[1:ll]'/ll*ones(padding+1,1)';
%     end
% elseif (angleo>44 & angleo<90)%134)
% F(IND)=cumsumresetEXTRAlateral1(A(IND),extrafetch/dx)/m*dx;
%     F(end-padding:end,:)=2*extrafetch;       
% elseif (angleo>270 & angleo<316)
% F(IND)=cumsumresetEXTRAlateral1(A(IND),extrafetch/dx)/m*dx;
%     F(end-padding:end,:)=extrafetch;        
% else
% F(IND)=cumsumreset(A(IND))/m*dx;
% end




% 
% figure
% imagesc(F)
% pause
% 

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

