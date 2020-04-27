function [E]=wavediffusionstep(E,i,c,cg,d,ho,sigma,dx,angle,Cbr,Cbed,alpha,angle0,N,M,A,AW,Awup,p,ii,jj,periodic,Holateral)


%if there is at least an active cell..., then calculate the diffusion
%if sum(Active==1)>0
if sum(ho>0)>0
    
    
%Use the fomrula of Mase
%use chain rule
%get this:  + 1/2*CCg*Eyy (CCg)y*Ey  
%solved the first one as diffusion
%solve second one as edvection

%when you stretch the coordinates, you also stretch the dx!!!
stretchdx=1/cos(angle/180*pi);

        DD=(cg.*c)';
        DD(d<0.1 | A(i,:)==0)=0;
        DDp=([DD(2:end); DD(end)]-[DD(1); DD(1:end-1)])/2;DDp(1)=0;DDp(end)=0;
        s=[];S=0*p;
        %%%%%%%%%%%%%%%%%%
        for k=[-1 1];
            
%         if periodic==0;     q=p+k;if k==1;a=[1:M-1];elseif k==-1;a=[2:M];end 
%         elseif periodic==1; q=1+mod(p+k-1,M);a=[1:M];
%         end
        q=1+mod(p+k-1,M);a=[1:M];%%KINDA force the periodic options in chosiing the cell (Dleaware wave diffusion)
        
        %Lateral diffraction (diffusion)     
        sigma=min(2*pi/1,max(2*pi/20,sigma));
        cg=min(100,max(1,cg));
        if length(sigma)>1
        value=alpha/dx^2 /2./sigma(p(a))'  *0.5.*(DD(p(a))+DD(q(a)))/2  *dx*stretchdx./cg(p(a))'*(cos(angle/180*pi))^2;   %*stretchdx
        else
        value=alpha/dx^2 /2./sigma  *0.5.*(DD(p(a))+DD(q(a)))/2  *dx*stretchdx./cg(p(a))'*(cos(angle/180*pi))^2;   %*stretchdx
        end    
            
        up=[];F=0;
        if (k==1);UR=DDp;up=find(UR>0);F=UR(up);end
        if (k==-1);UR=DDp;up=find(UR<0);F=-UR(up);end
        if length(sigma)>1
        value(up)=value(up)+alpha/dx^2 /2./sigma(up)'  .*F  *dx*stretchdx./cg(up)'*(cos(angle/180*pi))^2;
        else
        value(up)=value(up)+alpha/dx^2 /2./sigma  .*F  *dx*stretchdx./cg(up)'*(cos(angle/180*pi))^2;
        end
        %note that you gor a diveded by dx^2. One dx is ffrom DDp (dx not included there)
        %the other dx is from the derivate of E
        
        %NOTE:
        %in both case you have *dx/cg, which is the time elapsed
        %you need to strech the dx, becuase with the same cg it needs to
        %travle farther if it is at anlgle
        %then you get the cos^2 form the Mase formulae
        
            %REMOVE THE LATERAL DIFFUSION
            if periodic==0;%DO IT ONLY IF IT IS NOT PERIODIC!!!!!!
            value=value.*(AW(p(a))==0)';%zero diffusion at the lateral boundary if not periodic
            end
        
        S(p(a))=S(p(a))+value; %exit from that cell
        s=[s;-value]; %gain from the neigborh cell
        end

        

        %Shoaling, bed friction, and breaking
        s=[s;S+1];
        %s=[s;S+(  (cg(i,p)+dx*stretchdx*(bedfr(i,p))+brk+whitecap)./cg(i+1,p))'];

        A = sparse(ii(abs(s)>0),jj(abs(s)>0),s(abs(s)>0));%solve the matrix inversion
        %A = sparse(ii,jj,s);%solve the matrix inversion
        
        E=(A\E')';
        
        
        
%         %Set the lateral bouundaries if not periodics. Two options:
%         %no-gradient (+-1) or zero wave heigth (+-2)
%         if periodic==0 %DO IT ONLY IF IT IS NOT PERIODIC!!!!!!
%             if angle>0%angle<0
%                 a=find(AW(:)==-1); %find the incomign boundary
%                 if a>0;
%                 b=find(AW(1+mod(-1+a-1,M))==0);
%                 E(a)=E(a(b)); 
%                 end   
%                 a=find(AW(:)==-2); %find the incomign boundary
%                 if a>0;
%                 b=find(AW(1+mod(-1+a-1,M))==0); %the first free cel;
%                 E(a)=E(a(b))*0+Holateral^2/16; 
%                 end 
%             elseif angle<0%angle>0
%                 a=find(AW(:)==1); %find the incomign boundary
%                 if a>0;
%                 b=find(AW(1+mod(-1+a+1,M))==0);
%                 E(a)=E(a(b)); 
%                 end   
%                 a=find(AW(:)==2); %find the incomign boundary
%                 if a>0;
%                 b=find(AW(1+mod(-1+a+1,M))==0);
%                 E(a)=E(a(b))*0+Holateral^2/16;  
%                 end  
%             end
%         end
%         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                
        %value=alpha/dx^2 /(2*(2*pi/T)) *((DD(p(a))+DD(q(a)))/2-1/2*DD(p(a))) *dx*stretchdx./cg(i+1,p(a))';   
        %value=alpha/dx^2 /(2*(2*pi/T)) *0.5*(DD(p(a))+DD(q(a)))/2  *dx.*STRDX'./cg(i+1,p(a))';   %*stretchdx
        %value=alpha/dx^2 /(2*(2*pi/T)).*min(DD(p(a)),DD(q(a))) *dx*stretchdx./cg(i+1,p(a))'; %GOOD     
        %value=alpha/dx^2 /(2*(2*pi/T)) *min(DD(p(a)),DD(q(a))) *dx*stretchdx;   
        %value=alpha/dx^2 /(2*(2*pi/T)) *min(DD(p(a)),DD(q(a))) *dx*stretchdx./cg(i+1,p(a))';
        %value=alpha/dx^2 ./max(kwave(i,p(a)),kwave(i,q(a)))' *dx*stretchdx;   
        %value=alpha/dx^2 ./kwave(i+1,p(a))' *dx*stretchdx;  
           
        
        
        
        
end
 