function [Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr]=stratigraphy2D(A,Y1,Y2,Y3,Yb,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,plyr,nlyr,dlyr,tlyrU,tlyrD);;

%to temove the negatives
StY1=max(0,-Y1);
StY2=max(0,-Y2);
StY3=max(0,-Y3);
Y1=Y1+StY1;
Y2=Y2+StY2;
Y3=Y3+StY3;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=Y1+Y2+Y3;

exsE=find(Y>tlyrU & A==1); %accretion past the limit
exsD=find(Y<tlyrD & A==1); %erosion past the limit
    


%accretion past limit -> need to remove material from top cell
    flyrU1=Y1./(Y1+Y2+Y3); 
    flyrU2=Y2./(Y1+Y2+Y3); 
    flyrU3=Y3./(Y1+Y2+Y3); 
    deltaY1=dlyr*flyrU1(exsE);%how much is removed
    deltaY2=dlyr*(flyrU2(exsE));%how much is removed
    deltaY3=dlyr*(flyrU3(exsE));%how much is removed
    Y1(exsE)=Y1(exsE)-deltaY1;%remove it from the active layer
    Y2(exsE)=Y2(exsE)-deltaY2;%remove it from the active layer
    Y3(exsE)=Y3(exsE)-deltaY3;%remove it from the active layer
   
    
    %add the stuff back to the top layer         
        ao=find(plyr(exsE)<nlyr);
        a=find(plyr(exsE)==nlyr);
        
        %add a new layer on top
        plyr(exsE(ao))=plyr(exsE(ao))+1;
        
        %shift the layer compostion, so that you can fill in the top layer
        Ybnew=Yb(exsE(a))+dlyr;%the new thickness of the bottom layer
        firstlyr1=squeeze(flyr1(:,:,1));
        firstlyr2=squeeze(flyr2(:,:,1));
        firstlyr3=squeeze(flyr3(:,:,1));
        %flyrb(exsE(a))=(flyrb(exsE(a)).*Yb(exsE(a))+firstlyr(exsE(a))*dlyr)./Ybnew;%store the fraction in the bottom layer
        flyrb1(exsE(a))=(flyrb1(exsE(a)).*Yb(exsE(a))+firstlyr1(exsE(a))*dlyr)./Ybnew;%store the fraction in the bottom layer
        flyrb2(exsE(a))=(flyrb2(exsE(a)).*Yb(exsE(a))+firstlyr2(exsE(a))*dlyr)./Ybnew;%store the fraction in the bottom layer
        flyrb3(exsE(a))=(flyrb3(exsE(a)).*Yb(exsE(a))+firstlyr3(exsE(a))*dlyr)./Ybnew;%store the fraction in the bottom layer
        Yb(exsE(a))=Ybnew;%increase the thickness of the bed layer
        
        for i=1:nlyr-1;
        temp=squeeze(flyr1(:,:,i));temp2=squeeze(flyr1(:,:,i+1));
        temp(exsE(a))=temp2(exsE(a));
        flyr1(:,:,i)=temp;%make the actual shift in layer composition
        
        temp=squeeze(flyr2(:,:,i));temp2=squeeze(flyr2(:,:,i+1));
        temp(exsE(a))=temp2(exsE(a));
        flyr2(:,:,i)=temp;%make the actual shift in layer composition
        
        temp=squeeze(flyr3(:,:,i));temp2=squeeze(flyr3(:,:,i+1));
        temp(exsE(a))=temp2(exsE(a));
        flyr3(:,:,i)=temp;%make the actual shift in layer composition
        
        end
        
        %eiether way, update the top layer  
        for i=1:nlyr
            temp=squeeze(flyr1(:,:,i));
            a=find(plyr(exsE)==i);temp(exsE(a))=flyrU1(exsE(a));
            flyr1(:,:,i)=temp;
            
            temp=squeeze(flyr2(:,:,i));
            a=find(plyr(exsE)==i);temp(exsE(a))=flyrU2(exsE(a));
            flyr2(:,:,i)=temp;
            
            temp=squeeze(flyr3(:,:,i));
            a=find(plyr(exsE)==i);temp(exsE(a))=flyrU3(exsE(a));
            flyr3(:,:,i)=temp;
        end
        
        
        
        %YOU CAN ONLY TAKE SEDIEMNT FROM ONE LAYER AT THE TIME
%erosion past the limit: need to take from bottom cell

        a=find(plyr(exsD)>=1); %the level occupied is above the base layer
        %plyr(exsD)
        %pause
        ao=find(plyr(exsD)==0);
        
        deltaY1=0*Y;deltaY2=0*Y;deltaY3=0*Y;
        
        %take it from the uppermost layer
        for i=1:nlyr;            %go over all the slices
            aa=find(plyr(exsD(a))==i);
            
            temp=squeeze(flyr1(:,:,i));
            deltaY1(exsD(a(aa)))=temp(exsD(a(aa)))*dlyr;      
            temp(exsD(a(aa)))=NaN;%clear the cell that you eliminated
            flyr1(:,:,i)=temp;
            
            temp=squeeze(flyr2(:,:,i));
            deltaY2(exsD(a(aa)))=(temp(exsD(a(aa))))*dlyr;   
            temp(exsD(a(aa)))=NaN;%clear the cell that you eliminated
            flyr2(:,:,i)=temp;
            
            temp=squeeze(flyr3(:,:,i));
            deltaY3(exsD(a(aa)))=(temp(exsD(a(aa))))*dlyr;   
            temp(exsD(a(aa)))=NaN;%clear the cell that you eliminated
            flyr3(:,:,i)=temp;
            
%             temp(exsD(a(aa)))=NaN;%clear the cell that you eliminated
%             flyr(:,:,i)=temp;
        end
        plyr(exsD(a))=plyr(exsD(a))-1;%lower the elevation of the top layer
        
       
        %take it from the bottom layer (if your at the lowest level)
        deltaYb=min(dlyr,Yb(exsD(ao)));
        deltaY1(exsD(ao))=deltaYb.*flyrb1(exsD(ao));
        deltaY2(exsD(ao))=deltaYb.*(flyrb2(exsD(ao)));
        deltaY3(exsD(ao))=deltaYb.*(flyrb3(exsD(ao)));
        Yb(exsD(ao))=Yb(exsD(ao))-deltaYb;
        
        
        %either way (normal layer or bottom layer), add it!, becuase you already removed it
        Y1=Y1+deltaY1;%how much is added
        Y2=Y2+deltaY2;%how much is added
        Y3=Y3+deltaY3;%how much is added
        
       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%re-add the negative layers
Y1=Y1-StY1;Y2=Y2-StY2;Y3=Y3-StY3;