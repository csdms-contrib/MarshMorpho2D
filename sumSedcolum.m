function sumY=sumSedcolum(Yb,flyrb,flyr,dlyr,Yi)
%sumY=[Yb.*flyrb+nansum(flyr.*isfinite(flyr),3)*dlyr+Yi];

%PORCODIOO QUESTO ERA QUELLO ORIGINAL
sumlyr=nansum(flyr.*isfinite(flyr),3)*dlyr; %THE ORIGINAL

% sumlyr=sum(flyr,3)*dlyr;


sumlyr(isnan(sumlyr))=0;
sumY=[Yb.*flyrb+sumlyr+Yi];