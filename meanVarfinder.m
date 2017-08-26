%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: V struct 
% OUTPUT: [hx,mx,FAx,CRx] : Design matrix with all variables at
% Function that works using V struct. V struct is built from
% classifierWrapper. This will find the mean variable of interest within
% each trial and within each trial type. 
function [hx,mx,FAx,CRx] = meanVarfinder (array,var)

meanTtheta =[];
thetas = array.var.CR{var};
Ttmp=[array.touchNum.CR];
yesTouches = find(Ttmp>0);

%finding index of last touches within each trial
tmp2 = [Ttmp zeros(1,length(Ttmp))];
for i = 1:length(Ttmp)
    tmp0 = [zeros(1,i) Ttmp zeros(1,length(Ttmp)-i)];
    tmp2 = [tmp2 ;tmp0];
end
tIdx = sum(tmp2(:,1:length(Ttmp)));
    if tIdx(1) == 0 % added this case for when 1st trial = 0 touches, just replacing it w/ 1 instead 
        tIdx = [unique(tIdx)];
        tIdx(1) = 1;
    else 
         tIdx = [1 unique(tIdx)];
    end

%find mean of all touches within those trials
meanTtheta = NaN(1,length(Ttmp));
for i = 1:length(tIdx)-1
    meanTtheta(yesTouches(i)) = mean(thetas(tIdx(i):tIdx(i+1)));
end
meanTtheta(isnan(meanTtheta))=nanmean(meanTtheta);

CRx = meanTtheta;
NANTHETA = nanmean(meanTtheta);
%%
Ttype = {'hit','miss','FA'};
% BUILDING AVG THETA AT TOUCH
for d = 1:length(Ttype)
    meanTtheta =[];
    thetas = array.var.(Ttype{d}){var};
    Ttmp=[array.touchNum.(Ttype{d})];
    yesTouches = find(Ttmp>0);
    
    if ~isempty(thetas)
    %finding index of last touches within each trial
    tmp2 = [Ttmp zeros(1,length(Ttmp))];
    for i = 1:length(Ttmp)
        tmp0 = [zeros(1,i) Ttmp zeros(1,length(Ttmp)-i)];
        tmp2 = [tmp2 ;tmp0];
    end
    
    tIdx = sum(tmp2(:,1:length(Ttmp)));
    if tIdx(1) == 0 % added this case for when 1st trial = 0 touches, just replacing it w/ 1 instead 
        tIdx = [unique(tIdx)];
        tIdx(1) = 1;
    else 
         tIdx = [1 unique(tIdx)];
    end
   
    
    %find mean of all touches within those trials
    meanTtheta = NaN(1,length(Ttmp));
    for i = 1:length(tIdx)-1
        meanTtheta(yesTouches(i)) = mean(thetas(tIdx(i):tIdx(i+1)));
    end
    meanTtheta(isnan(meanTtheta))=NANTHETA; %found from CR trials.
    else
        meanTtheta = [];
    end
    thetatmps{d} = meanTtheta;
    
end

hx = thetatmps{1};
mx = thetatmps{2} ;
FAx = thetatmps{3} ;

