%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: V struct, variable (1=theta.. etc.), Uarray (needed for peakIDx),
% and sampling type (random/peakPro)
% OUTPUT: [hx,mx,FAx,CRx] : Design matrix with all variables at
% Function that works using V struct. V struct is built from
% classifierWrapper. This will find the mean variable of interest within
% each trial and within each trial type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hx,mx,FAx,CRx] = meanVarfinder_v2(array,var)



Ttype = {'hit','miss','FA','CR'};
% BUILDING AVG THETA AT TOUCH
for d = 1:length(Ttype)
    selvar = array.var.(Ttype{d}){var};
    
    %used to toss out retraction touches
    phasevar = array.var.(Ttype{d}){5};
    keepPro = phasevar<0;

    Ttmp=[array.touchNum.(Ttype{d})];
    yesTouches = find(Ttmp>0);
    meanTtheta = NaN(1,length(Ttmp)); %NaN matrix size = numTrials of trial type
    
    if ~isempty(selvar)
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
        for i = 1:length(tIdx)-1
            rawVals = selvar(tIdx(i):tIdx(i+1));
            keepVals = keepPro(tIdx(i):tIdx(i+1));
            meanTtheta(yesTouches(i)) = mean(rawVals(keepVals));
        end
               
    end
    
    thetatmps{d} = meanTtheta;
    
end  

hx = thetatmps{1}';
mx = thetatmps{2}' ;
FAx = thetatmps{3}' ;
CRx = thetatmps{4}';




