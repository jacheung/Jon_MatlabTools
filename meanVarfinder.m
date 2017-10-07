%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: V struct
% OUTPUT: [hx,mx,FAx,CRx] : Design matrix with all variables at
% Function that works using V struct. V struct is built from
% classifierWrapper. This will find the mean variable of interest within
% each trial and within each trial type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hx,mx,FAx,CRx] = meanVarfinder (array,var,varargin)

nanvals = varargin{1};

Ttype = {'hit','miss','FA','CR'};
% BUILDING AVG THETA AT TOUCH
for d = 1:length(Ttype)
    meanTtheta =[];
    thetas = array.var.(Ttype{d}){var};
    Ttmp=[array.touchNum.(Ttype{d})];
    yesTouches = find(Ttmp>0);
    meanTtheta = NaN(1,length(Ttmp)); %NaN matrix size = numTrials of trial type
    
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
        for i = 1:length(tIdx)-1
            meanTtheta(yesTouches(i)) = mean(thetas(tIdx(i):tIdx(i+1)));
        end
        

                   
%     else
% Filling NaN values with distribution of FA/CR
%         if d == 3 || 4 %if FA or CR has no touches in it at all, resample those NaN values with FA and CR theta touches 
%         ThetaGroup=[array.var.FA{var} array.var.CR{var}];
%         meanTtheta = datasample(ThetaGroup,length(meanTtheta));
%         end
       
    end
    
    thetatmps{d} = meanTtheta;
    

end  

thetas = [thetatmps{1} thetatmps{2} thetatmps{3} thetatmps{4}];
thetas(isnan(thetas))=[];
for d = 1:4
    nansamps = find(isnan(thetatmps{d}));
    for i = 1:sum(isnan(thetatmps{d}))
        thetatmps{d}(nansamps(i)) = datasample(thetas,1);
    end
end

hx = thetatmps{1};
mx = thetatmps{2} ;
FAx = thetatmps{3} ;
CRx = thetatmps{4};




