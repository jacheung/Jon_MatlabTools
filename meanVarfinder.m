%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: V struct, variable (1=theta.. etc.), Uarray (needed for peakIDx),
% and sampling type (random/peakPro)
% OUTPUT: [hx,mx,FAx,CRx] : Design matrix with all variables at
% Function that works using V struct. V struct is built from
% classifierWrapper. This will find the mean variable of interest within
% each trial and within each trial type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hx,mx,FAx,CRx] = meanVarfinder (array,var,Uarray,samp)



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
       
                   
    else
% Filling NaN values with distribution of FA/CR
%         if d == 3 || d == 4 %if FA or CR has no touches in it at all, resample those NaN values with FA and CR theta touches 
%         ThetaGroup=[array.var.FA{var} array.var.CR{var}];
%         meanTtheta = datasample(ThetaGroup,length(meanTtheta));
%         end
%        
    end
    
    thetatmps{d} = meanTtheta;
    

end  

if strcmp(samp,'random')
    %Filling all NaN values with randomly sampled theta vals
    thetas = [thetatmps{1} thetatmps{2} thetatmps{3} thetatmps{4}];
    disp(['Filling in ' num2str(sum(isnan(thetas))) ' values out of ' num2str(numel(thetas))])
    uniformthetas = min(thetas):.5:max(thetas);   
    thetas(isnan(thetas))=[];
    for d = 1:4
        nansamps = find(isnan(thetatmps{d}));
        for i = 1:sum(isnan(thetatmps{d}))
            thetatmps{d}(nansamps(i)) = datasample(thetas,1); %sampling from distribution of all touches and their data (skewed towrads go)  
%             thetatmps{d}(nansamps(i)) = datasample(uniformthetas,1); %sampling from uniform distribution of obtained theta values 
        end
    end
    
elseif strcmp(samp,'resampSameD')
    for d = 1:2 %resample FA/CR from nogo distribution
        thetas = [thetatmps{1} thetatmps{2}];
        thetas(isnan(thetas))=[];
        
        nansamps = find(isnan(thetatmps{d}));
        for i = 1:sum(isnan(thetatmps{d}))
            thetatmps{d}(nansamps(i)) = datasample(thetas,1);
        end
    end
    
    for d = 3:4 %resample FA/CR from nogo distribution
        thetas = [thetatmps{3} thetatmps{4}];
        thetas(isnan(thetas))=[];
        
        nansamps = find(isnan(thetatmps{d}));
        for i = 1:sum(isnan(thetatmps{d}))
            thetatmps{d}(nansamps(i)) = datasample(thetas,1);
        end
    end
        
elseif strcmp(samp,'peakPro')
    
    peakIDX = Uarray.peakIdx;
    thetas = squeeze(Uarray.S_ctk(1,:,:));
    
    trialNums =floor((peakIDX/4000)+1);
    tNums = unique(trialNums) ;
    for i = 1:length(tNums) 
        trialIdx = find(trialNums==tNums(i));
        peakPro(tNums(i)) = max(thetas(trialIdx));
    end
    
    hits = peakPro(find(array.trialNums.matrix(1,:)==1));
    miss = peakPro(find(array.trialNums.matrix(2,:)==1));
    FA = peakPro(find(array.trialNums.matrix(3,:)==1));
    CR = peakPro(find(array.trialNums.matrix(4,:)==1));
    
    peakProMat = {hits,miss,FA,CR};
    
    for d = 1:4
        nansamps = find(isnan(thetatmps{d}));
        for i = 1:sum(isnan(thetatmps{d}))
            thetatmps{d}(nansamps(i)) = peakProMat{d}(nansamps(i));
        end
    end
    
else 
    error('Choose what to fill NaN values with: random , resampSameD , or peakPro')
    
end

    
    
hx = thetatmps{1};
mx = thetatmps{2} ;
FAx = thetatmps{3} ;
CRx = thetatmps{4};




