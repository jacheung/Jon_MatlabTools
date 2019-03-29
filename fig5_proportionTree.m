preDecisionMaskFull = preDecisionTouchMat(U);

preD_protraction_Touches = cell(1,10);
for i = 1:length(U)
    array = U{i};
    preDecisionMask = preDecisionMaskFull{i};
    touchmat = nan(size(preDecisionMask));
    touchIdx = [find(array.S_ctk(9,:,:)==1) ;find(array.S_ctk(12,:,:)==1)];
    touchmat(touchIdx)=1;
    
    phaseMask = squeeze(array.S_ctk(5,:,:));
    proPhaseMask = double(phaseMask<=0);
    proPhaseMask(proPhaseMask==0)=nan;
    
    preD_protraction_Touch_mask = touchmat .* preDecisionMask .* proPhaseMask; %protraction touches only
%     preD_protraction_Touch_mask = touchmat .* preDecisionMask;

    preD_protraction_Touches{i} = nansum(preD_protraction_Touch_mask);
    
    
    hit = intersect(find(array.meta.trialType==1), find(array.meta.trialCorrect==1));
    miss = intersect(find(array.meta.trialType==1), find(array.meta.trialCorrect==0));
    FA = intersect(find(array.meta.trialType==0), find(array.meta.trialCorrect==0));
    CR = intersect(find(array.meta.trialType==0), find(array.meta.trialCorrect==1));
    
    lix = zeros(1,length(preD_protraction_Touches{i}));
    nolix = lix;
    lix([hit FA]) = 1;
    nolix([miss CR]) = 1;
    go = array.meta.trialType==1;
    nogo = array.meta.trialType==0;
    touch = find(preD_protraction_Touches{i}>0);
    notouch = find(preD_protraction_Touches{i}==0);
    
    outputPro.propTouchGo(i) = mean(go(touch));
    outputPro.propTouchNogo(i) = 1-mean(go(touch));
    
    outputPro.propNoTouchGo(i) = mean(go(notouch));
    outputPro.propNoTouchNogo(i) = 1-mean(go(notouch));
    
    outputPro.propTouchGoLick(i) = mean(lix(intersect(find(go==1),touch)));
    outputPro.propTouchGoNoLick(i) = 1-mean(lix(intersect(find(go==1),touch)));
    
    outputPro.propTouchNoGoLick(i) = mean(lix(intersect(find(nogo==1),touch)));
    outputPro.propTouchNoGoNoLick(i) = 1-mean(lix(intersect(find(nogo==1),touch)));
    
    outputPro.propNoTouchGoLick(i) = mean(lix(intersect(find(go==1),notouch)));
    outputPro.propNoTouchGoNoLick(i) = 1-mean(lix(intersect(find(go==1),notouch)));
    
    outputPro.propNoTouchNoGoLick(i) = mean(lix(intersect(find(nogo==1),notouch)));
    outputPro.propNoTouchNoGoNoLick(i) = 1-mean(lix(intersect(find(nogo==1),notouch)));
    
end



