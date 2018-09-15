%Function for finding the touch times from whisk onset to touch. Couple
%specifications is...
% 1) only looks within window of pole availability to lick
% 2) only looks at protraction touches
% 3) toss touch times >100ms
%Written by JC 180315 along with findMaxMinProtraction



function [ tperWhisksorted, trialNums,stats] = uber_whisksThatTouch(array)

%P struct for finding peak and trough protraction values
% [P] = findMaxMinProtraction(array,'avail2lick');
[P] = findMaxMinProtraction(array,5,'avail2lick');
[ALLwhisks] = findMaxMinProtraction(array,2.5);

%masks for touches 
[objmask]= assist_touchmasks(array);
mask = objmask.availtolick;

%Variables for classification
nogos =find(array.meta.trialType ==0);
gos = find(array.meta.trialType ==1);
corr = find(array.meta.trialCorrect ==1);
incorr = find(array.meta.trialCorrect == 0);
hits = intersect(gos,corr);
FA = intersect(nogos,incorr);
CR = intersect(nogos,corr);
miss = intersect(gos,incorr);


% FIND SPAN OF WHISKING Inherited from searchbiasV3.0

        %Protraction and Retraction 
        matchedIdxRet = nan(length(P.peakidx),1);
        matchedIdxPro = nan(length(P.peakidx),1);
        for k = 1:length(P.peakidx)
            matchedtmpRet  = find(ALLwhisks.troughIdx>P.peakidx(k),1);
            matchedtmpPro =  find(ALLwhisks.troughIdx<P.peakidx(k),1,'last');
            if ~isempty(matchedtmpRet)
                matchedIdxRet(k) = ALLwhisks.troughIdx(matchedtmpRet);
            end
            if ~isempty(matchedtmpPro)
                matchedIdxPro(k) = ALLwhisks.troughIdx(matchedtmpPro);
            end
        end
        
        whiskDurationPro = P.peakidx - matchedIdxPro ;
        tossIdxPro = find(whiskDurationPro>60);
        spanIdxPro = [matchedIdxPro P.peakidx];
        spanIdxPro(tossIdxPro,:) = [];
        spanIdxPro(sum(isnan(spanIdxPro), 2) >= 1, :) = []; 

        
        whiskDurationRet = matchedIdxRet-P.peakidx;
        tossIdxRet = find(whiskDurationRet>60);
        spanIdxRet = [P.peakidx matchedIdxRet];
        spanIdxRet(tossIdxRet,:) = [];
        spanIdxRet(sum(isnan(spanIdxRet), 2) >= 1, :) = [];        


% FIND TOUCHES WITHIN MASKPERIOD PERIOD
 touchIdx = [find(array.S_ctk(9,:,:)==1);find(array.S_ctk(12,:,:)==1)];
 touchIdx= intersect(touchIdx,find(mask==1));
 thetas = squeeze(array.S_ctk(1,:,:));
 
 whiskmatPro = nan(length(spanIdxPro),61);
 for q=1:length(spanIdxPro)
     tmpidx = spanIdxPro(q,1):spanIdxPro(q,2);
     whiskmatPro(q,1:length(tmpidx))=tmpidx;
 end
 
 whiskmatRet= nan(length(spanIdxRet),61);
 for q=1:length(spanIdxRet)
     tmpidx = spanIdxRet(q,1):spanIdxRet(q,2);
     whiskmatRet(q,1:length(tmpidx))=tmpidx;
 end
 

 % CALCULATIONS OF BASIC STATS 
 proTouches = numel(intersect(whiskmatPro,touchIdx));
 retTouches = numel(intersect(whiskmatRet,touchIdx));
 
 %total proportion of touches accounted for 
 proportionaccountedtouches =  (proTouches+retTouches)/numel(touchIdx);
 
 accountedtouchesPertrial = zeros(array.k,1);
 allaccountedTouches = [intersect(whiskmatPro,touchIdx) ; intersect(whiskmatRet,touchIdx)];
 touchintrialIdx = ceil(allaccountedTouches/array.t);
 for g = 1:array.k
     numTouches = sum(touchintrialIdx==g);
     accountedtouchesPertrial(g) = numTouches;
 end
 
 numWhiskPertrial = zeros(array.k,1);
 accountedWhisks = unique([spanIdxRet(:,1);spanIdxPro(:,2)]);
 whisksintrialIdx = ceil(accountedWhisks/array.t);
 for p=1:array.k
     numWhisks = sum(whisksintrialIdx ==p);
      numWhiskPertrial(p) = numWhisks;
 end
 
 touchesperwhisk = accountedtouchesPertrial./numWhiskPertrial;
 
 %Main outputs 
 stats.proportionTouches =  proportionaccountedtouches;
 stats.touchstats.totalNumTouches = numel(touchIdx);
 stats.touchstats.protractionTouches = proTouches;
 stats.touchstats.retractionTouches = retTouches; 
 
 tperWhisksorted{1} = touchesperwhisk(hits);
 tperWhisksorted{2} = touchesperwhisk(FA);
 tperWhisksorted{3} = touchesperwhisk(CR);
 tperWhisksorted{4} = touchesperwhisk(miss);
 
 trialNums{1} = hits;
 trialNums{2} = FA;
 trialNums{3} = CR;
 trialNums{4} = miss; 
    
%  tmp=((touchIdx./4000)-floor(touchIdx./4000))*4000;
% figure;histogram(tmp,0:250:4000)
 
%% test a trial to make sure that theta and touches are aligned right         
%             trial = randi([1 array.k],1,1); %shifted by 1 (ex. trial=0 ..> trial =1)
%             
%             figure(354793);clf;plot(thetas(:,trial+1),'k');
%             xlabel('Time from trial start (ms)');ylabel('Whisker position')
%    
%             validxpro=spanIdxPro(floor(spanIdxPro(:,2)/array.t)==trial,:);
%             validxret=spanIdxRet(floor(spanIdxRet(:,2)/array.t)==trial,:);
%             validtouchx = touchIdx(floor(touchIdx/array.t)==trial);
%             
%             
%             for k = 1:size(validxpro,1)
%                 xvals = (validxpro(k,1):validxpro(k,2))-(trial*array.t);
%                 hold on; plot(xvals,thetas(validxpro(k,1):validxpro(k,2)),'b','linewidth',2)
%             end
%             
%             for k = 1:size(validxret,1)
%                 xvals = (validxret(k,1):validxret(k,2))-(trial*array.t);
%                 hold on; plot(xvals,thetas(validxret(k,1):validxret(k,2)),'b','linewidth',2)
%             end
%             
%               touchexcurx=round(((validtouchx/array.t)-trial)*array.t);
%             
%              for i = 1:length(touchexcurx)
%              hold on; scatter(touchexcurx(i),thetas(touchexcurx(i),trial+1),'go','filled')
%             end
            
