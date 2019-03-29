function [numWhisks,whiskTimes] = numWhiskPreTouch(U)
for i = 1:length(U)
    [whisks] = findMaxMinProtraction(U{i},5,'avail2end');
    whisktnums = ceil(whisks.peakidx./U{i}.t);
    peaks= whisks.peakidx;
    ptheta = whisks.theta;
    
    
    touched = squeeze(nansum(U{i}.S_ctk(9,:,:)));
    idx_touch = find(touched==1);
    idx_notouch = find(touched==0);
    
    %finding med lick/median touch time
    clear lix
    clear touch
        pOnset = round(U{i}.meta.poleOnset(1)*1000);
        samplingPeriod = 750;
        
        for k = 1:U{i}.k
            lix(k) =min([find(U{i}.S_ctk(16,pOnset+samplingPeriod:4000,k)==1,1) 4000])+pOnset+samplingPeriod;
            touch(k) = min([find(U{i}.S_ctk(9,pOnset:4000,k)==1,1) 4000])+pOnset;
        end
        medlt = median(lix(lix<4000));
        medt = median(touch(touch<4000));
        
    
    %FIND ALL NUM OF WHISK POST FIRST TOUCH!
    ttpre = nan(length(idx_touch),20);
    ttpost = nan(length(idx_touch),20);
    clear preTwt postTwt
    
    for b = 1:length(idx_touch)
        tr = idx_touch(b);
        firstt = find(U{i}.S_ctk(9,:,tr)==1);
        
        tmp2=peaks(find(whisktnums==tr));
        thetapeaks = ptheta(find(whisktnums==tr));
        whisktimes = mod(tmp2,U{i}.t);
        
        pret = whisktimes<medt;
        postt = whisktimes>medt & whisktimes<medlt;
        
        preTwt{b} = whisktimes(pret);
        postTwt{b} = whisktimes(postt);
        
        ttpre(b,1:length(thetapeaks(pret))) = thetapeaks(pret)';
        ttpost(b,1:length(thetapeaks(postt))) = thetapeaks(postt)';
        
    end
    
    whiskTimes.touchTrials.trialNumbers{i} = idx_touch;
    whiskTimes.touchTrials.preFirstTouch{i} =  preTwt;
    whiskTimes.touchTrials.postFirstTouch{i} = postTwt;
    
    
    nttpre = nan(length(idx_notouch),20);
    nttpost = nan(length(idx_notouch),20);
    clear preTwt postTwt
    for b = 1:length(idx_notouch)
        tr = idx_notouch(b);
        
        tmp2=peaks(find(whisktnums==tr));
        thetapeaks = ptheta(find(whisktnums==tr));
        whisktimes = mod(tmp2,U{i}.t);
        
        pret = whisktimes<medt;
        postt = whisktimes>medt & whisktimes<medlt;
        
        preTwt{b} = whisktimes(pret);
        postTwt{b} = whisktimes(postt);
        
        nttpre(b,1:length(thetapeaks(pret))) = thetapeaks(pret)';
        nttpost(b,1:length(thetapeaks(postt))) = thetapeaks(postt)';
    end
    
    whiskTimes.nontouchTrials.trialNumbers{i} = idx_notouch;
    whiskTimes.nontouchTrials.preFirstTouch{i} =  preTwt;
    whiskTimes.nontouchTrials.postFirstTouch{i} = postTwt;
    
    numWhisks.touchTrials.preFirstTouch{i} = sum(~isnan(ttpre),2);
    numWhisks.touchTrials.postFirstTouch{i} = sum(~isnan(ttpost),2);
    numWhisks.nontouchTrials.preFirstTouch{i} = sum(~isnan(nttpre),2);
    numWhisks.nontouchTrials.postFirstTouch{i} = sum(~isnan(nttpost),2);
    
    
    
end