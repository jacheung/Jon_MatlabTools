touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
[rc] = numSubplots(length(selectedCells));

%%
close all;
window = (1:50);
v=1;
rtcomp = zeros(length(selectedCells),2);
for k = 12
    curr = U{selectedCells(k)};
    motors = curr.meta.motorPosition';
    
    
     touchIdx = [find(curr.S_ctk(9,:,:)==1) ;find(curr.S_ctk(12,:,:)==1)];
     touchIdxoff = [find(curr.S_ctk(10,:,:)==1); find(curr.S_ctk(13,:,:)==1)];
    spks = squeeze(curr.R_ntk(:,:,:));
    thetas = squeeze(curr.S_ctk(1,:,:));
    vel = squeeze(curr.S_ctk(2,:,:));
    
    %touch Duration
    durSpk = zeros(length(touchIdx),1);
    for g = 1:length(touchIdx)
        durSpk(g) = nansum(spks(touchIdx(g):touchIdxoff(g)));
    end
    
    %touch Count
    ct = zeros(curr.k,1);
    ctspks = zeros(curr.k,1);

    for b = 1:curr.k
        ttouchIdx = [find(curr.S_ctk(9,:,b)==1) find(curr.S_ctk(12,:,b)==1)];
        tspks = [squeeze(curr.R_ntk(:,:,b)) nan(1,50)];

        if ~isempty(ttouchIdx)

            ctspks(b) = nansum(nansum(tspks(ttouchIdx'+window),2));
            ct(b) = numel(ttouchIdx); 
        end
       
    end

    touchDur = touchIdxoff-touchIdx;
    thetaAtT = thetas(touchIdx);
    preTvel = nanmean(diff(thetas(touchIdx+(-5:0)),1,2),2);
    
    
    featSpkRaw = spks(touchIdx+window);
    featSpk = nanmean(spks(touchIdx+window),2)*1000;
    
    stmodel = fitlm(featSpkRaw,thetaAtT);
    stcorr=sqrt(stmodel.Rsquared.Ordinary);
    

    
    
    %correlation matrix 
    % NOTE HERE that correlation for spks for 1:50ms post touch are for
    % touch angle and preTvel. Duration of touch is correlated with number
    % of spks during whole touch window. 
    mat = [touchDur thetaAtT preTvel featSpk];
    durcor = corr([touchDur durSpk]);
    ctcorr = corr([ct ctspks]);
    tmpCor = corr(mat,'rows','complete');
    tmpCor(tmpCor == tmpCor(4))=durcor(2);
    popCor{k} = tmpCor;
    
    rtcomp(k,:) = [popCor{k}(4,2) stcorr];
    
    for d = 1:3
        figure(80+d);subplot(rc(1),rc(2),k);
        if abs(popCor{k}(size(popCor{k},1),d))>.5
            if d==1
                scatter(mat(:,d),durSpk,'b.')
            else
                scatter(mat(:,d),featSpk,'b.')
            end
            
        else
            if d==1
                scatter(mat(:,d),durSpk,'k.')
            else
                scatter(mat(:,d),featSpk,'k.')
            end
        end
        set(gca,'xlim',[min(mat(:,d)) max(mat(:,d))])
    end

    figure(80);
    subplot(rc(1),rc(2),k);
    if abs(ctcorr(2))>.5
        scatter(ct,ctspks,'b.')
    else
        scatter(ct,ctspks,'k.')
    end
    set(gca,'xlim',[0 25])
end


figure(39);clf
scatter(rtcomp(:,1).^2,rtcomp(:,2).^2,[],'k','filled')
hold on; plot([0 1],[0 1],'-.k')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.5:1,'ytick',0:.5:1)
axis square
xlabel('Rate lm Rsq');ylabel('Timing lm Rsq') 


