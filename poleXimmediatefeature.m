touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
[rc] = numSubplots(length(selectedCells));
%%
for k = 12
    curr = U{selectedCells(k)};
    motors = curr.meta.motorPosition';
    
    meanTRraw = nan(curr.k,numel(-25:50));
    velTR = nan(curr.k,1);
    ftvelTR = velTR;
    thetaTR = nan(curr.k,1);
    ampTR = thetaTR;
    phaseTR = thetaTR;
    midpointTR = thetaTR;
    for b = 1:curr.k
        ttouchIdx = [find(curr.S_ctk(9,:,b)==1) find(curr.S_ctk(12,:,b)==1)];
        ftouchIdx = [find(curr.S_ctk(9,:,b)==1)];
        tspks = [squeeze(curr.R_ntk(:,:,b)) nan(1,50)];
        tthetas = squeeze(curr.S_ctk(1,:,b));
        tphase =squeeze(curr.S_ctk(5,:,b));
        tmidpoint=squeeze(curr.S_ctk(4,:,b));
        tamp = squeeze(curr.S_ctk(3,:,b));
        
        if ~isempty(ttouchIdx)
            meanTRraw(b,:)=nanmean(tspks(ttouchIdx'+(-25:50)));
            velTR(b) = nanmean(nanmean(diff(tthetas(ttouchIdx'+(-5:0)),1,2)));
            ftvelTR(b) = nanmean(nanmean(diff(tthetas(ftouchIdx'+(-5:0)),1,2)));
            thetaTR(b) = nanmean(tthetas(ttouchIdx));
            ampTR(b) = nanmean(tamp(ttouchIdx));
            midpointTR(b) = nanmean(tmidpoint(ttouchIdx));
            phaseTR(b) = nanmean(tphase(ttouchIdx));    

        end
       
    end

    keepIdx = ~isnan(meanTRraw(:,1));
    phaseTR = phaseTR(keepIdx,:);
    ampTR = ampTR(keepIdx,:);
    midpointTR = midpointTR(keepIdx,:);
    thetaTR = thetaTR(keepIdx,:);
    velTR = velTR(keepIdx,:);
    ftvelTR = ftvelTR(keepIdx,:);
%     meanTR = smoothdata(meanTRraw(keepIdx,:)')';
    meanTR = meanTRraw(keepIdx,:);
    
    [~,velIdx] = sort(velTR); 
    [~,ftvelIdx] = sort(ftvelTR);
    [~,thetaIdx] = sort(thetaTR); 
    [~,motorIdx] = sort(motors(keepIdx));
    [~,ampIdx] = sort(ampTR); 
    [~,mpIdx] = sort(midpointTR);    
    [~,pIdx] = sort(phaseTR); 
    
    smdata = meanTR(motorIdx,:);
    stdata = meanTR(flipud(thetaIdx),:);
    svdata = meanTR(velIdx,:); %MEAN VELOCITY/TOUCH 
    sftvdata = meanTR(ftvelIdx,:); %FIRST TOUCH VELOCITY 
    
    pdata = meanTR(flipud(pIdx),:);
    adata = meanTR(flipud(ampIdx),:);
    mpdata = meanTR(flipud(mpIdx),:);
    
    
    %FT VEL ANALYSIS DIRECTIONALITY 
    retT = ftvelTR<0;
    proT = ftvelTR>0;
    
    meanbase = mean(mean(sftvdata(:,1:25),2));
    stdbase = std(mean(sftvdata(:,1:25),2));
    postresponses = sftvdata(:,26:75);
    xresp = mean(postresponses,2);
    zscore = (xresp - meanbase) ./ stdbase;
    
    comp{k} = nan(1000,2);
    comp{k}(1:sum(retT),1) = zscore(retT);
    comp{k}(1:sum(proT),2) = zscore(proT);
    [p(k),~,stats]=anova1(comp{k},[],'off');
%     vals =multcompare(stats);
    
    %FT VEL ANALYSIS STRENGTH OF TOUCH
    absftv = abs(ftvelTR);
    [~,aftvelIdx] = sort(absftv);
    
    figure(100);subplot(4,8,k);
    hold on;scatter(absftv,zscore,3,'filled')
    cor(k) = corr(absftv,zscore);
    
    %THETA TUNING
    meanbase = mean(mean(stdata(:,1:25),2));
    stdbase = std(mean(stdata(:,1:25),2));
    postresponses = stdata(:,26:75);
    xresp = mean(postresponses,2);
    zscore = (xresp - meanbase) ./ stdbase;
    
    
    
end
gmeans = cell2mat(cellfun(@nanmean,comp,'uniformoutput',0)');

figure(44);clf
scatter(gmeans(:,1),gmeans(:,2),[],'k','filled')
axis square
set(gca,'xlim',[0 5],'ylim',[0 5],'xtick',0:2.5:5,'ytick',0:2.5:5)
hold on; plot([0 5],[0 5],'-.k')
xlabel('retraction touch z score');ylabel('protraction touch z score')

tmp = gmeans(12,:);
hold on; scatter(tmp(:,1),tmp(:,2),[],'r','filled')

%ZSCORE PLOT
figure(8);clf;bar(1:length(zscore),zscore,'k')
set(gca,'xtick',[],'xlim',[0 length(zscore)+1])

    %%
        
    %PLOTTING
    figure(48);clf
    subplot(2,3,1)
    imagesc(smdata)
    set(gca,'xtick',1:25:76,'xticklabel',[-25:25:50],'ytick',[])
    hold on; plot([26 26],[0 size(meanTR,1)],'-.w')
    title('motor sorted') 
    
    subplot(2,3,2)
    imagesc(stdata)
    set(gca,'xtick',1:25:76,'xticklabel',[-25:25:50],'ytick',[])
    hold on; plot([26 26],[0 size(meanTR,1)],'-.w')
    title('angle sorted') 
    
    subplot(2,3,3)
    imagesc(svdata)
    set(gca,'xtick',1:25:76,'xticklabel',[-25:25:50],'ytick',[])
    hold on; plot([26 26],[0 size(meanTR,1)],'-.w')
    title('preTvel sorted') 
    
    subplot(2,3,4)
    imagesc(adata)
    set(gca,'xtick',1:25:76,'xticklabel',[-25:25:50],'ytick',[])
    hold on; plot([26 26],[0 size(meanTR,1)],'-.w')
    title('amp sorted')
    
     subplot(2,3,5)
    imagesc(mpdata)
    set(gca,'xtick',1:25:76,'xticklabel',[-25:25:50],'ytick',[])
    hold on; plot([26 26],[0 size(meanTR,1)],'-.w')
    title('midpoint sorted') 
    
     subplot(2,3,6)
    imagesc(pdata)
    set(gca,'xtick',1:25:76,'xticklabel',[-25:25:50],'ytick',[])
    hold on; plot([26 26],[0 size(meanTR,1)],'-.w')
    title('phase sorted')
    
    figure(889);clf
    imagesc(imgaussfilt(stdata,[1 1],'padding','replicate'))
    set(gca,'xtick',1:25:76,'xticklabel',[-25:25:50],'ytick',[])
    colorbar
    caxis([0 .2])
    hold on; plot([26 26],[0 size(meanTR,1)],'-.w')

    