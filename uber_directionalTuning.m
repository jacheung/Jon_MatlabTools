touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
[rc] = numSubplots(length(selectedCells));
%%
touchForcevals=cell(1,length(selectedCells));
for k = 1:length(selectedCells)
    curr = U{selectedCells(k)};
    motors = curr.meta.motorPosition';
    
    meanTRraw = nan(curr.k,numel(-25:50));
    velTR = nan(curr.k,1);
    ftvelTR = velTR;
    thetaTR = nan(curr.k,1);
    FTTRraw = meanTRraw;
    
    % POPULATION comparison pro vs ret
    ttouchIdx = [find(curr.S_ctk(9,:,:)==1) ;find(curr.S_ctk(12,:,:)==1)];
    tphase =squeeze(curr.S_ctk(5,:,:));
    tspks = squeeze(curr.R_ntk(1,:,:));
    
    ret = find(tphase(ttouchIdx)>0);
    pro = find(tphase(ttouchIdx)<0);
    % baseline subtracted spiking response
    sspks = sum(tspks(ttouchIdx + (0:50)),2) - sum(tspks(ttouchIdx + (-25:0)),2);
    
    
    [~,prp(k)] = ttest2(sspks(ret),sspks(pro));

    prcmeansraw(k,:) = [mean(sspks(ret)) mean(sspks(pro))] ;
    prccounts(k,:) = [numel(ret) numel(pro)];
    
   figure(102);subplot(4,8,k);
   if prp(k)<.05
   scatter(tphase(ttouchIdx),sspks,[],'r.')
   else 
   scatter(tphase(ttouchIdx),sspks,[],'k.')
   end
   hold on; plot([0 0],[0 20],'k-.')
   set(gca,'xlim',[-5 5],'ylim',[0 max(sspks)])
   
   

    %FT VEL ANALYSIS DIRECTIONALITY 
    
    ttouchIdx = [find(curr.S_ctk(9,:,:)==1) ;find(curr.S_ctk(12,:,:)==1)];
    ttheta =squeeze(curr.S_ctk(1,:,:));
    tspks = squeeze(curr.R_ntk(1,:,:));
    ftvelTR=nanmean(diff(ttheta(ttouchIdx+(-5:0)),1,2),2);
    FTTRraw = tspks(ttouchIdx + (-25:50));
    
    retT = ftvelTR<0;
    proT = ftvelTR>0;
    
    %ZSCORE = respVect - BLvect - std(BLvect)scalar. 
    meanbase = mean(FTTRraw(:,1:25),2);
    stdbase = nanstd(meanbase);
    postresponses = FTTRraw(:,26:75);
    xresp = mean(postresponses,2);
    zscore = (xresp - meanbase) ./ stdbase;
    
    tFvals{k} = [ftvelTR zscore];
    




%     %FT VEL ANALYSIS STRENGTH OF TOUCH
%     absftv = abs(ftvelTR);
%     [~,aftvelIdx] = sort(absftv);
%     
%     figure(100);subplot(4,8,k);
%     hold on;scatter(absftv,zscore,3,'filled')
%     cor(k) = corr(absftv,zscore)
    
    
end
%num cells sig diff in ret/pro
pthresh = .01;
prcmeans = (prcmeansraw./50)*1000;
keep = prcmeans(prp<pthresh,:);
[~,prpIdx ] = max(keep,[],2);

ut=(length(selectedCells)-length(prpIdx))./length(selectedCells);
pt=sum(prpIdx==2)./length(prpIdx);
rt=sum(prpIdx==1)./length(prpIdx);

figure(90);clf
scatter(prcmeans(prp>pthresh,1),prcmeans(prp>pthresh,2),[],'k','filled')
hold on;scatter(keep(prpIdx==2,1),keep(prpIdx==2,2),[],'b','filled');
hold on;scatter(keep(prpIdx==1,1),keep(prpIdx==1,2),[],'r','filled');
plot([0 100],[0 100],'-.k')
set(gca,'xlim',[0 100],'ylim',[0 100],'xtick',0:25:100,'ytick',0:25:100)
axis square
legend('n.s.','protraction tuned','retraction tuned','location','southeast')
xlabel('mean retraction spk/s (50ms post)');ylabel('mean protraction spk/s (50ms post)');

plotprc = prccounts(prp<pthresh,:);
figure(48);clf
scatter(prccounts(prp>pthresh,1),prccounts(prp>pthresh,2),[],'k','filled')
hold on;scatter(plotprc(prpIdx==2,1),plotprc(prpIdx==2,2),[],'b','filled');
hold on;scatter(plotprc(prpIdx==1,1),plotprc(prpIdx==1,2),[],'r','filled');
hold on; plot([0 2500],[0 2500],'-.k')
set(gca,'xlim',[0 2500],'ylim',[0 2500],'xtick',0:500:3000,'ytick',0:500:3000)
legend('non sig diff','protraction tuned','retraction tuned','location','southeast')
xlabel('ret touch num samples');ylabel('pro touch num samples')
axis square

ptvals = keep(prpIdx==2,:);
rtvals = keep(prpIdx==1,:);
figure(8);clf;
subplot(2,1,1);histogram(ptvals(:,1),0:10:100,'facealpha',1,'facecolor','b')
hold on; histogram(rtvals(:,1),0:10:100,'facealpha',1,'facecolor','r')
set(gca,'ylim',[0 10],'ytick',0:10:10)
ylabel('ret spks')

subplot(2,1,2);histogram(ptvals(:,2),0:10:100,'facealpha',1,'facecolor','b')
hold on; histogram(rtvals(:,2),0:10:100,'facealpha',1,'facecolor','r')
set(gca,'ylim',[0 10],'ytick',0:10:10)
ylabel('pro spks')

%% 
fvs = prp<.05;
[~,tidx ] = max(prcmeans,[],2);
tcs = fvs.*tidx';

utcells = find(tcs==0);
ptcells =find(tcs==2);
rtcells = find(tcs==1);

[pcrc] = numSubplots(length(ptcells));
figure(488);clf
for g = 1:length(ptcells)
    curr = ptcells(g);
    cvals = tFvals{curr};

    ptouches = cvals(:,1)<0;
    
    cvals(:,1) = cvals(:,1).*-1;
    subplot(pcrc(1),pcrc(2),g)
    hold on;scatter(cvals(ptouches,1),cvals(ptouches,2),'.b')
    
    corvals = corr(cvals(ptouches,:));
    
    cv(g) = corvals(2);
    
end
    

[pcrc] = numSubplots(length(rtcells));
figure(489);clf
for g = 1:length(rtcells)
    curr = rtcells(g);
    cvals = tFvals{curr};

    rtouches = cvals(:,1)>0;
    
    subplot(pcrc(1),pcrc(2),g)
    hold on;scatter(cvals(rtouches,1),cvals(rtouches,2),'.r')
    
    
    corvals = corr(cvals(rtouches,:));
    
    fitlm(cvals(rtouches,1),cvals(rtouches,2))
    cvrt(g) = corvals(2);
    
    
end

%best ex of each. 

   cvals = tFvals{rtcells(4)};
   rtouches = cvals(:,1)>0;
   rtvals =cvals(rtouches,:);
   figure(58);clf;subplot(2,1,1)
    hold on;scatter(rtvals(:,1),rtvals(:,2),'.r')
    set(gca,'xlim',[0 2.5],'ylim',[-10 10],'ytick',-10:5:10)
    axis square
    
    
    [~,idx] = sort(rtvals(:,1));
    [coeff, ~ , mu] = polyfit(rtvals(idx,1),rtvals(idx,2),1);
    f = polyval(coeff,0:2,[],mu);
    hold on; plot(0:2,f,'k');
    
    
    
    cvals = tFvals{ptcells(6)};  
    ptouches = cvals(:,1)<0;
    cvals(:,1) = cvals(:,1).*-1;
    subplot(2,1,2)
    hold on;scatter(cvals(ptouches,1),cvals(ptouches,2),'.b')
    set(gca,'xlim',[0 3],'ylim',[-10 10],'ytick',-10:5:10)
    axis square
    
    vals =cvals(ptouches,:);
    [~,idx] = sort(vals(:,1));
    [coeff, ~ , mu] = polyfit(vals(idx,1),vals(idx,2),1);
    f = polyval(coeff,0:5,[],mu);

    hold on; plot(0:5,f,'k');
    
    
    %%
figure(4);clf
histogram(cv,-1:.1:1,'facecolor','b','facealpha',.5)
hold on;histogram(cvrt,-1:.1:1,'facecolor','r','facealpha',.5)
set(gca,'xtick',-1:.5:1,'ytick',0:3:6,'ylim',[0 8])
legend('protraction tuned cells','retraction tuned cells')
title('pre-touch vel vs zscored spk response')
xlabel('correlation coeff');ylabel('num cells')

