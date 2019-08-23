

clearvars -except pop U ocells
% touchCells = touchCell(U);
% selectedCells = find(touchCells==1);
close all
% selectedCells=1:length(U);
[ns] = numSubplots(length(ocells));
thetabounds = [-100:2.5:100];

for d = 1:length(ocells)
    
    thetaBounds = -100:2.5:100;
    currCell = ocells(d);
     window = [-50:50]; %ms to look around touch
    
    [varspikes, ~,~] = assist_varAtTouch(U{currCell},window);

    %assigning variables
    vars = varspikes(:,1:6);
    spksraw = varspikes(:,7:end);
        
    postresponses = spksraw(:,55:75);
    baseresponses = spksraw(:,25:51); 
    % ONE THING TO CHANGE IS TO CAPTURE THE TOUCH RESPONSE FOR EACH NEURON.
    % SINCE EACH ONE IS DIFFERENT. WE CAN USE THE TOUCHCELL FUNCTION TO
    % FIND THIS WHICH BIN IS SIG AND USE THOSE BINS TO DO THIS COMPARISON. 
    
    tmpmb  = mean(baseresponses,2);
    tmpmp = mean(postresponses,2);
    
    [smb, ~ ,~]=binslin(vars(:,1),[tmpmb tmpmp],'equalE',numel(thetabounds)+1,-100,100);
     
    % FOR BINNING DURING VISUALIZATION 
    samps = cellfun(@numel,smb);
    selBins = find(samps>sum(samps)./100);
    
    %% MODEL USING MAP+BL
    prcts = sum(postresponses,2);
    brcts = sum(baseresponses,2);
    pam = vars(:,3:5);
    theta = vars(:,1); 
    
  
    
    selRows = ~sum(isnan([brcts pam]),2);
    thetaPlot = theta(selRows); 
    
    DmatX = [brcts pam];

    DmatX = DmatX(selRows,:);
%     DmatX = DmatX-mean(DmatX);
        DmatX = filex_whiten(DmatX); 
    DmatY = prcts(selRows,:); 
    rando = randperm(length(DmatY));
    
    qtrs = floor(length(DmatY)/4);
    clear samp;
    samp(1,:) = rando(1:qtrs);
    samp(2,:) = rando(qtrs+1:qtrs*2);
    samp(3,:) = rando(qtrs*2+1:qtrs*3);
    samp(4,:) = rando(qtrs*3+1:qtrs*4);

    pred = nan(size(samp));
    
    for k = 1:size(samp,1)
        test = samp(k,:); 
        train = samp(setdiff(1:4,k),:); 
        [mdl,~,mdlStats] = glmfit(DmatX(train(:),:),DmatY(train(:)),'poisson');
        testpred = glmval(mdl,DmatX(test,:),'log');
%         figure(80);subplot(2,2,k)
        pred(k,:) = testpred;
        
        mdlBCoeff{d}(:,k) = mdlStats.beta;
        mdlBp{d}(:,k) = mdlStats.p;
        
    end
   
     [ppred, ~ ,~]=binslin(thetaPlot(samp(:)),[DmatY(samp(:)) pred(:)],'equalE',numel(thetabounds)+1,-100,100);
     
     sortedPred = cell2mat(cellfun(@(x) mean(x,1),ppred,'uniformoutput',0));
    
      figure(850);subplot(ns(1),ns(2),d)
%       figure(850);clf
     scatter(thetabounds(selBins),sortedPred(selBins,1)*1000,[],'k','filled')
     hold on;  scatter(thetabounds(selBins),sortedPred(selBins,2)*1000,[],'r','filled')
        [evalmod] = fitlm(sortedPred(selBins,2),sortedPred(selBins,1));
        rsqmdl(d) = evalmod.Rsquared.Ordinary;
        set(gca,'xlim',[min(thetabounds(selBins))-5 max(thetabounds(selBins))+5])
        title(['Rsqd = ' num2str(evalmod.Rsquared.Ordinary)])
        
        
        
      frmeanstd(d,:) = [mean(DmatY) var(DmatY)];
          
%     %% MODEL THAT USES BOTH AVERAGE MAP AND AVERAGE BL TO PREDICT 
%      [sorted, ~ ,~]=binslin(vars(:,1),[tmpmb tmpmp vars],'equalE',numel(thetabounds)+1,-100,100);
%      nums = cellfun(@(x) size(x,1),sorted);
%      thetavect=[];
%      for k = 1:length(thetabounds)
%          if nums(k)>0
%          thetavect = [thetavect; ones(nums(k),1)*thetabounds(k)];
%          end
%      end
%      
%      tmp =  cellfun(@(x) nanmean(x,1),sorted,'uniformoutput',0);
%      tmp2= cell2mat(tmp(selBins));
%      DmatX = filex_whiten(tmp2(:,[1 5 6 7]));
%      DmatX2 = filex_whiten(tmp2(:,[1 3])); 
%      
%      fitmdl = fitlm(DmatX,tmp2(:,2));
%      fitmdl2= fitlm(DmatX2,tmp2(:,2));
%      
%      figure(850);subplot(ns(1),ns(2),d)
%      scatter(thetabounds(selBins),tmp2(:,2)*1000,[],'k','filled')
%      hold on; scatter(thetabounds(selBins),fitmdl.predict*1000,[],'r','filled')
%      hold on; scatter(thetabounds(selBins),fitmdl2.predict*1000,[],'b','filled')
%      set(gca,'xlim',[min(thetabounds(selBins))-5 max(thetabounds(selBins))+5])
%      title(['Rsqd = ' num2str(fitmdl.Rsquared.Ordinary)])
%      
%     coeffVals(:,d)=fitmdl.Coefficients.Estimate(2:end);
%     coeffpv(:,d) = fitmdl.Coefficients.pValue(2:end);
%     rsqd(1,d) =  fitmdl.Rsquared.Adjusted;
    
    
    
     
    %% testing to see if if baselines are significantly different
    braw = nan(length(smb),2000);
    for b = 1:length(smb)
        if ismember(b,selBins)
            currvals = smb{b};
            if ~isempty(smb{b})
                braw(b,1:length(currvals)) = currvals(:,1)';
            end
        end
    end
    [p,~,stats]=anova1(braw',[],'off');
    vals =multcompare(stats,[],'off');
    sigs =  vals(vals(:,end)<.01,:);
     
    zs = nanmean(braw,2);
    zs(isnan(zs))=[];
    sigs =  vals(vals(:,end)<.01,:);
    
    clear peaksig
    [peakFR,maxidx]=max(zs);
    midxsel = find(sum(sigs(:,[1 2])==maxidx,2)==1);
    
    if sum(midxsel>0)
        peakBase(d,:) = [peakFR*1000 mean(mean(baseresponses,2))*1000];
    else
        peakBase(d,:) = [nan mean(mean(baseresponses,2))*1000]; 
    end
    
    basemean = cell2mat(cellfun(@(x) mean(x,1),smb,'uniformoutput',0));

    figure(50);subplot(ns(1),ns(2),d)
    allspks = [basemean(selBins,1)*1000 ; basemean(selBins,2)*1000];
    scatter(basemean(selBins,1)*1000,basemean(selBins,2)*1000,'k','filled')
    set(gca,'xlim',[min(allspks) max(allspks)],'ylim',[min(allspks) max(allspks)])  
    hold on; plot([min(allspks) max(allspks)],[min(allspks) max(allspks)],'-.k')
%     ylabel('touch response')
%     xlabel('baseline response')
    linepred = fitlm(basemean(selBins,1)*1000,basemean(selBins,2)*1000);
    hold on;plot(basemean(selBins,1)*1000,linepred.predict,'-r')
    
    title(['Rsqd = ' num2str(linepred.Rsquared.Ordinary)])
    
    if linepred.Coefficients.pValue(2)<.05
   text(median(allspks),median(allspks),['corr = ' num2str(corr(basemean(selBins,1)*1000 , basemean(selBins,2)*1000))]);
    end
    
    figure(6580);subplot(ns(1),ns(2),d)
    scatter(thetaBounds(selBins),basemean(selBins,1)*1000,[],[.6 .6 .6],'filled');
    hold on;scatter(thetaBounds(selBins),basemean(selBins,2)*1000,[],'k','filled');


    
end

figure(6580)
suptitle('Does baseline response(gray) predict post touch response (black)')

figure(50)
suptitle('Relationship between baseline and postresponse')

figure(850)
suptitle('How well can baseline + MAP at touch predict post touch firing rate? raw(black) predict(red)')

%sig modulated baselines
blModNaive = sum(~isnan(peakBase(1:9,1)))./ numel(~isnan(peakBase(1:9,1)));
blModExp = sum(~isnan(peakBase(10:end,1)))./ numel(~isnan(peakBase(10:end,1)));
figure(580);clf
bar(1,blModNaive,'facecolor',[.6 .6 .6])
hold on; bar(2,blModExp,'k')
set(gca,'ylim',[0 .75],'ytick',0:.25:1,'xlim',[0 3],'xtick',1:2,'xticklabel',{'naive','expert'})
ylabel('Proportion of neurons with baselines sig. modulated')

%FOR LINEAR MODEL (zscore|params) ABSOLUTE NORMALIZED WEIGHTS (did this because we had negative effects.
%This way shows us the ABSOLUTE WEIGHT. 
% sigcoeff = coeffVals.*double(coeffpv<.05);
% normAbsCoeff = normalize_var(abs(sigcoeff),0,1);
% meanCoeff = nanmean(normAbsCoeff,2);
% sem=nanstd(normAbsCoeff,[],2)./sqrt(size(normAbsCoeff,2));
% [~,~,stats]=anova1(normAbsCoeff',[],'off');
% multcompare(stats)
% 
% figure(51560);clf
% x= repmat(1:size(normAbsCoeff,1),1,size(normAbsCoeff,2));
% hold on; scatter(x(:),normAbsCoeff(:),[],[.6 .6 .6],'filled')
% % hold on; plot([1 4],[0 0],'-.k')
% errorbar(1:4,meanCoeff,sem,'-ko','linewidth',2)
% set(gca,'xtick',1:4,'ytick',-1:.5:1,'xticklabel',{'baseline','amp','mp','phase'},'ylim',[0 1.1],'xlim',[0 5])
% ylabel('normalized absolute weights')
% title('Absolute weight of each feature in predicting post touch response')

%FOR LNP MODEL 
groupCoeffs = cell2mat(cellfun(@(x) mean(x,2),mdlBCoeff,'uniformoutput',0));

%USed to nan out n.s. pvals for weights. 
for n = 1:length(ocells)
    sigvals = mdlBCoeff{n}.*(mdlBp{n}<.05);
    sigvals(sigvals==0)=0; 
    wt(:,n) = nanmean(sigvals,2);
end
% poissEffect = exp(groupCoeffs(2:end,:));
poissEffect = exp(wt(2:end,:));
mpe = nanmean(poissEffect,2);
mpsem = nanstd(poissEffect,[],2)./sqrt(size(poissEffect,2));
[~,~,stats]=anova1(poissEffect',[],'off');
multcompare(stats,[],'off')

figure(51560);clf
x= repmat(1:size(poissEffect,1),1,size(poissEffect,2));
hold on; scatter(x(:),poissEffect(:),[],[.8 .8 .8],'filled')
hold on; plot(x(:),poissEffect(:),'color',[.8 .8 .8])
hold on; plot([1 4],[1 1],'-.k')
errorbar(1:4,mpe,mpsem,'-ko','linewidth',2)
set(gca,'xtick',1:4,'ytick',0:.25:2,'xticklabel',{'baseline','amp','mp','phase'},'ylim',[.7 1.8],'xlim',[0 5])
ylabel('% fr change per one std increase')







figure(57);clf
histogram(rsqmdl,0:.1:1,'facecolor','k')
set(gca,'xtick',0:.25:1,'ylim',[0 12],'ytick',0:2:10)
xlabel('RSquared');ylabel('number of neurons')
