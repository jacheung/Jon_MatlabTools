%% quantifying late touch features

touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
[rc] = numSubplots(length(selectedCells));

%%
close all;
clearvars -except U selectedCells rc
rsqThresh = .5;

for k = 1:length(selectedCells)
    curr = U{selectedCells(k)};
    motors = curr.meta.motorPosition';
    theta = curr.S_ctk(1,:,:); 
    touchIdx = [find(curr.S_ctk(9,:,:)==1) ;find(curr.S_ctk(12,:,:)==1)];
    touchIdxoff = [find(curr.S_ctk(10,:,:)==1); find(curr.S_ctk(13,:,:)==1)];
    dkap = curr.S_ctk(6,:,:);
    spks = squeeze(curr.R_ntk(:,:,:));
    
    %touch Duration
    durSpk = zeros(length(touchIdx),1);
    durMotor = durSpk;
    dkapMax = durSpk;
    dkapInteg = dkapMax;
    sta = zeros(length(touchIdx),11);
    for g = 1:length(touchIdx)
        curIdx = touchIdx(g):touchIdxoff(g);
        rawSpks = spks(touchIdx(g):touchIdxoff(g));
        durMotor(g) = motors(ceil(touchIdx(g)./curr.t));
        durSpk(g) = nansum(rawSpks);
        dkapMax(g) = max(dkap(curIdx));
        dkapInteg(g) = trapz(dkap(touchIdx(g):touchIdxoff(g)));
    end
    touchDur = touchIdxoff-touchIdx;
        %regressing out effects of motor position on duration using a linear
    %model
    durMdl = fitlm(durMotor,touchDur);
    durResid = durMdl.Residuals.Raw;
    
    %     mdl = fitlm(touchDur,durSpk);
    mdl = fitlm(durResid,durSpk);
    pvals = mdl.Coefficients.pValue;
    if pvals(2)<.01
        rsq(k) = mdl.Rsquared.Ordinary;
        coeff = mdl.Coefficients.Estimate;
        plotV(k,:) = coeff(1) + coeff(2)*(1:10:500);
    end
    
    
    dkapVal = dkapMax;
    durVal = durResid; % durResid or touchDur

    [sortedDur,idx] = sort(durVal);
    [sorteddk,dkidx] = sort(dkapVal);
    ninetyptouchdur = sortedDur(round(length(idx)*.90));
    ninetydkap = sorteddk(round(length(dkidx)*.90));
    fsdk = flipud(sorteddk);
    ninetydkapmin = fsdk(round(length(dkidx)*.90));
    
    % INDIVIDUAL FEATURE PLOTS (durationXspks and max(dkap)Xspks)
    %     figure(88);subplot(rc(1),rc(2),k)
    %     if rsq(k)>rsqThresh
    %         figure(88);subplot(rc(1),rc(2),k);scatter(touchDur,durSpk,'r.')
    %         set(gca,'xlim',[0 ninetyptouchdur])
    %         corVal = corr(touchDur,durSpk);
    %         title(num2str(corVal))
    %
    %         figure(91);subplot(rc(1),rc(2),k);scatter(dkapMax,durSpk,'r.')
    %         set(gca,'xlim',[0 ninetydkap])
    %         figure(90);subplot(rc(1),rc(2),k);scatter(touchDur,(durSpk./touchDur).*1000,[],'r.')
    %     else
    %         figure(88);subplot(rc(1),rc(2),k);scatter(touchDur,durSpk,[],'.','markeredgecolor',[.6 .6 .6])
    %          set(gca,'xlim',[0 ninetyptouchdur])
    %          corVal = corr(touchDur,durSpk);
    %             title(num2str(corVal))
    %
    %         figure(90);subplot(rc(1),rc(2),k);scatter(touchDur,(durSpk./touchDur).*1000,[],'.','markeredgecolor',[.6 .6 .6])
    %     end
    %     axis square
    
  
    mdl=fitlm([dkapVal durVal],durSpk);
    mrrsq = mdl.Rsquared.Ordinary;
    cortmp=corr([dkapVal durVal durSpk]);
    srrsq = (cortmp([1 2],3).^2).*(mdl.Coefficients.pValue(2:end)<.01);
    prsqd(k,:) = mrrsq-srrsq;
    jointrsqd(k,1) = mrrsq-sum(prsqd(k,:),2);
    
    
    [sds,sdsidx] = sort(durSpk);
    ninesds = sds(round(length(sdsidx)*.90));
    durSpk(durSpk>ninesds) = ninesds;
    
    figure(89);subplot(rc(1),rc(2),k)
    scatter(durVal,dkapVal,10,durSpk,'filled')
    %     if rsq(k)>rsqThresh
    %        scatter(touchDur,dkapVal,10,durSpk,'filled')
    %     else
    %        scatter(touchDur,dkapVal,[],'.k')
    %     end
    %      set(gca,'xlim',[0 ninetyptouchdur],'ylim',[0 ninetydkap])
    set(gca,'xlim',[0 ninetyptouchdur],'ylim',[ninetydkapmin ninetydkap])
    
    if k == 14
        figure(43);clf
        scatter3(dkapVal,durVal,durSpk)
        hold on;scatter3(dkapVal,durVal,mdl.predict)
    end
    
end
rsqdVals = [prsqd jointrsqd];
totalRsqd = sum(rsqdVals,2);
figure(9);clf;subplot(2,1,1);histogram(totalRsqd,[0:.1:1],'facecolor','k','facealpha',1)
set(gca,'xtick',0:.25:1,'ylim',[0 8])

propCont = rsqdVals./totalRsqd;
selProps = propCont(totalRsqd>rsqThresh,:);
propmeanTop = mean(selProps);
propstdTop = std(selProps);
propsemTop = std(selProps)./size(selProps,1);
figure(9);subplot(2,1,2);filex_barwitherr(propsemTop,propmeanTop,'k')
set(gca,'xticklabel',{'duration of touch','dkap','joint'},'ylim',[0 1],'ytick',0:.5:1)

[p,~,stats]=anova1(propCont,[],'off');
    vals =multcompare(stats,[],'off');

% propmean = mean(propCont);
% propsem = std(propCont)./size(propCont,1);
% figure(9);subplot(2,1,2);filex_barwitherr(propsem,propmean,'k')
% set(gca,'xticklabel',{'duration of touch','dkap','joint'},'ylim',[0 1],'ytick',0:.5:1)
% figure(8880);clf
% selCells = find(rsq>rsqThresh);
% plot(1:10:500,plotV(rsq>rsqThresh,:),'r-')
% hold on; plot(1:10:500,plotV(rsq<=rsqThresh,:),'k-','color',[.6 .6 .6])
% hold on; plot(1:10:500,plotV([14 24],:),'g-')
% set(gca,'xtick',0:100:500,'xlim',[0 250],'ylim',[0 20],'ytick',0:5:20)
%
% propPsig = numel(rsq)./numel(selectedCells);
% propRsqSig = sum(rsq>rsqThresh)./numel(selectedCells);




