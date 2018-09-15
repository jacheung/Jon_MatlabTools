
%% F1 Scatter
colors = {'Gold','DarkTurquoise'};
figure(100);clf

for g = 1:2
    gos = [countsF1{g}(:,1) thetaF1{g}(:,1) timingF1{g}(:,1) motorF1{g}(:,1)];
    nogos = [countsF1{g}(:,2) thetaF1{g}(:,2) timingF1{g}(:,2) motorF1{g}(:,2)];
    gos(isnan(gos))=0;
    nogos(isnan(nogos))=0;
    
    %     gos = [motorF1(:,1)  thetaF1{2}(:,1)];
    %     nogos = [motorF1(:,2)  thetaF1{2}(:,2)];
    %     gos(isnan(gos))=0;
    %     nogos(isnan(nogos))=0;
    
    [~,vgidx(g,:)]=max(gos,[],2);
    [~,vnidx(g,:)]=max(gos,[],2);
    
    figure(100)
    subplot(2,2,1)
    hold on;
    plot(1:4,gos,'-o','Markerfacecolor',rgb(colors{g}),'color',rgb(colors{g}));
    set(gca,'xlim',[.5 4.5],'ylim',[0 1],'ytick',0:.5:1,'xtick',[])
    subplot(2,2,3)
    hold on;
    plot(1:4,nogos,'-o','Markerfacecolor',rgb(colors{g}),'color',rgb(colors{g}));
    set(gca,'xlim',[.5 4.5],'ylim',[0 1],'ytick',0:.5:1,'xtick',[])
    
end


subplot(2,2,2)
hold on;
bar([histcounts(vgidx(1,:),1:5)./length(vgidx) ;histcounts(vgidx(2,:),1:5)./length(vgidx)]);
set(gca,'ylim',[0 1],'ytick',0:.5:1,'xtick',[])
subplot(2,2,4)
hold on;
bar([histcounts(vnidx(1,:),1:4)./length(vnidx) ;histcounts(vnidx(2,:),1:4)./length(vnidx)]);
set(gca,'ylim',[0 1],'ytick',0:.5:1,'xtick',[])


suptitle('go F1s, nogo F1s')
%% RMSE plots
colors = {'Gold','DarkTurquoise'};
% contvals = POPv2.logRMSEDROP{3};
contvals = RMSEpop{2};
% semivals = POPv2.logRMSEDROP{2};
semivals = RMSEpop{1};
figure(480);clf

hold on;
plot(1:4,semivals,'-o','Markerfacecolor',rgb(colors{1}),'Markeredgecolor',rgb(colors{1}),'color',rgb(colors{1}));
set(gca,'xtick',[],'ytick',0:.25:.5,'ylim',[0 .5])
plot(1:4,contvals,'-o','Markerfacecolor',rgb(colors{2}),'Markeredgecolor',rgb(colors{2}),'color',rgb(colors{2}));
set(gca,'xlim',[.5 4.5],'xtick',[],'ytick',0:.25:.5,'ylim',[0 .5])
axis square

[~,sidx]=min(semivals,[],2);
[~,cidx]=min(contvals,[],2);
figure(20);subplot(1,2,1)
histogram(sidx,0:5,'normalization','probability')
set(gca,'xlim',[0 5],'ylim',[0 1])
figure(20);subplot(1,2,2)
histogram(cidx,0:5,'normalization','probability')
set(gca,'xlim',[0 5],'ylim',[0 1])


%% FEATURES WEIGHT

colors = {'Gold','DarkTurquoise'};
popfeatdom = [];
for d = 2
    
    THETAS = uberedfeat{d}(1,:)';
    COUNTS = uberedfeat{d}(2,:)';
    
    featureDom = (abs(THETAS)-abs(COUNTS))./(abs(THETAS)+abs(COUNTS));%multiply by negative 1 so we trend up insteda of down. LOL
    
    popfeatdom = [popfeatdom ;featureDom];
    
end

%EXTENDED TO USE ODDSRATIOS
keep = abs(pasfeat{2}(1:3,:));
normKeep = keep./max(keep);

plotval = normKeep;

xvals = repmat(1:3,1,10);
figure(31);subplot(2,1,1)
scatter(xvals(:),plotval(:),[],[.8 .8 .8],'filled')
set(gca,'xlim',[0 4],'xtick',[],'ylim',[-.1 1.1],'ytick',0:.5:1)
hold on; errorbar(1:3,mean(plotval,2),std(plotval,[],2),'ko')


keep = abs(dtimefeat{2}(1:3,:));
normKeep = keep./max(keep);
plotval = normKeep;
figure(31);subplot(2,1,2)
scatter(xvals(:),plotval(:),[],[.8 .8 .8],'filled')
set(gca,'xlim',[0 4],'xtick',[],'ylim',[-.1 1.1],'ytick',0:.5:1)
hold on; errorbar(1:3,mean(plotval,2),std(plotval,[],2),'ko')


%% True X Predicted
%calculating MCC
colors=jet(5);
for d = 2
%     gcmat = {timingtxp{d},countstxp{d},thetatxp{d},motortxp{d},uberedtxp{d}};
    gcmat = {thetatxp{d},pastxp{d},dtimetxp{d}};
    for k = 1:length(gcmat)
        cmatcurr = gcmat{k};
        for rec = 1:length(cmatcurr)
            cmat = confusionmat(cmatcurr{rec}(:,1),cmatcurr{rec}(:,2));
            TP = cmat(1);FP = cmat(3);
            TN = cmat(4);FN = cmat(2);
            %MCC 
            top  = TP*TN - FP*FN;
            bottom = sqrt((TP+FP) * (TP+FN) * (TN+FP)* (TN+FN));
            mcc{d}(k,rec) = top./bottom;
            if isnan(mcc{d}(k,rec))
                mcc{d}(k,rec) = 0 ;
            end
            %F1 
            precision = TP./(TP+FP);
            recall = TP./(TP+FN);
            F1{d}(k,rec) = 2.*((precision*recall)./(precision+recall));
            
            %Odds Ratio
            or{d}(k,rec) = (TP./FP)./(FN./TN);
        end
    end
end

figure(2);clf
subplot(2,1,1)
bar(mcc{2}([1:3],:)');
set(gca,'ylim',[-.25 1],'ytick',0:.5:1)
legend('timing','counts','theta')
subplot(2,1,2)
bar(mcc{2}([2 4],:)');
legend('theta','motor')
set(gca,'ylim',[-.25 1],'ytick',0:.5:1)


figure(100);clf
scatter(repmat(1:5,1,10)',mcc{2}(:),[],[.8 .8 .8],'filled')
hold on; errorbar(1:5,mean(mcc{d}'),std(mcc{d}'),'ko','linewidth',2)
set(gca,'xlim',[.5 5.5],'ytick',[-.5:.5:1],'xtick',[])


[~,fidx]=max(mcc{2}(1:4,:));
[~,midx]=max(mcc{2}(3:4,:));

figure;bar(histcounts(fidx,1:5)./10)
set(gca,'ylim',[0 1],'ytick',0:.5:1,'xtick',[])

% figure(3);clf
% bar(F1{2}([2 4],:)');
% legend('counts','theta','timing','motor','ubered')







%PLOTTING FOR THETATXP vs PASTXP
xvals = repmat(1:3,10,1);
yvals = mcc{d}';
figure(20);clf
scatter(xvals(:),yvals(:),[],[.8 .8 .8],'filled')
hold on;errorbar(1:3,mean(mcc{d},2),std(mcc{d}'),'ko')
set(gca,'xlim',[0 4],'xtick',1:3,'xticklabel',{'angle','PAM','timeDecomp'},'ylim',[0 1],'ytick',0:.5:1)
ylabel('MCC')













