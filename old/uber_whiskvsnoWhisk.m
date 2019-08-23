%% SIMPLE whisking analysis Amp>5 vs amp<5
clearvars -except U
for k = 1:length(U)
    curr = U{k};
    masks = assist_touchmasks(curr);
    
    
    amps = squeeze(curr.S_ctk(3,:,:)).*masks.touch;
    spks = squeeze(curr.R_ntk(1,:,:)).*masks.touch;
    
    
    nw=spks(amps<5);
    w=spks(amps>5);
    
    [~,p(k)] = ttest2(nw,w);
    
    
    nwfr(k) = mean(nw)*1000;
    wfr(k) = mean(w)*1000;
    
end

whiskTunedcells = find(p<.01);
unsig = find(p>.01);

figure(9);clf
nwc = nan(length(U),1);
wc = nan(length(U),1);

for g = 1:length(whiskTunedcells)
    curr = whiskTunedcells(g);
    if nwfr(curr) > wfr(curr)
        hold on;scatter(nwfr(curr),wfr(curr),[],'filled','r')
        nwc(curr) = 1;
    elseif nwfr(curr)<wfr(curr)
        hold on;scatter(nwfr(curr),wfr(curr),[],'filled','b')
        wc(curr) = 1;
    end
end
scatter(nwfr(unsig),wfr(unsig),[],'filled','markerfacecolor',[.7 .7 .7])
set(gca,'xlim',[0 30],'ylim',[0 30],'xtick',0:10:30,'ytick',0:10:30)
hold on; plot([0 30],[0 30],'-.k')
xlabel('no whisk spks/s');ylabel('whisk spks/s')
axis square

vals = [nwfr' wfr'];
nwcells = vals(find(nwc==1),:);
wcells = vals(find(wc==1),:);


figure(8);clf;
subplot(2,1,1);histogram(wcells(:,1),0:2.5:30,'facealpha',1,'facecolor','b')
hold on; histogram(nwcells(:,1),0:2.5:30,'facealpha',1,'facecolor','r')
set(gca,'ylim',[0 20],'ytick',0:20:20)
ylabel('nw spks')

subplot(2,1,2);histogram(wcells(:,2),0:2.5:30,'facealpha',1,'facecolor','b')
hold on; histogram(nwcells(:,2),0:2.5:30,'facealpha',1,'facecolor','r')
set(gca,'ylim',[0 20],'ytick',0:20:20)
ylabel('w spks')


%prop tuned blocking out touch periods.
length(whiskTunedcells)./length(U)

wccells = find(wc==1);
nwcells = find(nwc==1);

%% WHISKING ON CELL ANALYSIS PRO/RET analysis
wccells = find(wc==1);
rVSp = zeros(length(wccells),2);

pthresh = .01;

figure(100);clf
for d = 1:length(wccells)
    curr = U{wccells(d)};
    masks = assist_touchmasks(curr);
    
    amps = squeeze(curr.S_ctk(3,:,:)).*masks.touch;
    phase = squeeze(curr.S_ctk(5,:,:)).*masks.touch;
    spks = squeeze(curr.R_ntk(1,:,:)).*masks.touch;
    
    retWidx = intersect(find(amps>5),find(phase>0));
    proWidx = intersect(find(amps>5),find(phase<0));
    
    [~,rppval(d)] = ttest2(spks(retWidx),spks(proWidx));
    
    rVSp(d,:) = [mean(spks(retWidx))*1000 mean(spks(proWidx))*1000];
    
    [~,rORp]  = max(rVSp(d,:));
    
    figure(100);
    if rppval(d)<=pthresh
        if rORp ==2
            hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','b')
        elseif rORp ==1
            hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','r')
        end
    elseif rppval(d)>pthresh
        hold on; scatter(rVSp(d,1),rVSp(d,2),[],'filled','markerfacecolor',[.7 .7 .7])
    end
    
end

hold on; plot([0 30],[0 30],'-.k')
axis square
xlabel('retraction spks/s'); ylabel('protraction spks/s')
set(gca,'xtick',0:10:30,'ytick',0:10:30)

propWhiskDtuned = sum(rppval<=pthresh)./numel(rppval);
dirCells = (rppval<=pthresh);
rVSp(dirCells,:)

%% WHISKING cells tuned to peak or trough
wccells = find(wc==1);

for d = 1:length(wccells)
    curr = U{wccells(d)};
    findMaxMinProtraction(curr,5)
end

%% WHISKING ON CELLS TUNED TO PAM or Theta
wccells = find(wc==1);
nwcells = find(nwc==1);
untuned = find(nansum([wc nwc],2)==0);
all = 1:length(U); 
tsp = numSubplots(length(U));
v = 5
sw = 50;
peakBins = linspace(-pi,pi,7);
figure(548);clf
% datasample(1:length(U),1)
sel = all;
for d = 1:length(sel)
    rec =sel(d);
    figure(548);clf
    
    
    
    [objmask]= assist_touchmasks(U{rec});
    mask = objmask.touch;
    amp = squeeze(U{rec}.S_ctk(3,:,:));
    ampmask = nan(size(amp));
    ampmask(amp>3)=1;
    %         var = squeeze(U{rec}.S_ctk(v,:,:)).*mask.*ampmask;
    %         spks = squeeze(U{rec}.R_ntk(:,:,:)).*mask.*ampmask;
    %
    rawvar = squeeze(U{rec}.S_ctk(v,:,:));
    peakWhiskIdx = intersect(find(rawvar<3.1+.05),find(rawvar>3.1-.05));
    keep = find(diff([0; peakWhiskIdx])>10);
    selpos = peakWhiskIdx(keep);
    
    peakWhiskIdx = intersect(find(rawvar<-3.1+.05),find(rawvar>-3.1-.05));
    keep = find(diff([0; peakWhiskIdx])>10);
    selneg = peakWhiskIdx(keep);
    
    tmp = (selpos - selneg');
    tidx = find(isnan(mask)==1);
    ps = nan(size(tmp,1),2);
    for b = 1:size(tmp,1)
        currTmp = tmp(b,:);
        tmpVal = min(currTmp(currTmp>0));
        if ~isempty(tmpVal)
            ps(b,:) =[selneg(find(currTmp == tmpVal)) selpos(b)];
        end
        pwindow = ps(b,1):ps(b,2);
        
        if ~isempty(intersect(pwindow,tidx))
            ps(b,:) = [nan nan];
        end
    end
    
    ps(isnan(ps(:,1)),:) = [];
    completeMask = nan(size(rawvar));
    for u = 1:size(ps,1)
        completeMask(ps(u,1):ps(u,2))=1;
    end
    
    var = squeeze(U{rec}.S_ctk(v,:,:)).*completeMask.*ampmask;
    spks = squeeze(U{rec}.R_ntk(:,:,:)).*completeMask.*ampmask;
    
    
%               var = squeeze(U{rec}.S_ctk(v,:,:)).*mask.*ampmask;
%             spks = squeeze(U{rec}.R_ntk(:,:,:)).*mask.*ampmask;
    
    
    %          figure(58);clf;
    %         plot(rawvar(:,1),'.k')
    %         hold on;plot(selneg(selneg<4000),rawvar(selneg(selneg<4000),1),'ro')
    %         hold on;plot(selpos(selpos<4000),rawvar(selpos(selpos<4000),1),'go')
    
    for g = 1:length(peakBins)
        
        peakWhiskIdx = intersect(find(var<peakBins(g)+.05),find(var>peakBins(g)+-.05));
        keep = find(diff([0; peakWhiskIdx])>5);
        selPhase = peakWhiskIdx(keep);
        
        clear spksPhase
        for b = 1:length(selPhase)
            if (selPhase(b)-sw)>0 & (selPhase(b)+sw)<=numel(spks)
                spksPhase(b,:) = spks(selPhase(b)-sw : selPhase(b)+sw);
            end
        end
        
        
        

%         figure(548);
%         subplot(7,1,g)
%         %     subplot(tsp(1),tsp(2),d)
%         bar(-sw:sw,nanmean(spksPhase)*1000,1,'facecolor',[.6 .6 .6])
%         hold on; plot(-sw:sw,smooth(nanmean(spksPhase)*1000,10),'g','linewidth',2)
%         set(gca,'xlim',[-sw sw],'xtick',linspace(-sw,sw,5))
%         if g == 7
%             xlabel('time from peak phase (ms)');ylabel('firing rate (hz)')
%         end
%         %
        
        
        % CALCULATION OF TUNING 
        [p,tbl,stats] = anova1(spksPhase,[],'off');
        frs = smooth(nanmean(spksPhase)*1000,10);
        trimfrs = frs(5:end-5);
        [~,idx] = max(trimfrs);
        modidx = max(trimfrs)./mean(trimfrs);
        
%         zresp = (modidx-mean(datasample(trimfrs,50)./mean(trimfrs)))./std(datasample(trimfrs,50)./mean(trimfrs));
        

        frp(g,:) = [p modidx idx-sw+5];
        
        rawSpk(g,:) = smooth(nanmean(spksPhase)*1000,10); 
        
    end
    ptune(d).varNames = {'ANOVA p-Val','modIdx','peakFR'};
    ptune(d).stats = frp;
    ptune(d).fr = rawSpk; 
    %     [~,midx] = max(frp(:,1));
%         print(figure(548),'-dpng',['C:\Users\jacheung\Dropbox\LocationCode\Figures\dump\v2cellNum_' num2str(rec)])
    
end

%%
sel = all;
tuned = nan(length(sel),size(frp,2)+1);
nsp = numSubplots(numel(sel));
pvalThresh=.05;
close all
for k = 1:length(sel)
    curr = sel(k);
    
    vals = [ptune(curr).stats (1:length(ptune(curr).stats))'];
    
    thresd = double(vals(:,1)<pvalThresh) .* double(nanmean(ptune(curr).fr,2)>1);

    sigvals = vals(logical(thresd),:);
    sigvals = sigvals(sigvals(:,3)>0,:); %choose only ones with positive peaks
    
    
    if ~isempty(sigvals)
        
        for o = 1:size(sigvals,1)
            figure(843);subplot(nsp(1),nsp(2),k)
            pspks = normalize_var(ptune(curr).fr(sigvals(o,end),:),0,1);
            hold on;plot(-sw:sw,pspks,'color',[.6 .6 .6],'linewidth',1)
        end
        
        [~,idx ] = max(sigvals(:,2));
        tuned(k,:) = sigvals(idx,:);
        
        
        chosen = normalize_var(ptune(curr).fr(sigvals(idx,end),:),0,1);
        
        if sum(curr == wccells)
            hold on; plot(-sw:sw,chosen,'b','linewidth',2)
        set(gca,'xlim',[-sw+5 sw-5],'xtick',linspace(-sw,sw,5))
            
        figure(69);hold on; scatter(sigvals(idx,end),sigvals(idx,2),'b','filled')
        elseif sum(curr == nwcells)
            hold on; plot(-sw:sw,chosen,'r','linewidth',2)
            set(gca,'xlim',[-sw+5 sw-5],'xtick',linspace(-sw,sw,5))
            figure(69);hold on; scatter(sigvals(idx,end),sigvals(idx,2),'r','filled')
        else
            hold on; plot(-sw:sw,chosen,'k','linewidth',2)
            set(gca,'xlim',[-sw+5 sw-5],'xtick',linspace(-sw,sw,5))
            figure(69);hold on; scatter(sigvals(idx,end),sigvals(idx,2),'k','filled')
        end
        
        figure(580);subplot(1,7,sigvals(idx,end))
        hold on; plot(-sw:sw,chosen,'color',color(sigvals(idx,end),:),'linewidth',1)
        set(gca,'xlim',[-sw+5 sw-5],'xtick',linspace(-sw,sw,5),'ylim',[0 1])
        
        
    end
end

figure(9);clf;
cts = histcounts(tuned(:,end),.8:8);
bar(peakBins,cts)
set(gca,'ylim',[0 5],'xlim',[-pi-.5 pi+.5],'xtick',peakBins,'xticklabel',peakBins,'ytick',0:5)


set(gca,'ylim',[0 5],'xlim',[0 8],'xtick',1:7,'xticklabel',peakBins)