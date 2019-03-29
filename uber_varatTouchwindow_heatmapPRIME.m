
%%
% clearvars -except U

for rec=1:length(U)
    countThresh = 10; %min touch in each bin to be considered
    gaussFilt = [1]; %sigma filter for imgaussfilt (can be single value for filtering across rows and columns or two values [x y] for filter of row/column
    window = [-25:50]; %ms to look around touch
    
    [varspikes, preDvarspikes,postDvarspikes] = assist_varAtTouch(U{rec},window);
    % First 6 columns will be values for the variables
    % 1) THETA
    % 2) PRE TOUCH VELOCITY
    % 3) AMP
    % 4) SETPOINT
    % 5) PHASE
    % 6) MAX KAPPA
    % Last columns will be the spikes around your given window
    
    
    % Run this line below if you want to elim all protraction touches and
    % just look at ret touches.
    %     varspikes(varspikes(:,5)<0,:)=[];
    
    
    
    %     comp={varspikes, preDvarspikes, postDvarspikes};
    comp={varspikes};
    for g = 1:length(comp)
        %%
        %Plot theta at touch
        fields = [      {'theta'}       {'velocity'}      {'amplitude'}    {'setpoint'}          {'phase'}           {'kappa'}];
        V.bounds = [{[-100:2:100]} {[-9750:125:9750]} {[-99.5:2.5:99.5]} {[-99.5:2.5:99.5]} {linspace(pi*-1,pi,13)} {[-.95:.05:.95]}];
        for d = [1 2 3 4 5 6] %for variables 1:6
            
            if d == 2
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalE',numel(V.bounds{d})+1,-10000,10000);
            elseif d == 5
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalX',numel(V.bounds{d})+1);
            elseif d == 6
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalE',numel(V.bounds{d})+1,-1,1);
            else
                [sorted, sortedBy ,binBounds]=binslin(comp{g}(:,d),comp{g}(:,7:end),'equalE',numel(V.bounds{d})+1,-100,100);
            end
            
            binrange = V.bounds{d};
            
            % Trimming unused ranges
            trims=[binrange' cell2mat(cellfun(@mean,sortedBy,'Uniformoutput',0))];
            indexHasTheta = ~isnan(trims(:,2));
            trims = trims(indexHasTheta, :);
            counttmp=cell2mat(cellfun(@size,sorted,'uniformoutput',0));
            
            % Populating fields with values
            V.(fields{d}).counts =counttmp(indexHasTheta==1,1);
            V.(fields{d}).range = trims(:,1);
            V.(fields{d}).spikes=cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
            V.(fields{d}).spikes=V.(fields{d}).spikes(indexHasTheta==1,:);
            V.(fields{d}).raw = sorted(indexHasTheta==1);
            
            %Trimming bins with touchcounts below thresholds
            mintouches=find(V.(fields{d}).counts<countThresh);
            V.(fields{d}).counts(mintouches,:)=[];
            V.(fields{d}).spikes(mintouches,:)=[];
            V.(fields{d}).range(mintouches,:)=[];
            V.(fields{d}).raw(mintouches,:) = [];
            
            %Plotting features
            %             figure(30+rec);subplot(2,3,d);
            %
            %
            %
            % %             imagesc(imgaussfilt(V.(fields{d}).spikes,gaussFilt,'padding','replicate')); %GAUSS PLOTTING
            %             imagesc(V.(fields{d}).spikes) %OLD PLOTTING
            %
            %             colormap(gca,parula);
            %             set(gca,'Ydir','normal','ytick',(1:length(V.(fields{d}).range)),'yticklabel',[V.(fields{d}).range],...
            %                 'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
            %             for k=1:size(V.(fields{d}).counts,1)
            %                 text(20,k,num2str(V.(fields{d}).counts(k)),'FontSize',8,'Color','white')
            %             end
            %             hold on;plot([sum(window<0) sum(window<0)],[length(V.(fields{d}).range) 0],'w:')
            %             axis('square')
            %             xlabel('time from touch onset (ms)')
            %             title([fields{d}])
            %
            
            pop{rec}=V;
            
            
        end
        
    end
end

%% POPULATION PLOTTING TUNING

% PLOT TOUCH CELLS+
window = [-25:50];
gaussFilt = [1];
fields = [    {'theta'}  ];

% touchCells = touchCell(U,2.5);
touchCells = touchCell(U);
selectedCells = find(touchCells==1);
[rc] = numSubplots(length(selectedCells));

plotrow = rc(1);
plotcolumn = rc(2);

figure(333);clf
figure(322);clf
for d = 1:length(selectedCells)
    
    currCell = selectedCells(d);
    array = U{currCell};
    
    %toss touches out that are >10 degree difference from nearest bin. 
    toKeep = diff(pop{currCell}.theta.range)<10; 
    if ~isempty(find(toKeep==0))
        toKeep(find(toKeep==0):end) = 0;
    end
    
    toKeep = logical(ones(length(pop{currCell}.theta.range),1)); 
    
    rawHeat = pop{currCell}.(fields{1}).spikes(toKeep,:);
    smoothedHeat = imgaussfilt(pop{currCell}.(fields{1}).spikes(toKeep,:),gaussFilt,'padding','replicate'); 
    heatToPlot = smoothedHeat; 
    
    figure(322);subplot(plotrow,plotcolumn,d);
    imagesc(heatToPlot);
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',1:sum(toKeep),'yticklabel',[pop{currCell}.(fields{1}).range(toKeep)],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{currCell}.(fields{1}).counts(toKeep),1)
        text(20,k,num2str(pop{currCell}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{currCell}.(fields{1}).range) 0],'w:')
    axis square
    caxis([0 prctile(heatToPlot(:),99)]) %99th percentile cut off so that outliers don't alias signal
    sampledRange(d) = max([pop{currCell}.(fields{1}).range(toKeep)]) - min([pop{currCell}.(fields{1}).range(toKeep)]);
    
    selIdx = [find(window==array.meta.responseWindow(1)):find(window==array.meta.responseWindow(2))];
    hold on; plot([selIdx(1) selIdx(1)],[0 30],'-.w')
    hold on; plot([selIdx(end) selIdx(end)],[0 30],'-.w')

end


%% ZSCORE CONVERSION and tuning width
clearvars -except pop U ocellidx
window = [-25:50];
touchCells = touchCell(U);
selectedCells = find(touchCells==1);
% selectedCells=1:length(U);
[rc] = numSubplots(length(selectedCells));
fields = [    {'theta'}  ];
% thetabounds = linspace(pi*-1,pi,13);
thetabounds = [-100:2:100];
plotrow = rc(1);
plotcolumn = rc(2);
ytune = nan(length(thetabounds),numel(selectedCells));

clear popzscore

plotrow = rc(1);
plotcolumn = rc(2);

bw= nan(length(selectedCells),1);
fr= nan(length(selectedCells),2);
pw=fr;
pwzs = nan(length(selectedCells),50);
pwtheta = nan(length(selectedCells),2);
figure(390);clf
numInterpPts = 16;

for d = 1:length(selectedCells)
    
    %ZSCORREEEEE
    currCell = selectedCells(d);
    
    pOn = round(U{currCell}.meta.poleOnset(1)*1000);
    restW(d) = nanmean(nanmean(U{currCell}.S_ctk(1,1:pOn,:)));
    
    counts = pop{currCell}.(fields{1}).counts;
    ranges = pop{currCell}.(fields{1}).range ;
    thetavect = [];
    for k = 1:length(ranges)
        thetavect = [thetavect ;ones(counts(k),1).*ranges(k)];
    end
    tvectNorm = normalize_var(thetavect,-1,1);
    
    allspikes = cell2mat(pop{currCell}.(fields{1}).raw);
    meanbase = mean(mean(allspikes(:,1:find(window==0)),2));
    stdbase = std(mean(allspikes(:,1:find(window==0)),2));
%     postresponses = allspikes(:,find(window==0)+5:find(window==0)+35);
    postresponses = allspikes(:,find(window==0)+ U{currCell}.meta.responseWindow(1) : find(window==0)+U{currCell}.meta.responseWindow(2)); 

    xresp = mean(postresponses,2);
    zscore = (xresp - meanbase) ./ stdbase;
    
    for k = 1:size(postresponses,1)
        spkT = find(postresponses(k,:)==1);
        if ~isempty(spkT) & numel(spkT)>1
            ISI = diff(spkT);
            CV(k) = std(ISI)./mean(ISI);
        else
            CV(k) = 0;
        end
    end
    
    [sorted,~,~] = binslin(thetavect,zscore,'equalE',numel(thetabounds),thetabounds(1),thetabounds(end));
    [CVsorted,~,~] = binslin(thetavect,CV','equalE',numel(thetabounds),thetabounds(1),thetabounds(end));

    
    samps = cellfun(@numel,sorted);
    selBins = find(samps>sum(samps)./100);
    
    for p=1:length(selBins)
        CVpop(p,:) = [thetabounds(selBins(p)) mean(CVsorted{selBins(p)})];
    end
    
    
    zraw = nan(length(sorted),2000);
    for b = 1:length(sorted)
        if sum(b == selBins)>0
            currvals = sorted{b};
            if ~isempty(sorted{b})
                zraw(b,1:length(currvals)) = currvals';
            end
        end
    end
    
    
    [p,~,stats]=anova1(zraw',[],'off');
    vals =multcompare(stats,[],'off');
    popzscore{d} = cellfun(@mean,sorted);
    otune(d)=p;
    

    x=cellfun(@str2num,stats.gnames);
    
    %CI BINS
    cibins = nan(size(x,1),1);
    
    SEMg = nanstd(zraw,[],2) ./ sqrt(sum(~isnan(zraw),2));
    for i = 1:length(x)
        rawx = zraw(x(i),:);
        SEM = SEMg(x(i));
        ts = tinv(.95,sum(~isnan(rawx),2)-1);      % T-Score
        cibins(i,:) = ts.*SEM;   %confidence intervals
    end
    
    
    y = stats.means;
%     ytune(x,d) = y;  %FOR PHASE
    
    x=thetabounds(x)';
    
    
    zstd = nanstd(zraw,[],2);
    ystd = zstd(~isnan(zstd))';
    
    xy = [x y'];
    
    kxy{d} = xy;

   
    figure(390);subplot(plotrow,plotcolumn,d);
    filex_shadedErrorBar(xy(:,1),xy(:,2),cibins,'-k'); 
%     scatter(xy(:,1),xy(:,2),[],[.7 .7 .7],'filled')
    set(gca,'ytick',round(min(xy(:,2))-1):round(max(xy(:,2))+1),'ylim',[round(min(xy(:,2))-1) round(max(xy(:,2))+1)],...
        'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
    
    nxy = [normalize_var(x,0,1),y'];
    iy(:,d) = interp1(nxy(:,1),nxy(:,2),linspace(0,1,numInterpPts));
    
    modely(:,d) = interp1(normalize_var(x,0,1),y,linspace(0,1,numInterpPts));
    modelystd(:,d) = interp1(normalize_var(x,0,1),ystd,linspace(0,1,numInterpPts));
    
    axis square
    
    zs = nanmean(zraw,2);
    zs(isnan(zs))=[];
    
    sigs =  vals(vals(:,end)<.01,:);
    
    clear peaksig
    [~,maxidx]=max(zs);
    midxsel = find(sum(sigs(:,[1 2])==maxidx,2)==1);
    
    if sum(midxsel>0)
        peaksig = sigs(midxsel,1:2);
        
        [~,midx] = min(abs(diff(peaksig')));
        pw(d,:) = peaksig(midx,:);
        figure(390);
        hold on; scatter(xy(peaksig(midx,:),1),xy(peaksig(midx,:),2),[],'filled','g')
        
        raws=zs(pw(d,1):pw(d,2));
        pwzs(d,1:length(raws)) = raws;
        pwtheta(d,:) = x(pw(d,:));
    end
    
end    


ocellidx = find(~isnan(pw(:,1)));
ocells = selectedCells(ocellidx);


%% LOCATION TUNING DISTRIBUTION PLOTTING

ocellidx = find(~isnan(pw(:,1)));
propTuned = numel(ocellidx)./numel(selectedCells)
pwstats = [nanmean(diff(pw')) nanstd(diff(pw'))];
pwstats*2


shapeD = normalize_var(pwzs',0,1);
figure(488);clf
for g = 1:size(shapeD,2)
    currVals = shapeD(:,g);
    plotvals = sort(currVals(~isnan(currVals)),'descend');
    hold on;plot(0:length(plotvals)-1,plotvals,'color',[.7 .7 .7])
    
end
set(gca,'xtick',0:4:16,'xticklabel',0:10:40,'ytick',0:.5:1)
ylabel('normalized z scored responses')
xlabel('width of tuning')
hold on; plot([pwstats(1) pwstats(1)],[0 1],'k-','linewidth',2)
hold on; plot([pwstats(1)-pwstats(2) pwstats(1)-pwstats(2)],[0 1],'k-.','linewidth',2)
hold on; plot([pwstats(1)+pwstats(2) pwstats(1)+pwstats(2)],[0 1],'k-.','linewidth',2)
axis square

%heat of location tuned cells
heat = normalize_var(iy,0,1);
heat = heat(:,ocellidx);
[~,hidx] = max(heat);
[~,sidx] = sort(hidx);

fheat = heat(:,sidx);
figure(990);clf;imagesc(flipud(fheat'))
set(gca,'ytick',[],'xtick',[])
caxis([0 1])

% range sampled for each cell
rangeExplored=sum(~isnan(cell2mat(popzscore)))*2;
reVals = [mean(rangeExplored) std(rangeExplored)];



%RAW ANGLE TUNING WITH Z SCORE VECTORS
maxt = cell2mat(cellfun(@max,kxy,'uniformoutput',0)');
mint = cell2mat(cellfun(@min,kxy,'uniformoutput',0)');

minmaxT = [[270-min(mint(:,1)) 270-30 270-max(maxt(:,1))]'];
bounds = [cosd(minmaxT) sind(minmaxT)].*[10 10 10]';

ztheta = [pwtheta(:,2) max(pwzs,[],2)];
% ztheta = [pwtheta(:,2)-restW' max(pwzs,[],2)]; %for normalizing to
% resting whisker position;

angs=[cosd(270-ztheta(:,1)) sind(270-ztheta(:,1))];
angsxy = ztheta(:,2).*angs;
figure(88);clf;subplot(2,1,1)
scatter(angsxy(find(selectedCells<=24),1),angsxy(find(selectedCells<=24),2),[],[.6 .6 .6],'filled')
hold on;scatter(angsxy(find(selectedCells>24),1),angsxy(find(selectedCells>24),2),'k','filled')
set(gca,'xlim',[-6 3],'ylim',[-6 0])
axis square
hold on;plot([0 bounds(1)],[0 bounds(3)],'-.k')
hold on;plot([0 0],[0 -10],'-k')
hold on;plot([0 -10],[0 0],'-k')
hold on;plot([0 bounds(3)],[0 bounds(6)],'-.k')

figure(88);subplot(2,1,2);
histogram(ztheta(find(selectedCells<=24),1),[round(min(mint(:,1))):2.5:round(max(maxt(:,1)))],'facecolor',[.6 .6 .6],'facealpha',1)
hold on;histogram(ztheta(find(selectedCells>24),1),[round(min(mint(:,1))):2.5:round(max(maxt(:,1)))],'facecolor',[0 0 0],'facealpha',1)
set(gca,'xlim',[-10 55])
axis square
xlabel('angle at touch');ylabel('number of cells')


%modeled tuning curves
widths = diff(pw');
widths = widths(~isnan(widths));
rawsigmas = ((widths*2)./2.58)';
modelsigmas = normrnd(pwstats(1)*2.5,pwstats(2)*2,length(rawsigmas),1);

figure(8);clf
rawTuning = ztheta(~isnan(ztheta(:,1)),1);
angTuning = datasample(min(mint(:,1)):2:max(maxt(:,1)),length(rawsigmas));
modelRaw = normrnd(mean(rawTuning),std(rawTuning),length(rawsigmas),1);

modeled = nan(301,length(rawsigmas));

for b = 1:length(rawsigmas)
    simSpks = normrnd(1,rawsigmas(b),10000,1);
    simTunecurves = histc(simSpks,[-50:50])./10000;
     subplot(3,1,1)
    hold on;plot(-50+rawTuning(b):50+rawTuning(b),simTunecurves,'color',[.7 .7 .7])
    modeled(150+rawTuning(b)-50:150+50+rawTuning(b),b) = simTunecurves;
    title('raw')
    set(gca,'xlim',[min(mint(:,1)) max(maxt(:,1))])
    
    subplot(3,1,2)
    hold on;plot(-50+modelRaw(b):50+modelRaw(b),simTunecurves,'g')
    title('raw normal simulated')
    set(gca,'xlim',[min(mint(:,1)) max(maxt(:,1))])
    
    subplot(3,1,3)
    hold on;plot(-50+angTuning(b):50+angTuning(b),simTunecurves,'r')
    title('uniform simulated')
end
set(gca,'xlim',[min(mint(:,1)) max(maxt(:,1))])
suplabel('whisker angle at touch','x')
suplabel('proportion of spikes','y')

figure(8);subplot(3,1,1)
summedPop = nansum(modeled,2);
selPlot = summedPop(min(mint(:,1))+150:max(maxt(:,1))+150);
hold on;plot(min(mint(:,1)):max(maxt(:,1)),selPlot,'k','linewidth',2)
% set(gca,'xtick',1:25:301,'xticklabel',-150:25:150,'xlim',[min(mint(:,1))+150 max(maxt(:,1))+150])
% xlabel('angle at touch');ylabel('summed norm. population rate activity')

%% PLotting out rasters
for d = 1:length(ocells)
    figure(8);clf
    
    currCell = ocells(d);
    allSpks = squeeze(U{currCell}.R_ntk);
    [~,idx] = sort(U{currCell}.meta.motorPosition);
    allSpks = allSpks(:,idx);
    for k = 1:size(allSpks,2)
        st = find(allSpks(:,k)==1);
        if ~isempty(st)
        figure(8);hold on
        scatter(st,ones(length(st),1).*k,[],'.k')
        end
    end
    pause
end

   
%% LOCATION TUNING MODELING
clearvars -except pop U
touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
% selectedCells=1:length(U);
[rc] = numSubplots(length(selectedCells));

plotrow = rc(1);
plotcolumn = rc(2);
%

clear popzscore


thetabounds = [-100:2.5:100];
numInterpPts = 16;

bins = [-50:10:50];

for g = 1:length(bins)
    for d = 1:length(selectedCells)
        
        %ZSCORREEEEE
        currCell = selectedCells(d);
        
        pOn = round(U{currCell}.meta.poleOnset(1)*1000);
        restW(d) = nanmean(nanmean(U{currCell}.S_ctk(1,1:pOn,:)));
        
        counts = pop{currCell}.theta.counts;
        ranges = pop{currCell}.theta.range ;
        thetavect = [];
        for k = 1:length(ranges)
            thetavect = [thetavect ;ones(counts(k),1).*ranges(k)];
        end
        tvectNorm = normalize_var(thetavect,-1,1);
        
        allspikes = cell2mat(pop{currCell}.theta.raw);
        meanbase = mean(mean(allspikes(:,1:51),2));
        stdbase = std(mean(allspikes(:,1:51),2));
        
        if bins(g) <0 
       postresponses = allspikes(:,51+bins(g):51+bins(g)+10);
        else
        postresponses = allspikes(:,51:51+bins(g));
        end
        xresp = mean(postresponses,2);
        zscore = (xresp - meanbase) ./ stdbase;
        
        [sorted,~,~] = binslin(thetavect,zscore,'equalE',numel(thetabounds),thetabounds(1),thetabounds(end));
        
        samps = cellfun(@numel,sorted);
        selBins = find(samps>sum(samps)./100);
        
        
        zraw = nan(length(sorted),2000);
        for b = 1:length(sorted)
            if sum(b == selBins)>0
                currvals = sorted{b};
                if ~isempty(sorted{b})
                    zraw(b,1:length(currvals)) = currvals';
                end
            end
        end
        
        
            [p,~,stats]=anova1(zraw',[],'off');
            vals =multcompare(stats,[],'off');
        %     popzscore{d} = cellfun(@mean,sorted);
        %     otune(d)=p;
        
        x = str2num(cell2mat(stats.gnames));
        x=thetabounds(x)';
        y = stats.means;
        
        zstd = nanstd(zraw,[],2);
        ystd = zstd(~isnan(zstd))';
        
        xy = [x y'];
        
        %     kxy{d} = xy;
        %     figure(390);subplot(plotrow,plotcolumn,d);
        %     scatter(xy(:,1),xy(:,2),[],[.7 .7 .7],'filled')
        %     set(gca,'ytick',round(min(xy(:,2))-1):round(max(xy(:,2))+1),'ylim',[round(min(xy(:,2))-1) round(max(xy(:,2))+1)],...
        %         'xlim',[min(xy(:,1)) max(xy(:,1))],'xtick',-20:20:60)
        
        %     nxy = [normalize_var(x,0,1),y'];
        %     iy(:,d) = interp1(nxy(:,1),nxy(:,2),linspace(0,1,numInterpPts));
        %
        modelY{g}(:,d) = interp1(normalize_var(x,0,1),y,linspace(0,1,numInterpPts));
        modelYstd{g}(:,d) = interp1(normalize_var(x,0,1),ystd,linspace(0,1,numInterpPts));
        
        
    end
end
%%
% modely
% modelystd
for b = 1:length(modelY)
    modely = modelY{b};
    modelystd = modelYstd{b} ;
    resampNum = 500;
    resampX = nan(size(modely,1)*resampNum,size(modely,2));
    for i = 1:resampNum
        resampX(size(modely,1)*(i-1)+1:size(modely,1)*(i-1)+size(modely,1),:) =  normrnd(modely,modelystd);
    end
    DmatX = resampX;
    
    DmatYnorm = repmat([1:size(modely,1)]',resampNum,1);
    randshuff = randperm(length(DmatYnorm));
    DmatYshuff = DmatYnorm(randshuff);
    
    DmatY = {DmatYnorm,DmatYshuff};
    
%     for d = 1:length(DmatY)
    for d = 1
        rando = randperm(length(DmatX));
        tmpDmatX=DmatX(rando,:);
        tmpDmatY=DmatY{d}(rando,:);
        txp = [];
        
        [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY{d})),0);
        [pred,opt_thresh,prob]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
            ,tmpDmatY(end*.7:end,:),'Max');
        txp = [txp ; tmpDmatY(end*.7:end) pred];
        
        cmat = confusionmat(txp(:,1),txp(:,2));
        ncmat = cmat./sum(cmat);
        figure(9+d);clf;
        imagesc(ncmat);
        caxis([0 .40])
        axis square
        colorbar
        title(['chance = ' num2str(1/size(modely,1))])
        set(gca,'xtick',[],'ytick',[])
        xlabel('truth');ylabel('predicted')
        
        predDiff(:,d) = txp(:,2) - txp(:,1);
    end
    model{b}=histcounts(predDiff(:,1),-16.5:16.5,'normalization','probability');
    stats(b,:) = [model{b}(17) std(predDiff)];
%     shuff{b}=histcounts(predDiff(:,2),-16.5:16.5,'normalization','probability');
    figure(40);clf
    plot(-16:16,model{b},'k')
%     hold on; plot(-16:16,shuff{b},'color',[.7 .7 .7])
    
    set(gca,'xtick',-16:8:16,'ytick',0:.1:.5,'ylim',[0 .35])
%     decodeResolution = std(predDiff)* (mean(rangeExplored)./size(modely,1))
end
%% NAIVE VS TRAINED

for g = 1:length(U)
    if strcmp('NL5b',U{g}.meta.layer)
        exp(g) = 1;
    else
        exp(g) = 2;
    end
end

ncells = find(exp==1);
ecells = find(exp==2);

propntouch = numel(intersect(ncells,selectedCells))./numel(ncells);
propetouch = numel(intersect(ecells,selectedCells))./numel(ecells);
figure(888);clf
subplot(2,2,[1 2]);scatter(propntouch,propetouch,[],'c','filled')


otcells = selectedCells(ocellidx);

propnobj = numel(intersect(ncells,otcells))./numel(intersect(ncells,selectedCells) );
propeobj =numel(intersect(ecells,otcells))./numel(intersect(ecells,selectedCells) );
subplot(2,3,[1 2 3]);hold on; scatter(propnobj,propeobj,[],'b','filled')
hold on; plot([0 1],[0 1],'-.k')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.5:1,'ytick',0:.5:1)
axis square
legend('prop cells TOUCH','prop touch cells OBJECT tuned')
xlabel('naive');ylabel('expert')

hold on;subplot(2,3,4);
bar(1:2,[propntouch propetouch],'c')
axis square
set(gca,'ylim',[0 1],'ytick',0:.5:1,'xtick',[])
subplot(2,3,5);bar(1:2,[propnobj propeobj],'k')
set(gca,'ylim',[0 1],'ytick',0:.5:1,'xtick',[])
axis square



[~,idx] = intersect(ncells,selectedCells)

sharpness = diff(pw');


nsharp = nanmean(sharpness(selectedCells<=24));
nstd = nanstd(sharpness(selectedCells<=24));
esharp = nanmean(sharpness(selectedCells>24));
estd = nanstd(sharpness(selectedCells>24));

subplot(2,3,6);errorbar(1:2,[nsharp esharp]*2.5,[nstd estd]*2.5,'ok')
hold on; scatter(ones(length(sharpness(selectedCells<=24)),1),sharpness(selectedCells<=24)*2.5,[],'k','filled')
hold on; scatter(ones(length(sharpness(selectedCells>24)),1)*2,sharpness(selectedCells>24)*2.5,[],'k','filled')
set(gca,'xlim',[0 3],'xtick',[])
axis square


%% PLOT ALL
[rc] = numSubplots(length(pop));

plotrow = rc(1);
plotcolumn = rc(2);

for d = 1:length(pop)
    
    
    figure(30);subplot(plotrow,plotcolumn,d);
    %       imagesc(imgaussfilt(pop{d}.(fields{1}).spikes,gaussFilt,'padding','replicate'));
    imagesc(pop{d}.(fields{1}).spikes)
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',(1:length(pop{d}.(fields{1}).range)),'yticklabel',[pop{d}.(fields{1}).range],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{d}.(fields{1}).counts,1)
        text(20,k,num2str(pop{d}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{d}.(fields{1}).range) 0],'w:')
    
    
    figure(31);subplot(plotrow,plotcolumn,d);
    imagesc(imgaussfilt(pop{d}.(fields{1}).spikes,gaussFilt,'padding','replicate'));
    %      imagesc(pop{d}.(fields{1}).spikes)
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',(1:length(pop{d}.(fields{1}).range)),'yticklabel',[pop{d}.(fields{1}).range],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{d}.(fields{1}).counts,1)
        text(20,k,num2str(pop{d}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{d}.(fields{1}).range) 0],'w:')
    %       axis('square')
    %             xlabel('time from touch onset (ms)')
    %             title(['theta'])
    
end
