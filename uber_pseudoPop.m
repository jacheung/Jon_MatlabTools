
numberSamples = 500;
boundaries = pop{1}.bounds{1};
binRaw = cell(1,length(boundaries));
pseudoPopRaw = cell(1,length(boundaries));
numcellsUsed = zeros(1,length(boundaries));

for k = 1:length(boundaries)
    
    currentBin = boundaries(k) ;
    
    for d = 1:length(pop)    
        if ~isempty(intersect(pop{d}.theta.range,currentBin))
            [~,idx] = intersect(pop{d}.theta.range,currentBin);
            binRaw{k}{d} = pop{d}.theta.raw{idx};
        end
    end
    
    resampled =[];
    
    if ~isempty(binRaw{k})
        filledIdx = find(~cellfun(@isempty,binRaw{k}));
        resampled = zeros(numberSamples,101);
        for j = 1:length(filledIdx)
            currentData = binRaw{k}{filledIdx(j)};
            resampled = resampled+datasample(currentData,500);
        end
        resampled = resampled./length(filledIdx);
        numcellsUsed(k) = numel(filledIdx);
    end
    
    pseudoPopRaw{k} = resampled;
end

sampledBounds = boundaries(~cellfun(@isempty,pseudoPopRaw));
numcellsUsed(cellfun(@isempty,pseudoPopRaw))=[];
pseudoPopRaw(cellfun(@isempty,pseudoPopRaw))=[];

sampledBounds(numcellsUsed<4) = [];
pseudoPopRaw(numcellsUsed<4)=[];

%Plotting Heatmap 
gaussFilt = [1]; %sigma filter for imgaussfilt (can be single value for filtering across rows and columns or two values [x y] for filter of row/column
popHeatmap = cell2mat(cellfun(@mean,pseudoPopRaw,'uniformoutput',0)')*1000;
figure(800);clf;imagesc(popHeatmap)
% figure(800);clf;imagesc(imgaussfilt(popHeatmap,gaussFilt,'padding','replicate'))


hold on; plot([51 51],[0 size(popHeatmap,1)],'w-.')
set(gca,'ytick',1:2:size(popHeatmap,1),'yticklabel',sampledBounds(1:2:end),'xtick',[1 51 101],'xticklabel',[-50 0 50])
xlabel('Time from touch onset (ms)');ylabel('Theta at touch') 
axis square
colorbar

%Plotting PSTH
colors= parula(length(pseudoPopRaw));
figure(580);clf
for g = 1:length(pseudoPopRaw)-1
    spkrate = smooth(mean(pseudoPopRaw{g}),5)*1000; %*1000 so that we get spks/s
%     spkratestd = smooth(std(pseudoPopRaw{g}));
    figure(580);
    hold on;plot(1:length(spkrate),spkrate,'color',colors(g,:))
end
plot([51 51],[0 50],'-.k')
set(gca,'xlim',[0 101],'xtick',[0 51 101],'xticklabel',[-50 0 50],'ytick',[0:10:50],'ylim',[0 40])
xlabel('Time from touch onset (ms)')
ylabel('Spike rate (spks/s)')

%% Decoding
numIterations    = 100;
DmatX = [];
DmatY = [];
for k = 1:length(sampledBounds)
    DmatX = [DmatX ;pseudoPopRaw{k}(:,51:end)];
    DmatY = [DmatY ;ones(size(pseudoPopRaw{k},1),1)*k];
end

% Logistic Classifier 
guesses = [];
rateDmatX = mean(DmatX,2); %predicting using average rate post touch

for reit = 1:numIterations   
   rando = randperm(length(DmatX));
   
   tmpDmatX=DmatX(rando,:);tmpDmatY=DmatY(rando,:);
    
    [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY)),0);
    Bfeat.theta{reit}=thetas;
    [pred,opt_thresh(reit),prob]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
        ,tmpDmatY(end*.7:end,:),'Max');
    guesses = [guesses ; [pred tmpDmatY(end*.7:end)]];
end

test = guesses(:,1) == guesses(:,2);
correctvals = guesses(test==1,:);
for k = 1:length(unique(guesses))
    correctPreds =  find(correctvals(:,1)==k);
    if ~isempty(correctPreds)
        accPercent(k) = numel(correctPreds) ./ sum(guesses(:,2)==k) ;
    else
        accPercent(k) = 0;
    end
end
figure(5320);clf;bar(1:length(accPercent),accPercent)
set(gca,'ylim',[0 1],'ytick',[0:.25:1],'xticklabel',sampledBounds)
xlabel('Theta at touch');ylabel('Proportion accurately predicted') 

%Confusion Matrix
[CFsorted,sortedBy,~] = binslin(guesses(:,2),guesses(:,1),'equalE',14,.5,13.5);
finalCF=zeros(numel(unique(guesses)),numel(unique(guesses)));
for b = 1:length(CFsorted) 
   [CF,sortedBy,~]= binslin(CFsorted{b},CFsorted{b},'equalE',14,.5,13.5);
   finalCF(:,b) = cellfun(@numel,CF) ./ sum( cellfun(@numel,CF));
end
figure(100);imagesc(finalCF);
set(gca,'xtick',[1:2:numel(unique(guesses))],'xticklabel',sampledBounds(1:2:end),'ytick',[1:2:numel(unique(guesses))],'yticklabel',sampledBounds(1:2:end))
caxis([0 .5])
colorbar
axis square
xlabel('Theta actual');ylabel('Theta predicted') 

% Maybe trying a nonlinear classifier using rate code may improve
% prediction

%Clustering... need to think about this for a little bit. How do we recover
%the identity? 
opts = statset('Display','final');
[idx,~] = kmeans(DmatX,13,'Distance','cityblock','Replicates',5,'Options',opts);
guesses = [idx DmatY];

%% population z scoring
thetabounds = [-100:5:100]; 
for i = 1:length(pop)
    
    counts = pop{i}.theta.counts;
    ranges = pop{i}.theta.range ;
    thetavect = [];
    for d = 1:length(ranges)
    thetavect = [thetavect ;ones(counts(d),1).*ranges(d)];
    end
    
allspikes = cell2mat(pop{i}.theta.raw);
meanbase = mean(mean(allspikes(:,1:51),2));
stdbase = std(mean(allspikes(:,1:51),2));
postresponses = allspikes(:,52:end);
xresp = mean(postresponses,2);
zscore = (xresp - meanbase) ./ stdbase;

uniquethetas = unique(thetavect);

[sorted,sortedBy,~] = binslin(thetavect,zscore,'equalE',numel(thetabounds),thetabounds(1),thetabounds(end));
popzscore{i} = cellfun(@mean,sorted);
end

zs = cell2mat(popzscore)';
xticks = thetabounds(:,~nansum(zs)==0);
zs(:,nansum(zs)==0)=[];
[~,idx] = max(zs,[],2);

[~,peakidx] = sort(idx);


figure(480);clf;imagesc(zs(peakidx,:))
h=colorbar;
colormap(redbluecmap)
caxis([-2 2])
set(gca,'xtick',[1:2:size(zs,2)],'xticklabel',xticks(1:2:end),'ytick',[])
xlabel('Theta at touch') 
%% tuning quant using cdfs and firing rate (need to consider choosing optimal window to count spikes).

for d = 1:length(pop)
binnedrates=mean(pop{d}.theta.spikes(:,51:76),2)*1000;


    current=sort(binnedrates./sum(binnedrates));
    
    previousVal = 0;
    cdfvals=zeros(length(current),1);
    for k = 1:length(current)
        previousVal = previousVal+current(k);
        cdfvals(k)=previousVal;
    end  

cdfvals = [0; cdfvals];

    figure(4030);clf;
    plot(0:length(cdfvals)-1,cdfvals)
    set(gca,'xlim',[0 length(cdfvals)-1],'ylim',[0 1])
    hold on; plot([0 length(cdfvals)-1],[0 1],'-.k')
    title(['n=' num2str(d) ' ' num2str(sum(cdfvals))])
    pause
end 



