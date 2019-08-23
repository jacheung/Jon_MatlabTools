
%% Resampling dataset to build training examples
%touch onset 
tOnset = 51;

%sampling from touch onset
tps = [5:5:50];

%coding scheme can be either spikecount OR spiketime
codingProperty = 'spikecount'; 

numIterations = 5;


propExplained = nan(numel(tps),length(pop));
codingacc= nan(numel(tps),length(pop));

for j = 1:numel(tps)
    numresamples = 500;
    
    popDmatX = cell(1,length(pop));
    popDmatY = cell(1,length(pop));
    for k = 1:length(pop)
        resampledDataX = [];
        resampledDataY = [];
        for g = 1:length(pop{k}.theta.raw)
            currentData = pop{k}.theta.raw{g}(:,tOnset:tOnset+tps(j));
            resampledDataX = [resampledDataX ;datasample(currentData,500)];
            resampledDataY = [resampledDataY ; ones(numresamples,1)*g];
        end
        
        popDmatX{k} = resampledDataX;
        popDmatY{k} = resampledDataY;
    end
    
    %% Logistic Classifier for single neuron

    %to study rate only and not time
    switch codingProperty
        case 'spikecount'
            popDmatX = cellfun(@(x) mean(x,2),popDmatX,'uniformoutput',0);
        case 'spiketime'
            popDmatX = popDmatX;
    end
    
    for f = 1:length(popDmatX)
        DmatX = popDmatX{f};
        DmatY = popDmatY{f};
        guesses = [];
        sampledBounds = pop{f}.theta.range;
        for reit = 1:numIterations
            rando = randperm(size(DmatX,1));
            
            tmpDmatX=DmatX(rando,:);tmpDmatY=DmatY(rando,:);
            
            [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY)),0);
            Bfeat.theta{reit}=thetas;
            [pred,opt_thresh(reit),prob]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
                ,tmpDmatY(end*.7:end,:),'Max');
            
            guesses = [guesses ; [pred tmpDmatY(end*.7:end)]];
        end
        
        %Confusion Matrix
        [CFsorted,sortedBy,~] = binslin(guesses(:,2),guesses(:,1),'equalE',numel(unique(guesses))+1,min(unique(guesses))-.5,max(unique(guesses))+.5);
        finalCF=zeros(numel(unique(guesses)),numel(unique(guesses)));
        for b = 1:length(CFsorted)
            [CF,sortedBy,~]= binslin(CFsorted{b},CFsorted{b},'equalE',numel(unique(guesses))+1,min(unique(guesses))-.5,max(unique(guesses))+.5);
            finalCF(:,b) = cellfun(@numel,CF);
        end
        
        accPercent = sum(finalCF.*eye(length(finalCF)))';
        
        
        groupdata{f}.accuracy = accPercent;
        groupdata{f}.bounds = sampledBounds;
        groupdata{f}.CF = finalCF;
    end
    %Mutual information calculator
    for d = 1:length(groupdata)
        currentCF = groupdata{d}.CF;
        
        tmp=currentCF./sum(currentCF);
        accuracy(d,:)=[mean(sum(tmp.*eye(length(currentCF)))) 1/length(currentCF)];
        
        vect=sum(currentCF);
        pstim = vect./sum(vect);
        % pstim = 1/length(currentCF);
        
        HS(d) = sum(((pstim.*log2(pstim))*-1));
        QDvect = sum(currentCF,2);
        QD = repmat(QDvect./sum(QDvect),1,length(currentCF));
        
        currentCF = currentCF./sum(currentCF);
        MItmp  = (currentCF.*repmat(pstim,length(currentCF),1)) .* log2(currentCF./QD);
        MI(d)=nansum(nansum(MItmp ));
    end
    
    propExplained(j,:) = MI./HS;
    codingacc(j,:) = accuracy(:,1)';
end
parcolors = parula(size(propExplained,1));
for p=1:size(propExplained,1)
figure(2320);hold on;h=scatter(propExplained(p,:),codingacc(p,:),'filled');
h.CData = parcolors(p,:);
end

%%
figure(2380);clf
hold on; plot(0:10,[0 ;mean(propExplainedRTD,2)]*100,'k');
hold on; plot(0:10,[0 ;mean(propExplainedR,2)]*100,'-.k');
set(gca,'xtick',2:2:length(tps),'xticklabel',tps(2:2:end),'xlim',[0 10],'ytick',[0:5:10],'ylim',[0 7.5])
xlabel('Integration window (ms) from touch onset')
ylabel('Percent of total information explained') 
legend('Spike time','Spike count','location','northwest') 