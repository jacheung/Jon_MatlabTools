%Doing normalized resample... resampling to 8 bins of touch
numresamples = 500;
maxBins = 6;
allBounds = nan(15,length(pop));
allCounts = nan(15,length(pop));

for k = 1:length(pop)
    allBounds(1:length(pop{k}.theta.range),k) = pop{k}.theta.range;
    allCounts(1:length(pop{k}.theta.counts),k) = pop{k}.theta.counts;
end

numBinsCurrently=sum(~isnan(allBounds));
shrinkThese = find(numBinsCurrently>maxBins);
countsSelectedIdx = nan(maxBins,length(pop));

for t=1:length(pop)
    currIdx = t;
    delTmp = allCounts(~isnan(allCounts(:,currIdx)),currIdx);
    finalComparison = delTmp;
    
    if length(delTmp) > maxBins
        stepsAway = length(delTmp)-maxBins;
        tossed = zeros(stepsAway,1);
        for g=1:stepsAway
            tossVals = [delTmp(1) delTmp(end)];
            delVals = delTmp == min(tossVals);
            if sum(delVals)>1
                middleVal = mean(1:length(delVals));
                tmp=(1:length(delVals)) - middleVal;
                [~,idx] = max(abs(tmp'.*delVals));
                updatedDelvals = zeros(length(delVals),1);
                updatedDelvals(idx)=1;
                delTmp(logical(updatedDelvals))=[];

            else
                updatedDelvals = delTmp == min(tossVals);
                delTmp(updatedDelvals)=[];

            end

        end
          [vals,idx]=intersect(finalComparison,delTmp);
          eventualIdx = [];
          if numel(idx)<maxBins
             for u= 1:length(vals)
                 eventualIdx = [eventualIdx ; find(finalComparison == vals(u))];
             end
             countsSelectedIdx(1:maxBins,t) = sort(eventualIdx);
          else
              countsSelectedIdx(1:maxBins,t) = sort(idx);
          end
          
     else
         countsSelectedIdx(1:maxBins,t) =1:maxBins;
     end
end

%% Using bins with most counts in them...
popDmatX = cell(1,length(pop));
popDmatY = cell(1,length(pop));
thetaInBins = nan(maxBins,length(pop));
for k = 1:length(pop)
    resampledDataX = [];
    resampledDataY = [];
    selectedBins = countsSelectedIdx(:,k);
    for g = 1:length(selectedBins)
        bin = selectedBins(g);
        currentData = pop{k}.theta.raw{bin}(:,51:end);
        resampledDataX = [resampledDataX ;datasample(currentData,500)];
        %         resampledDataY = [resampledDataY ; ones(numresamples,1)*pop{k}.theta.range(bin)];
        resampledDataY = [resampledDataY ; ones(numresamples,1).*g];
    end
    thetaInBins(1:maxBins,k)= pop{k}.theta.range(selectedBins);
    popDmatX{k} = resampledDataX;
    popDmatY{k} = resampledDataY;
end

%% PCA of all 50 bins

tmp=cell2mat(popDmatX); %each row is a single theta with ALL cells 
firstPC=zeros(51,length(tmp));
secondPC = zeros(51,length(tmp));
thirdPC = zeros(51,length(tmp));
for k = 1:length(tmp)
test=reshape(tmp(k,:),length(popDmatX),51);
[coeff,score,latent,tsquared,explained,mu] = pca(test');
firstPC(:,k) = score(:,[1]);
secondPC(:,k) = score(:,[2]);
thirdPC(:,k) = score(:,[3]);
explainedAll(k) = sum(explained(1:3));
end

colorsPar = parula(6);

figure(30480); clf
hold on; h = scatter3(mean(firstPC(:,1:500),2),mean(secondPC(:,1:500),2),mean(thirdPC(:,1:500),2));
h.CData = colorsPar(1,:);
hold on;h=scatter3(mean(firstPC(:,501:1000),2),mean(secondPC(:,501:1000),2),mean(thirdPC(:,501:1000),2))
h.CData = colorsPar(2,:);
hold on;h=scatter3(mean(firstPC(:,1001:1500),2),mean(secondPC(:,1001:1500),2),mean(thirdPC(:,1001:1500),2))
h.CData = colorsPar(3,:);
hold on;h=scatter3(mean(firstPC(:,1501:2000),2),mean(secondPC(:,1501:2000),2),mean(thirdPC(:,1501:2000),2))
h.CData = colorsPar(4,:);
hold on;h=scatter3(mean(firstPC(:,2001:2500),2),mean(secondPC(:,2001:2500),2),mean(thirdPC(:,2001:2500),2))
h.CData = colorsPar(5,:);
hold on;h=scatter3(mean(firstPC(:,2501:3000),2),mean(secondPC(:,2501:3000),2),mean(thirdPC(:,2501:3000),2))
h.CData = colorsPar(6,:);


