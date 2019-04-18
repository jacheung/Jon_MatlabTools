%% Using equation below to calculate information from confusion matrix:
% I = P(stim).*Q(si|sd) * log2( Q(si|sd)*Q(sd))
% groupdata = singleCellRateD;
groupdata = singleCellRTD;

accuracy = nan(length(groupdata),2);
% but whatever this is outputting, it's providing me a negative value..
% might need to fix CF. Understand if CF must sum to 1 if probabilities? 
for d = 1:length(groupdata)
currentCF = groupdata{d}.CF;


tmp=currentCF./sum(currentCF);
accuracy(d,:)=[mean(sum(tmp.*eye(length(currentCF)))) 1/length(currentCF)];

vect=sum(currentCF);
pstim = vect./sum(vect);
% pstim = 1/length(currentCF);

HS(d) = sum(((pstim.*log2(pstim))*-1))
QDvect = sum(currentCF,2);
QD = repmat(QDvect./sum(QDvect),1,length(currentCF));

currentCF = currentCF./sum(currentCF);
MItmp  = (currentCF.*repmat(pstim,length(currentCF),1)) .* log2(currentCF./QD);
MI(d)=nansum(nansum(MItmp ))
end

figure(80);clf;histogram(MI,'normalization','probability','binwidth',.1)
set(gca,'xtick',0:.5:1.5,'ytick',0:.25:5,'ylim',[0 .5],'xlim',[0 1.5])
xlabel('Mutual information') 
ylabel('Proportion of cells') 

figure(81);clf;h=histogram(MI./HS,'normalization','probability','binwidth',.005)
set(gca,'xtick',0:.1:1,'ytick',0:.25:5,'ylim',[0 .35],'xlim',[0 .35])
xlabel('Proportion of Shannon Entropy explained') 
ylabel('Proportion of cells') 


propExplained = MI./HS;
[sorted,sortedBy,edges] = binslin(propExplained',propExplained','equalE',65,0,.32);
counts = cellfun(@numel,sorted);
proportion = counts./numel(counts) ;

for b = 1:length(edges)-1

    newX(b) = mean([edges(b) edges(b+1)]);
end

figure(4800);clf;bar(newX,proportion,'k');
set(gca,'xtick',0:.1:1,'ytick',0:.25:5,'ylim',[0 .35],'xlim',[0 .35])
xlabel('Proportion of Shannon Entropy explained') 
ylabel('Proportion of cells') 


figure(380); clf;scatter(MI./HS,accuracy(:,1)*100,'filled','k')
hold on; plot([0 1],[1/8 1/8]*100,'-.k')
set(gca,'xlim',[0 .35],'ylim',[0 50],'xtick',[0:.1:3],'ytick',[0:25:100])
xlabel('Proportion of Shannon Entropy explained');ylabel('Model accuracy')

%% Population information pooling
meanMaxInformation = mean(HS);
iterations = 25;
numSamples=100;
popDeviations = zeros(numSamples+1,iterations);
for g = 1:iterations 

growth = 0;
popContribution = zeros(numSamples+1,1);
for b = 1:numSamples
    growth = growth+datasample(MI,1);
    popContribution(b+1) = growth; 
end

popDeviations(:,g) = popContribution; 
end


popMeans=mean(popDeviations,2);
popSTDs = std(popDeviations,[],2);
popSTDs(popMeans>meanMaxInformation)=popSTDs(find(popMeans>meanMaxInformation==1,1));
popMeans(popMeans>meanMaxInformation)=meanMaxInformation;
popDeviations(popDeviations>meanMaxInformation) = meanMaxInformation;



figure(280);clf
for k = 1:iterations
   hold on; plot(0:numSamples,popDeviations(:,k),'color',[.8 .8 .8])
end
plot(0:numSamples,popMeans,'k','linewidth',6);
hold on;plot([0 numSamples],[meanMaxInformation meanMaxInformation],'-.k')
set(gca,'ylim',[0 round(meanMaxInformation + 1)],'xlim',[0 30],'ytick',0:1:4,'xtick',0:10:50)
xlabel('Number of cells');ylabel('Information (bits)') 


prcl = currentCF.*pstim
pr = sum(prcl,2) %marginal probability 
hr = -1 * (sum(pr.*log(pr)))
hl = -1 * (sum(pstim.*log(pstim)))
hrcl = -1 * nansum(nansum(prcl.*log(prcl)))
irl = hr + hl - hrcl
