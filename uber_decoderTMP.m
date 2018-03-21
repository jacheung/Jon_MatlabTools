
%% GROUPING BY BINS OF THETAS 
window = [8 30];
countThresh = 10;

figure(50);clf;
for rec = 1:length(U)
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    
    fields = [      {'theta'}       {'velocity'}      {'amplitude'}    {'setpoint'}          {'phase'}           {'kappa'}];
    V.bounds = [{[-20:10:60]} {[-9750:125:9750]} {[-99.5:5:99.5]} {[-99.5:5:99.5]} {linspace(pi*-1,pi,12)} {[-.95:.05:.95]}];
    
    
    for d = 1
        var = squeeze(U{rec}.S_ctk(d,:,:));
        varAtTouch = var(touchIdx);
        varAtTouch = varAtTouch(~isnan(varAtTouch));
        touchIdx = touchIdx(~isnan(varAtTouch));
        spikesAtTouch = sum(spikes(repmat(touchIdx,1,numel(window(1):window(2)))+repmat([window(1):window(2)],numel(varAtTouch),1)),2);
        
        [sorted, sortedBy ,binBounds]=binslin(varAtTouch,spikesAtTouch,'equalE',numel(V.bounds{d})+1,V.bounds{d}(1),V.bounds{d}(end));
        binrange = V.bounds{d};
        
        trims=[binrange' cell2mat(cellfun(@mean,sorted,'Uniformoutput',0)) cell2mat(cellfun(@std,sorted,'Uniformoutput',0))];
        counttmp=cell2mat(cellfun(@size,sorted,'uniformoutput',0));
        mintouches=find(counttmp(:,1)>countThresh);
        
        
        %Plotting bins of FR distributions
        x = [0:.01:ceil(max(trims(mintouches,2)))];
        norm = normpdf(x,trims(:,2),trims(:,3));
%         norm = poisspdf(repmat(x,size(trims(:,2),1),1),trims(:,2));
        figure(50);subplot(5,4,rec)
        bloop = parula(size(binrange,2));
        
        for k = 1:length(mintouches)
            hold on; plot(x,norm(mintouches(k),:),'color',bloop(mintouches(k),:))
        end
        set(gca,'xtick',linspace(0,ceil(max(trims(mintouches,2))),5))
        xlabel(['Spikes/Touch within ' num2str(window(1)) ':' num2str(window(2)) 'ms window'])
        ylabel('Probability')
        % legend(num2str(trims(:,1)))
    end
    
    
    gaussmeans{rec} = trims(:,2);
    gaussstds{rec} = trims(:,3);
    gcounts{rec} = counttmp(:,1); 
    
end

collatcounts=[gcounts{:}];
thresh = collatcounts<countThresh;
abovethresh = ones(size(collatcounts));
abovethresh(thresh) = NaN;
collatgauss=[gaussmeans{:}];
collatstds = [gaussstds{:}];


for g = 13
    means = collatgauss(g,:).*abovethresh(g,:);
    stds = collatstds(g,:).*abovethresh(g,:);
    x = [0:.01:ceil(max(means))];
    figure(560);clf
    normvals = normpdf(repmat(x,size(means,2),1),means',stds');

    selectedcells=find(~isnan(means));
    
    randselection = randi(numel(selectedcells));
    
    hold on; plot(x,normvals(randselection,:))
    
    title(['prediction for ' num2str(binrange(g)) ' degree theta'])
end
% % for g = 1:length(U)
%   


