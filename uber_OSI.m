%% build touch onset curves 

for rec=1:length(U)
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    touchIdx = touchIdx(touchIdx<(numel(spikes)-151)); %elim last touches
    spikesAligned = zeros(numel(touchIdx),201);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-50:150]);
    end
    
    touchONN{rec}=spikesAligned;
end


%% identify window after touch onset from when spiking returns to baseline

cellmean=cellfun(@mean, touchONN,'uniformoutput',0);
baselineFR=zeros(length(cellmean),1);
stdFR=baselineFR;
b2b=stdFR;

for d = 1:length(baselineFR)
    baselineFRtmp=mean(cellmean{d}(25:50));
    stdFR(d)=std(cellmean{d}(25:50));
    baselineFR(d)=baselineFRtmp+stdFR(d); %baselineFR+std of FR
    Rpeak=51+find(cellmean{d}(51:150)==max(cellmean{d}(51:150)),1); %max response peak  after touch onset
    b2b(d)=Rpeak+(find(cellmean{d}(Rpeak:end)<=baselineFR(d),1)); %time point from -50 of when FR returns back to baseline
end

%% OPTIONAL plot window to look at spikes window
figure(39);clf;
plotrow=5;
plotcol=7;

for k = 1:length(touchONN)
    subplot(plotrow,plotcol,k);
    bar(-50:150,sum(touchONN{k})/size(touchONN{k},1),'k');
    hold on; plot([b2b(k)-50 b2b(k)-50] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))])
    %hold on; text(75,max(sum(touchONN{k})/size(touchONN{k},1)),num2str(osi(k)));
    set(gca,'xlim',[-50 150]);
    if k==1
        xlabel('Time from all touch onset (ms)')
        ylabel('spks / ms')
    end
end

%% build PSTH around touch index based on time from touch to back to baseline
values = cell(1,length(U));
grp=values;

for rec = 1:length(U)
    tptouch=b2b(rec)-50;
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    
    touchnID=ceil(touchIdx/4000);  
    goT=find(U{rec}.meta.trialType==1);
    [tf, loc] = ismember(touchnID, goT);
    %touchIdx=touchIdx(tf); %only look at touches in go trials so can easily compare between cont and semi-cont
   
    for var = [1]
        theta = squeeze(U{rec}.S_ctk(var,:,:));
        thetaAtTouch = theta(touchIdx);
        thetaAtTouch = thetaAtTouch(~isnan(thetaAtTouch));
        touchIdx = touchIdx(~isnan(thetaAtTouch));
        spkswndow=spikes(repmat(touchIdx,1,tptouch)+repmat([1:tptouch],numel(thetaAtTouch),1));%build matrix of spikes post T depending on time to return back to baseline post touch
        ISI=cell(1,size(spkswndow,1));
        for i = 1:size(spkswndow,1)
            spks=find(spkswndow(i,:)==1);
            spksdif=spks-[0 spks(1:end-1)];
            ISI{i}=spksdif(2:end);
        end
        ISI=cellfun(@nanmean,ISI);
        spikesAtTouch = sum(spkswndow,2);
        baselinePreTouch = sum(spikes(repmat(touchIdx,1,25)+repmat([-25:-1],numel(touchIdx),1)),2);
        grp{rec}=[thetaAtTouch spikesAtTouch baselinePreTouch ISI'];
    end
end

%% Touch Quant - mod IDX and ISI 
% Builds a modulation index that measures change in AP rate from 25ms pre
% touch to b2b val post touch
%final product modIdx
% col 1 = theta 
% col 2 = avg spk count post touch (determined by b2b)
% col 3 = avg spk count pre touch 
% col 4 = avg ISI post touch
% col 5 = avg spk count 
% col 6 = modulation idx (AP rate post - AP rate pre)/(AP rate post + AP
% rate pre).
Q.varNames = {'theta', 'AP rate post touch', 'AP rate pre touch', 'ISI post touch', 'touch count', 'mod idx'};
Q.TW = b2b';
modIdx=cell(length(U),1);

for rec=1:length(U)
    [sorted, sortedBy, binBounds]=binslin(grp{rec}(:,1),grp{rec},'equalE',41,-50,50);%1 degree bins
    touchspks=cell2mat(cellfun(@(x) nanmean(x,1),sorted,'uniformoutput',0));
    touchspks(:,5)=cellfun(@numel, sorted);
    touchspks(:,1)=linspace(mean(binBounds(1:2)),mean(binBounds(end-1:end)),length(binBounds)-1);
    
    modIdxEQ = (touchspks(:,2)-touchspks(:,3))./(touchspks(:,2)+touchspks(:,3));%spks post touch - spks 25ms pre touch/ spks post touch + spks25mspretouch
    modIdx{rec}=[touchspks modIdxEQ];
    
    
    selected=find(modIdx{rec}(:,5)>nansum(modIdx{rec}(:,5))*.015); %selecting only for thetas that have at least 1.5% of total touches.
    figure(10);subplot(5,7,rec)
    plot(modIdx{rec}(selected,1),modIdx{rec}(selected,6))
    set(gca,'xlim',[min(modIdx{rec}(selected,1)) max(modIdx{rec}(selected,1))],'ylim',[-1 1]);
    
    figure(11);subplot(5,7,rec)
    plot(modIdx{rec}(selected,1),modIdx{rec}(selected,4));
    set(gca,'xlim',[min(modIdx{rec}(selected,1)) max(modIdx{rec}(selected,1))]);
    
end
Q.val = modIdx;



%% OSI: bin spike data (based on 1 degree bin at time of writing) and then do OSI
osi=zeros(length(U),1);
for rec=1:length(U)
    [sorted, sortedBy, binBounds]=binslinTMP(grp{rec}(:,1),grp{rec},'equalE',26,-50,50);%1 degree bins
    touchspks=cell2mat(cellfun(@(x) mean(x,1),sorted,'uniformoutput',0));
    touchspks(:,1)=floor(touchspks(:,1));
    figure(9);bar(touchspks(:,1),touchspks(:,2))
    [~,spkpref]=max(touchspks);
    Rpref=touchspks(spkpref(2),:);
    Rorth=(max(touchspks(:,1))-min(touchspks(:,1)))/2;
    osi(rec)=(Rpref(1)-Rorth)/Rpref(1);
end



