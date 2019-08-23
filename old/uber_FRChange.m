%L3
%% Touch v no touch FR
L3TotalPre=zeros(1,length(L3));
L3TotalPost=zeros(1,length(L3));
L3TotalPreOff=zeros(1,length(L3));
L3TotalPostOff=zeros(1,length(L3));

for j = 1:length(L3)
    %ID Touch Onset FR
    touchIdx = [find(L3{j}.S_ctk(9,:,:)==1);find(L3{j}.S_ctk(12,:,:)==1)]; %all periods of touch
    touchIdx = touchIdx;
    spikesAligned = zeros(numel(touchIdx),201); %empty matrix for size
    spikes = squeeze(L3{j}.R_ntk);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-100:100]);
    end
    
    PreTouchSpikes = spikesAligned(:,1:100); %first 100ms BEFORE touch
    AvgPreTouchSpikes = mean(mean(PreTouchSpikes)); %population average pretouch
    L3TotalPre(j)=AvgPreTouchSpikes;    
    
    PostTouchSpikes=spikesAligned(:,101:150); %first 50ms AFTER touch
    PostTouchSpikesMean=mean(PostTouchSpikes);
    touchBin=zeros(1,size(PostTouchSpikes,2)/2);
    touchBin=(PostTouchSpikesMean(:,1:2:end-1)+PostTouchSpikesMean(:,2:2:end))/2;
    MaxPostTouchSpikes = (max(touchBin));
    L3TotalPost(j)=MaxPostTouchSpikes;
    
    %figure;
    %bar(-100:100,sum(spikesAligned)/numel(touchIdx),'k')    
    
    %ID Touch Offset FR  
    offtouchIdx = [find(L3{j}.S_ctk(10,:,:)==1);find(L3{j}.S_ctk(13,:,:)==1)]; %all periods of touch
    offtouchIdx = offtouchIdx;
    spikesAlignedoff = zeros(numel(offtouchIdx),201); %empty matrix for size
    spikes = squeeze(L3{j}.R_ntk);
    
    for i = 1:size(spikesAlignedoff,1)
        spikesAlignedoff(i,:) = spikes(offtouchIdx(i)+[-100:100]);
    end
       
    PreTouchOffSpikes = spikesAlignedoff(:,1:100); %first 100ms BEFORE touch offset
    AvgPreTouchOffSpikes = mean(mean(PreTouchOffSpikes)); %population average pretouch
    L3TotalPreOff(j)=AvgPreTouchOffSpikes;    
    
    PostTouchOffSpikes=spikesAlignedoff(:,101:150); %first 50ms AFTER touch offset
    PostTouchOffSpikesMean=mean(PostTouchOffSpikes);
    touchBin=zeros(1,size(PostTouchOffSpikes,2)/2);
    touchBin=(PostTouchOffSpikesMean(:,1:2:end-1)+PostTouchOffSpikesMean(:,2:2:end))/2;
    MaxPostTouchOffSpikes = (max(touchBin));
    L3TotalPostOff(j)=MaxPostTouchOffSpikes;
    
end

    L3TotalPost=L3TotalPost*1000; %convert all to spks/s
    L3TotalPre=L3TotalPre*1000;
    L3TotalPreOff=L3TotalPreOff*1000;
    L3TotalPostOff=L3TotalPostOff*1000;
%% Whisk v no Whisk FR (excluding periods of touch)
L3TotalHighWhisk=zeros(1,length(L3));
L3TotalLowWhisk=zeros(1,length(L3));

for j=1:length(L3)
        
        touchOnIdx = [find(L3{j}.S_ctk(9,:,:)==1); find(L3{j}.S_ctk(12,:,:)==1)];
        %touchOnIdx = touchOnIdx+70 %do I need to add this to mask out spikes for
        %periods post touch?
        touchOffIdx = [find(L3{j}.S_ctk(10,:,:)==1); find(L3{j}.S_ctk(13,:,:)==1)];
        touchOffIdx = touchOffIdx+70;

        touchEx_mask = ones(size(squeeze(L3{j}.S_ctk(1,:,:))));

        for i = 1:length(touchOnIdx)
            touchEx_mask(touchOnIdx(i):touchOffIdx(i)) = NaN;
        end
        touchEx_mask(1:100,:) = 1;
        
    highamp_mask=squeeze(L3{j}.S_ctk(3,:,:)>2.5); %ID amplitudes of high whisking 
    selectedHighamp = touchEx_mask.*highamp_mask; %mask out touches for periods of contact during high whisk
%     highamp_idx=find(~isnan(selectedHighamp));
    highamp_idx=find(~isnan(selectedHighamp));
    highamp_idxone = find(selectedHighamp(highamp_idx) ~= 0);
    highampSpikes = selectedHighamp(highamp_idxone);
    highampSpikes(:,2) = L3{j}.R_ntk(highamp_idxone);

    lowamp_mask=squeeze(L3{j}.S_ctk(3,:,:)<1.5); %ID amplitudes of no whisking 
    selectedLowamp = touchEx_mask.*lowamp_mask; %mask out touches
%     lowamp_idx=find(~isnan(selectedLowamp));
    lowamp_idx=find(~isnan(selectedLowamp));
    lowamp_idxone = find(selectedLowamp(lowamp_idx) ~= 0);
    lowampSpikes = selectedLowamp(lowamp_idxone);
    lowampSpikes(:,2) = L3{j}.R_ntk(lowamp_idxone);

    L3TotalHighWhisk(j)=mean(highampSpikes(:,2));
    L3TotalLowWhisk(j)=mean(lowampSpikes(:,2));

end
   L3TotalHighWhisk=L3TotalHighWhisk*1000;%conversion to spks/s
   L3TotalLowWhisk=L3TotalLowWhisk*1000;
%%    
%%%%%%%%%%%%%%%%%Repeat above but w/ L5b%%%%%%%%%%%%%%%%%%%%%%%%%
%% Touch v no touch FR
L5bTotalPre=zeros(1,length(L5b));
L5bTotalPost=zeros(1,length(L5b));
L5bTotalPreOff=zeros(1,length(L5b));
L5bTotalPostOff=zeros(1,length(L5b));

for j = 1:length(L5b)
    %ID Touch FR
    touchIdx = [find(L5b{j}.S_ctk(9,:,:)==1);find(L5b{j}.S_ctk(12,:,:)==1)]; %all periods of touch
    touchIdx = touchIdx;
    spikesAligned = zeros(numel(touchIdx),201); %empty matrix for size
    spikes = squeeze(L5b{j}.R_ntk);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-100:100]);
    end
    
    %figure;
    %bar(-100:100,sum(spikesAligned)/numel(touchIdx),'k')
    
    PreTouchSpikes = spikesAligned(:,1:100); %first 100ms BEFORE touch
    AvgPreTouchSpikes = mean(mean(PreTouchSpikes));
    L5bTotalPre(j)=AvgPreTouchSpikes;
    
    
    PostTouchSpikes=spikesAligned(:,101:150); %first 50ms AFTER touch
    PostTouchSpikesMean=mean(PostTouchSpikes); %average all each time point across all touch instances
    touchBin=zeros(1,size(PostTouchSpikes,2)/2);
    touchBin=(PostTouchSpikesMean(:,1:2:end-1)+PostTouchSpikesMean(:,2:2:end))/2;%bin average every 2 ms
    MaxPostTouchSpikes = (max(touchBin)); %take max value 
    L5bTotalPost(j)=MaxPostTouchSpikes;
    
    %ID Touch Offset FR  
    offtouchIdx = [find(L5b{j}.S_ctk(10,:,:)==1);find(L5b{j}.S_ctk(13,:,:)==1)]; %all periods of touch
    offtouchIdx = offtouchIdx(offtouchIdx<(L5b{j}.t*L5b{j}.k-101)); %makes sure that offset at very end of trials are eliminated to prevent error out
    spikesAlignedoff = zeros(numel(offtouchIdx),201); %empty matrix for size
    spikes = squeeze(L5b{j}.R_ntk);
    
    for i = 1:size(spikesAlignedoff,1)
        spikesAlignedoff(i,:) = spikes(offtouchIdx(i)+[-100:100]);
    end
       
    PreTouchOffSpikes = spikesAlignedoff(:,1:100); %first 100ms BEFORE touch
    AvgPreTouchOffSpikes = mean(mean(PreTouchOffSpikes));
    L5bTotalPreOff(j)=AvgPreTouchOffSpikes;    
    
    PostTouchOffSpikes=spikesAlignedoff(:,101:150); %first 50ms AFTER touch
    PostTouchOffSpikesMean=mean(PostTouchOffSpikes); %average all each time point across all touch instances
    touchBin=zeros(1,size(PostTouchOffSpikes,2)/2);
    touchBin=(PostTouchOffSpikesMean(:,1:2:end-1)+PostTouchOffSpikesMean(:,2:2:end))/2;%average across every 2 ms
    MaxPostTouchOffSpikes=(max(touchBin)); %take max value 
    L5bTotalPostOff(j)=MaxPostTouchOffSpikes;    
end
    L5bTotalPost=L5bTotalPost*1000; %convert all to spks/s
    L5bTotalPre=L5bTotalPre*1000;
    L5bTotalPreOff=L5bTotalPreOff*1000;
    L5bTotalPostOff=L5bTotalPostOff*1000;
%% Whisk v no Whisk FR (excluding periods of touch)
L5bTotalHighWhisk=zeros(1,length(L5b));
L5bTotalLowWhisk=zeros(1,length(L5b));

for j=1:length(L5b)
        
        touchOnIdx = [find(L5b{j}.S_ctk(9,:,:)==1); find(L5b{j}.S_ctk(12,:,:)==1)];
        %touchOnIdx = touchOnIdx+70 %do I need to add this to mask out spikes for
        %periods post touch?
        touchOffIdx = [find(L5b{j}.S_ctk(10,:,:)==1); find(L5b{j}.S_ctk(13,:,:)==1)];
        if touchOffIdx<(L5b{j}.t*L5b{j}.k-71)
           touchOffIdx = touchOffIdx+70;
        end
        touchEx_mask = ones(size(squeeze(L5b{j}.S_ctk(1,:,:))));

        for i = 1:length(touchOnIdx)
            touchEx_mask(touchOnIdx(i):touchOffIdx(i)) = NaN;
        end
        touchEx_mask(1:100,:) = 1;
        
    highamp_mask=squeeze(L5b{j}.S_ctk(3,:,:)>2.5); %ID amplitudes of high whisking 
    selectedHighamp = touchEx_mask.*highamp_mask; %mask out touches for periods of contact during high whisk
%     highamp_idx=find(~isnan(selectedHighamp));
    highamp_idx=find(~isnan(selectedHighamp));
    highamp_idxone = find(selectedHighamp(highamp_idx) ~= 0);
    highampSpikes = selectedHighamp(highamp_idxone);
    highampSpikes(:,2) = L5b{j}.R_ntk(highamp_idxone);

    lowamp_mask=squeeze(L5b{j}.S_ctk(3,:,:)<1.5); %ID amplitudes of no whisking 
    selectedLowamp = touchEx_mask.*lowamp_mask; %mask out touches
%     lowamp_idx=find(~isnan(selectedLowamp));
    lowamp_idx=find(~isnan(selectedLowamp));
    lowamp_idxone = find(selectedLowamp(lowamp_idx) ~= 0);
    lowampSpikes = selectedLowamp(lowamp_idxone);
    lowampSpikes(:,2) = L5b{j}.R_ntk(lowamp_idxone);

    L5bTotalHighWhisk(j)=mean(highampSpikes(:,2));
    L5bTotalLowWhisk(j)=mean(lowampSpikes(:,2));

end
   L5bTotalHighWhisk=L5bTotalHighWhisk*1000; %conversion to spks/s
   L5bTotalLowWhisk=L5bTotalLowWhisk*1000;
%% Analysis
% Population Geometric Mean =
    % log of each value, mean of log values, base 10^x=mean of log value. 
    % x=geometric mean
   L3PreOnMean=exp(mean(log(L3TotalPre)));
   L3PostOnMean=exp(mean(log(L3TotalPost)));
   L3PreOffMean=exp(mean(log(L3TotalPreOff)));
   L3PostOffMean=exp(mean(log(L3TotalPostOff)));
   L3HighWhiskMean=exp(mean(log(L3TotalHighWhisk)));
   L3LowWhiskMean=exp(mean(log(L3TotalLowWhisk)));
   
   L5bPreOnMean=exp(mean(log(L5bTotalPre)));
   L5bPostOnMean=exp(mean(log(L5bTotalPost)));
   L5bPreOffMean=exp(mean(log(L5bTotalPreOff)));
   L5bPostOffMean=exp(mean(log(L5bTotalPostOff)));
   L5bHighWhiskMean=exp(mean(log(L5bTotalHighWhisk)));
   L5bLowWhiskMean=exp(mean(log(L5bTotalLowWhisk)));
% Plot pre v post touch onset    
    figure;
    scatter(L3TotalPre,L3TotalPost,'y');
    hold on
    scatter(L5bTotalPre,L5bTotalPost,'r')
    hold on
    scatter(L3PreOnMean,L3PostOnMean,50,'y','x') %50 to correspond to size of marker
    hold on
    scatter(L5bPreOnMean,L5bPostOnMean,50,'r','x')
    set(gca,'XScale','log')
    set(gca,'XTickLabel',num2str(get(gca,'XTick').')) %scientific notation to integers
    set(gca,'YScale','log')
    set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
    hline=refline(1,0);
    axis square
    set(hline,'LineStyle',':') %dashed refline
    xlabel('Pre Touch (spks/s)')
    ylabel('Peak Touch (spks/s)')
    title('Change in Firing Rate from Touch Onset')   

% Plot pre v post touch onset    
    figure;
    scatter(L3TotalPreOff,L3TotalPostOff,'y');
    hold on
    scatter(L5bTotalPreOff,L5bTotalPostOff,'r')
    hold on
    scatter(L3PreOffMean,L3PostOffMean,50,'y','x')
    hold on
    scatter(L5bPreOffMean,L5bPostOffMean,50,'r','x')
    set(gca,'XScale','log')
    set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
    set(gca,'YScale','log')
    set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
    hline=refline(1,0);
    axis square
    set(hline,'LineStyle',':')
    xlabel('Pre Touch Offset (spks/s)')
    ylabel('Peak Touch Offset (spks/s)')
    title('Change in Firing Rate from Touch Offset')       
% Plot peak onset v peak offset
    figure;
    scatter(L3TotalPost,L3TotalPostOff,'c');
    hold on
    scatter(L5bTotalPost,L5bTotalPostOff,'r')
    hold on
    scatter(L3PostOnMean,L3PostOffMean,50,'c','x')
    hold on
    scatter(L5bPostOnMean,L5bPostOffMean,50,'r','x')
    set(gca,'XScale','log')
    set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
    set(gca,'YScale','log')
    set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
    hline=refline(1,0);
    axis square
    set(hline,'LineStyle',':')
    xlabel('Onset Peak(spks/s)')
    ylabel('Offset Peak(spks/s)')
    title('Peak Onset v Peak Offset')   
% Plot low whisk v high whisk    
    figure;
    scatter(L3TotalLowWhisk,L3TotalHighWhisk,'c')
    hold on
    scatter(L5bTotalLowWhisk,L5bTotalHighWhisk,'r');
    hold on
    scatter(L3LowWhiskMean,L3HighWhiskMean,50,'c','x')
    hold on
    scatter(L5bLowWhiskMean,L5bHighWhiskMean,50,'r','x')
    set(gca,'XScale','log')
    set(gca,'XTickLabel',num2str(get(gca,'XTick').'))
    set(gca,'YScale','log')
    set(gca,'YTickLabel',num2str(get(gca,'YTick').'))
    hline=refline(1,0);
    axis square
    set(hline,'LineStyle',':')
    xlabel('No Whisk (spks/s)')
    ylabel('Whisk (spks/s)')   
    title('Change in Firing Rate from Whisk')    
    
%Wilcoxon signed-rank test
[L3touchp L3touchres]=signrank(L3TotalPre,L3TotalPost)
[L5btouchp L5btouchres]=signrank(L5bTotalPre,L5bTotalPost)
[L3Offtouchp L3Offtouchres]=signrank(L3TotalPreOff,L3TotalPostOff)
[L5bOfftouchp L5bOfftouchres]=signrank(L5bTotalPreOff,L5bTotalPostOff)
[L3Whiskp L3whiskres]=signrank(L3TotalLowWhisk,L3TotalHighWhisk)
[L5bWhiskp L5bwhiskres]=signrank(L5bTotalLowWhisk,L5bTotalHighWhisk)
%touch res 0 = do NOT reject null; res 1 = reject null
%p value = probablity of obtaining result equal to or more extreme than
%what was observed.