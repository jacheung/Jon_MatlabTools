touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
xaxlim = 100;
STAlookbackwndow = -30;
[rc] = numSubplots(length(selectedCells));

%%
for rec = 1:length(selectedCells)

    % Sorted duration with curvature values. Along with plotting all spikes within those regions
    array=U{selectedCells(rec)};
    
    masks =  assist_touchmasks(array);
    %find touchidx
    eventsON = [find(array.S_ctk(9,:,:)==1);find(array.S_ctk(12,:,:)==1)];
    eventsOFF = [find(array.S_ctk(10,:,:)==1);find(array.S_ctk(13,:,:)==1)];
    
    tnum = ceil(eventsON./array.t);
    motors = array.meta.motorPosition;
    selMotors = motors(tnum);
    %if you want to use first touches only
    %    eventsON = [find(array.S_ctk(9,:,:)==1)];
    %    eventsOFF = [find(array.S_ctk(10,:,:)==1)];
    
    %Hmmm do I really need to get rid of those long touches? Probably not
    %since plotting can just get rid of those
    touchDurations = eventsOFF-eventsON;
    %     elimOutliersThreshold = mean(touchDurations)+2*std(touchDurations);
    %     eventsON(touchDurations>elimOutliersThreshold)=[];
    %     eventsOFF(touchDurations>elimOutliersThreshold)=[];
    %     touchDurations(touchDurations>elimOutliersThreshold)=[];
    
    filldis = nan(length(touchDurations),round(max(touchDurations),-2));
    spikefill = nan(length(touchDurations),round(max(touchDurations),-2));
    %Next step is get indices of varying ranges of touch durations and find
    %Dkappa.
    stim = squeeze(array.S_ctk(6,:,:)); %deltaKappa
    spikes = squeeze(array.R_ntk(:,:,:)); %spikes
    
    for d = 1:length(eventsON)
        windowsIdx = eventsON(d) : eventsOFF(d);
        filldis(d,1:length(windowsIdx)) = stim(windowsIdx);
        spikefill(d,1:length(windowsIdx)) = spikes(windowsIdx);
    end
    
    [~,idx] = sort(touchDurations);
    spkcount = nansum(spikefill,2);
    sortedspikefill=spikefill(idx,:);
    
    % STA for spikes within touch windows and finding if there is a DKAPPA related to it
    %Find all spike indices within touch window
    spikeIdx = find(spikes.*isnan(masks.touch)==1);
    STAvals = nan(length(spikeIdx),abs(STAlookbackwndow)+1);
    for k = 1:length(spikeIdx)
        STAvals(k,:) = stim(spikeIdx(k)+STAlookbackwndow:spikeIdx(k));
    end
    
    td{rec} = [touchDurations(touchDurations<500) spkcount(touchDurations<500) selMotors(touchDurations<500)'];
    
    %% plotting stuff
%     figure(380)
%     subplot(rc(1),rc(2),rec)
%     imagesc(filldis(idx,:));
%     set(gca,'ytick',[],'xlim',[0 xaxlim])
%     hold on ;
%     %plot total spk count on opposite side by doing xlim-spikecount val
%     plot(xaxlim-spkcount(idx)*5,1:length(spkcount),'r')
%     
%     % populate plot with all dots... very computationally expensive.
%     %     for i = 1:size(spikefill(idx,:),1)
%     %         spkIdx = find(sortedspikefill(i,:)==1);
%     %         if ~isempty(spkIdx)
%     %             hold on;
%     %             scatter(spkIdx,ones(length(spkIdx),1)*i,'filled','k')
%     %         end
%     %     end
%     
%     figure(480);
%     subplot(rc(1),rc(2),rec)
%     boundedline(1:31,nanmean(STAvals),nanstd(STAvals))
%     set(gca,'xlim',[0 abs(STAlookbackwndow)],'xtick',[0:10:abs(STAlookbackwndow)],'xticklabel',STAlookbackwndow:10:0)
%     
    
    
    
end

figure(679);clf
figure(689);clf
for k = 1:length(td)
    test = fitlm(td{k}(:,1),td{k}(:,2));
    corr(k) =sqrt(test.Rsquared.Ordinary);
    
    test2 = fitlm(td{k}(:,3),td{k}(:,1));
    corr1(k) =sqrt(test.Rsquared.Ordinary);
    
    
    figure(679);
    subplot(rc(1),rc(2),k)
    if corr(k)>.5
        hold on; scatter(td{k}(:,1),td{k}(:,2),'b.')
    else
        hold on; scatter(td{k}(:,1),td{k}(:,2),'k.')
    end
    
    figure(689);
    subplot(rc(1),rc(2),k)
    if corr1(k)>.5
    hold on;scatter(td{k}(:,3),td{k}(:,1),'.b')
    else 
        hold on;scatter(td{k}(:,3),td{k}(:,1),'.k')
    end
    set(gca,'xdir','reverse','xlim',[min(td{k}(:,3)) max(td{k}(:,3))],'xtick',[] )
    
end

figure(689);
xlabel('pole pos');ylabel('touch duration (ms)')
figure(679);
xlabel('touch duration');ylabel('spks during duration')
