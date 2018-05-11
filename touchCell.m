function touchORnaw = touchCell(U,stdVal)
%Input U array and std threshold above baseline to classify as touch cell
%
baselinewindows = [1:50];
touchQuantWindows = [55:80];
%% build touch onset curves
clear Q

for rec=1:length(U)
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    spikes = squeeze(U{rec}.R_ntk);
    touchIdx = touchIdx(touchIdx<(numel(spikes)-151)); %elim last touches
    spikesAligned = zeros(numel(touchIdx),201);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-50:150]);
    end
    
    touchONN{rec}=spikesAligned;
    Q.touchraw{rec}=spikesAligned;
end


%% identify window after touch onset from when spiking returns to baseline

cellmean=cellfun(@mean, touchONN,'uniformoutput',0);
baselineFR=zeros(length(cellmean),1);
stdFR=baselineFR;
APstart=stdFR; %GREEN time post touch that AP rate reaches max ap rate/x(set at 3)
b2b=stdFR; %RED time post peak touch that AP rate returns back to baseline+sd
HRpeak=stdFR;%BLUE: time from post peak touch for AP rate to return to peak/x(set to 3 right now)
touchORnaw=zeros(1,length(U));

for d = 1:length(baselineFR)
    
    baselineFRtmp=mean(cellmean{d}(baselinewindows));
    stdFR(d)=std(cellmean{d}(baselinewindows));
    baselineFR(d)=baselineFRtmp+stdVal*stdFR(d); %baselineFR+std of FR 1 std
    
    if mean(cellmean{d}(touchQuantWindows))>baselineFR(d)
        touchORnaw(d) = 1;
        
        %If touch cell by our metric then quantify the below:
        Rpeak=50+find(cellmean{d}(51:150)==max(cellmean{d}(51:150)),1); %max response peak  after touch onset
        APstart(d)=50+find(cellmean{d}(51:Rpeak)>(max(cellmean{d}(51:150))/2),1);
        HRpeak(d)=Rpeak+find(cellmean{d}(Rpeak:150)<(max(cellmean{d}(51:150))/2),1);
        b2b(d)=Rpeak+(find(cellmean{d}(Rpeak:end)<=baselineFR(d),1)); %time point from -50 of when FR returns back to baseline
    else
    end
    
end

%%
figure(39);clf;
plotrow=8;
plotcol=8;


for k = 1:length(touchONN)
    subplot(plotrow,plotcol,k);
    avgSpikes = sum(touchONN{k})/size(touchONN{k},1);
    if touchORnaw(k)==1
        bar(-50:150,avgSpikes,'b');
    else
        b=bar(-50:150,avgSpikes,'k');
        b.FaceColor = [.5 .5 .5];
    end
    
    hold on; plot([0 0],[0 max(sum(touchONN{k})/size(touchONN{k},1))*1.5],'-.k','linewidth',2)
%     hold on; plot([b2b(k)-51 b2b(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'r','linewidth',2)
%     hold on; plot([HRpeak(k)-51 HRpeak(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'linewidth',2)
%     hold on; plot([APstart(k)-51 APstart(k)-51] ,[0 max(sum(touchONN{k})/size(touchONN{k},1))],'g','linewidth',2)
    
    set(gca,'xlim',[-50 100],'xtick',0,'ytick',[]);
%     if k==1
%         xlabel('Time from all touch onset (ms)')
%         ylabel('spks / ms')
%     end
end