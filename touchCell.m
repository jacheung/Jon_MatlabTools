%Input U array and std threshold above baseline to classify as touch cell
%Finds whether 5:30ms post touch is greater than the mean 50ms pre touch.
%Future edit could involve just grabbing peak post touch within that window
%and using that to define 

function touchORnaw = touchCell(U,stdVal,exceedSTDthresh)

baselinewindows = [1:51];
touchQuantWindows = [52:101];
rowColumns = numSubplots(length(U));
slidingwin = 10;
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
abovebaselineFR=zeros(length(cellmean),1);
belowbaselineFR=zeros(length(cellmean),1);
stdFR=abovebaselineFR;
APstart=stdFR; %GREEN time post touch that AP rate reaches max ap rate/x(set at 3)
b2b=stdFR; %RED time post peak touch that AP rate returns back to baseline+sd
HRpeak=stdFR;%BLUE: time from post peak touch for AP rate to return to peak/x(set to 3 right now)
touchORnaw=zeros(1,length(U));

for d = 1:length(abovebaselineFR)
    
    baselineFRtmp=mean(cellmean{d}(baselinewindows));
    stdFR(d)=std(cellmean{d}(baselinewindows));
    abovebaselineFR(d)=baselineFRtmp+stdVal*stdFR(d); %baselineFR+std of FR 1 std
    belowbaselineFR(d)=baselineFRtmp-stdVal*stdFR(d); 
    if belowbaselineFR(d) < 0 
        belowbaselineFR(d) = 0;
    end
    
    %SLIDING WINDOW 
    postexceedbase=double(cellmean{d}(touchQuantWindows)>abovebaselineFR(d));
    postbelowbase=double(cellmean{d}(touchQuantWindows)<=belowbaselineFR(d));

    
    postvalsUP = nan(length(postexceedbase)-slidingwin+1,1);
    postvalsDOWN = nan(length(postexceedbase)-slidingwin+1,1);
    for k1 = 1:length(postexceedbase)-slidingwin+1
     datawin = k1:k1+slidingwin-1;
        postvalsUP(k1) = mean(postexceedbase(datawin));
        postvalsDOWN(k1) = mean(postbelowbase(datawin));
    end
    
    if sum(postvalsUP>=exceedSTDthresh)>=1 %TOUCH ON 
        touchORnaw(d) = 1;
    elseif sum(postvalsDOWN>=exceedSTDthresh)>=1 %TOUCH OFF
        touchORnaw(d) = -1;
    else 
        touchORnaw(d) = 0;
    end
    
    tmp=(postvalsUP>=exceedSTDthresh);
    ftval(d,:) = max([find(postexceedbase==1,1) 0]);
    
    
    
%     elepass = sum(double(cellmean{d}(touchQuantWindows)>baselineFR(d))==1) ./ numel(touchQuantWindows); % percent of tps post baseline > baseline
%     
%     if  elepass>.2
%         touchORnaw(d) = 1;
%         
%         %If touch cell by our metric then quantify the below:
%         Rpeak=50+find(cellmean{d}(51:150)==max(cellmean{d}(51:150)),1); %max response peak  after touch onset
%         APstart(d)=50+find(cellmean{d}(51:Rpeak)>(max(cellmean{d}(51:150))/2),1);
%         HRpeak(d)=Rpeak+find(cellmean{d}(Rpeak:150)<(max(cellmean{d}(51:150))/2),1);
%         b2b(d)=Rpeak+(find(cellmean{d}(Rpeak:end)<=baselineFR(d),1)); %time point from -50 of when FR returns back to baseline
%     else
%     end
    
end

%%
figure(39);clf;
plotrow=rowColumns(1);
plotcol=rowColumns(2);


for k = 1:length(touchONN)
    subplot(plotrow,plotcol,k);
    avgSpikes = sum(touchONN{k})/size(touchONN{k},1);
    if touchORnaw(k)==1
        bar(-50:150,avgSpikes,'b');
    elseif touchORnaw(k) == -1 
        bar(-50:150,avgSpikes,'r');
    else
        b=bar(-50:150,avgSpikes,'k');
        b.FaceColor = [.5 .5 .5];
    end
    
    hold on; plot([0 0],[0 max(sum(touchONN{k})/size(touchONN{k},1))*1.5],'-.k','linewidth',2)
    set(gca,'xlim',[-50 100],'xtick',0,'ytick',[]);
    
    
    
    
end

tc=find(touchORnaw==1);
heatSpks = nan(length(tc),76);
for g = 1:sum(touchORnaw==1)
    avgSpikes = sum(touchONN{tc(g)})/size(touchONN{tc(g)},1);
    heatSpks(g,:) = normalize_var(avgSpikes(25:100),0,1);
end
[~,idx]  = max(heatSpks');
[~,ftidx] = sort(idx);


ntc = find(~(touchORnaw==1));
for g = 1:sum(~(touchORnaw==1))
    avgSpikes = sum(touchONN{ntc(g)})/size(touchONN{ntc(g)},1);
    ntheatSpks(g,:) = normalize_var(avgSpikes(25:100),0,1);
end


figure(48);clf;
subplot(3,1,[1 2])
imagesc(heatSpks(ftidx,:))
hold on; plot([26 26],[0 length(tc)],'-.w','linewidth',3)
colorbar
colormap(jet)
caxis([0 1])
set(gca,'xtick',1:25:76,'xticklabel',-25:25:75,'ytick',[])

subplot(3,1,[3])
imagesc(ntheatSpks)
hold on; plot([26 26],[0 length(ntc)],'-.w','linewidth',3)
colorbar
colormap(jet)
caxis([0 1])
set(gca,'xtick',1:25:76,'xticklabel',-25:25:75,'ytick',[])







