%%ELIM{rec} b/c only one cell
for rec = 1
    touchIdx = [find(d.S_ctk(9,:,:)==1);find(d.S_ctk(12,:,:)==1)];
    touchIdx = touchIdx
    spikesAligned = zeros(numel(touchIdx),201);
    spikes = squeeze(d.R_ntk);
 
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-50:150]);
    end
    figure(1);clf;
    bar(-50:150,sum(spikesAligned)/numel(touchIdx),'k')
    set(gca,'xlim',[-50 150])
    xlabel('Time from all touch onset (ms)')
    ylabel('Spks / ms')
   % print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\vpmlayerfour\figures\150417\PCTHall_' num2str(rec)])
end
% %% ORIGINAL 
% for rec = 1:11
%     touchIdx = [find(d{rec}.S_ctk(9,:,:)==1);find(d{rec}.S_ctk(12,:,:)==1)];
%  touchIdx = touchIdx
%     spikesAligned = zeros(numel(touchIdx),201);
%     spikes = squeeze(d{rec}.R_ntk);
%  
%     for i = 1:size(spikesAligned,1)
%         spikesAligned(i,:) = spikes(touchIdx(i)+[-50:150]);
%     end
%     figure(1);clf;
%     bar(-50:150,sum(spikesAligned)/numel(touchIdx),'k')
%     set(gca,'xlim',[-50 150])
%     xlabel('Time from all touch onset (ms)')
%     ylabel('Spks / ms')
%    % print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\vpmlayerfour\figures\150417\PCTHall_' num2str(rec)])
% end
%%
baseRate = zeros(21,1);
sampleRate = zeros(21,1);
peakRate = zeros(21,1);
spikes = {};
ft_spikesAligned={};
lt_spikesAligned = {};
L=U
for rec = [1]
    
    ftIdx = find(L{rec}.S_ctk(9,:,:)==1);
    ltIdx = find(L{rec}.S_ctk(12,:,:)==1);
    ltIdx = ltIdx(ltIdx < L{rec}.t*L{rec}.k-150);
    
    % clean up ftIdx to exclude second touch within 40ms
    
    diffMat = repmat(ftIdx,1,size(ltIdx,1))-repmat(ltIdx',size(ftIdx,1),1);
    
    doubleTouch = find(sum(diffMat>=-40 & diffMat <0,2));
    exd_ftIdx = ftIdx(setdiff(1:length(ftIdx),doubleTouch));
    
    
    ft_spikesAligned{rec} = zeros(numel(ftIdx),201);
    lt_spikesAligned{rec} = zeros(numel(ltIdx),201);
    exd_ft_spikesAligned{rec} = zeros(numel(exd_ftIdx),201);
    spikes{rec} = squeeze(L{rec}.R_ntk);
    
    
    % sampleSpikes = []
    % for i = 1:L{rec}.k
    %     sampleSpikes = cat(2,sampleSpikes, spikes(find(L{rec}.S_ctk(15,:,i))));
    %
    %
    % end
    
    for i = 1:size(ft_spikesAligned{rec},1)
        ft_spikesAligned{rec}(i,:) = spikes{rec}(ftIdx(i)+[-50:150]);
    end
    
    for i = 1:size(lt_spikesAligned{rec},1)
        
        lt_spikesAligned{rec}(i,:) = spikes{rec}(ltIdx(i)+[-50:150]);
    end
    
    for i = 1:size(exd_ft_spikesAligned{rec},1)
        exd_ft_spikesAligned{rec}(i,:) = spikes{rec}(exd_ftIdx(i)+[-50:150]);
    end
    
    %baseRate(rec) = sum(sum(spikes{rec}(200:450,:)))/numel(spikes{rec}(200:450,:))*1000
    %sampleRate(rec)  = sum(resL(rec).spkCountExploration)/sum(resL(rec).timeExploration)*1000;
    
    peakRate(rec) = max(smooth(sum(ft_spikesAligned{rec}),3)/length(ft_spikesAligned{rec}))*1000;
end
%% Outside of C2 (obsolete)

L4farCells = [1];
for rec = 1
    
    ftIdx_far = find(L{L4farCells(rec)}.S_ctk(9,:,:)==1)
    ltIdx_far = find(L{L4farCells(rec)}.S_ctk(12,:,:)==1)
    ltIdx_far = ltIdx_far(mod(ltIdx_far,L{L4farCells(rec)}.t)<(L{L4farCells(rec)}.t-150) & mod(ltIdx_far,L{L4farCells(rec)}.t)>800);
    ft_spikesAligned_far{rec} = zeros(numel(ftIdx_far),201);
    lt_spikesAligned_far{rec} = zeros(numel(ltIdx_far),201);
    spikes_far{rec} = squeeze(L{L4farCells(rec)}.R_ntk);
    
    
    % sampleSpikes = []
    % for i = 1:L4{rec}.k
    %     sampleSpikes = cat(2,sampleSpikes, spikes(find(L4{rec}.S_ctk(15,:,i))));
    %
    %
    % end
    
    for i = 1:size(ft_spikesAligned_far{rec},1)
        ft_spikesAligned_far{rec}(i,:) = spikes_far{rec}(ftIdx_far(i)+[-50:150]);
    end
    
    for i = 1:size(lt_spikesAligned_far{rec},1)
        lt_spikesAligned_far{rec}(i,:) = spikes_far{rec}(ltIdx_far(i)+[-50:150]);
    end
    
    
    baseRate_far(rec) = sum(sum(spikes_far{rec}(200:450,:)))/numel(spikes_far{rec}(200:450,:))*1000
    sampleRate_far(rec)  = sum(resL4CX(rec).spkCountExploration)/sum(resL4C2(rec).timeExploration)*1000;
    
    peakRate_far(rec) = max(smooth(sum(ft_spikesAligned_far{rec}),3)/length(ft_spikesAligned_far{rec}))*1000;
end

%plot(smooth(sum(ft_spikesAligned{rec}),5)/length(ft_spikesAligned{rec}),'r')


%max(smooth(sum(ft_spikesAligned{rec}),5)/length(ft_spikesAligned{rec}))*1000

%%
figure(1); clf;set(gcf,'paperposition',[0 0 2 2])
hold on

for i = 1:17
    plot([baseRate sampleRate peakRate]','ko-','markersize',4,'linewidth',1)
end
set(gca,'xlim',[0.8 3.2])
%ylabel('Spike Rate')
%set(gca,'xtick',[1 2 3],'xticklabel',{'Baseline','Exploration','Peak Touch'})
set(gca,'xtick',[1 2 3],'xticklabel',[])
%print(gcf,'-depsc', 'C:\Users\shires\Dropbox\NoiseProject\Manuscript\Fig2\Materials\BaseSamplePeakRate')


%Outtside of C2
figure(2); clf;set(gcf,'paperposition',[0 0 2 2])
hold on

for i = 1:11
    plot([baseRate_far; sampleRate_far; peakRate_far],'ko-','markersize',4,'linewidth',1)
end
set(gca,'xlim',[0.8 3.2])
%ylabel('Spike Rate')
%set(gca,'xtick',[1 2 3],'xticklabel',{'Baseline','Exploration','Peak Touch'})
set(gca,'xtick',[1 2 3],'xticklabel',[])
%print(gcf,'-depsc', 'C:\Users\shires\Dropbox\NoiseProject\Manuscript\Fig2\Materials\BaseSamplePeakRate_far')


% Overlay
figure(3); clf;set(gcf,'paperposition',[0 0 2 2])
hold on

for i = 1:21
    
    plot([baseRate sampleRate peakRate]','-','color',[240 177 177]/255,'markersize',4,'linewidth',1)
    plot([baseRate sampleRate peakRate]','o','color',[240 100 100]/255,'markersize',4,'linewidth',1)
end

for i = 1:11
    plot([baseRate_far; sampleRate_far; peakRate_far],'-','color',[.75 .75 .75],'markersize',4,'linewidth',1)
    
    plot([baseRate_far; sampleRate_far; peakRate_far],'o','color',[.5 .5 .5],'markersize',4,'linewidth',1)
end

plot(mean([baseRate_far; sampleRate_far; peakRate_far]'),'k+')
plot(mean([baseRate sampleRate peakRate]),'r+')
set(gca,'xlim',[0.8 3.2])
set(gca,'xtick',[1 2 3],'xticklabel',[])
%print(gcf,'-depsc', 'C:\Users\shires\Dropbox\NoiseProject\Manuscript\Fig2\Materials\BaseSamplePeakRate_overlay')

%%

for i = 1:21
    figure(1);clf;
    set(gcf,'paperposition',[0 0 2.5 1])
    bar([-50:150],1000*sum(ft_spikesAligned{i})/size(ft_spikesAligned{i},1),'facecolor',[.5 .5 1],'barwidth',1,'linestyle','none')
    hold on
    bar([-50:150],1000*sum(lt_spikesAligned{i})/size(lt_spikesAligned{i},1),'facecolor',[1 .5 .5],'barwidth',1,'linestyle','none')
    plot([-50:150],smooth(sum(ft_spikesAligned{i}),3)/size(ft_spikesAligned{i},1)*1000,'b','linewidth',.5)
    
    plot([-50:150],smooth(sum(lt_spikesAligned{i}),3)/size(lt_spikesAligned{i},1)*1000,'r','linewidth',.5)
    box off
    
    
    %xlabel('Time from touch')
    
    %ylabel('Spikes / s')
    %legend('First touch','Late touch')
    set(gca,'xlim',[-50 150],'xticklabel',[])
   % print(gcf,'-depsc', ['C:\Users\shires\Dropbox\NoiseProject\Manuscript\sFigures\SFig2\Materials\FirstLatePCTHsmall_Cell_' num2str(i)])
    
end

%% Fig 1 Normalized PCTH heatmap
figure(2);clf;set(gcf,'paperposition',[0 0 3 2])
norm_ft_spikesAligned = zeros(21,201);
cmap = jet(256);
colormap(cmap(1:224,:));
nearCells = [1:17 38:51];
for i = [1:32 38:51]
    norm_ft_spikesAligned(i,:) = sum(ft_spikesAligned{i})/size(ft_spikesAligned{i},1);
    norm_ft_spikesAligned(i,:) = norm_ft_spikesAligned(i,:)/max(norm_ft_spikesAligned(i,50:98));
    norm_lt_spikesAligned(i,:) = sum(lt_spikesAligned{i})/size(lt_spikesAligned{i},1);
    norm_lt_spikesAligned(i,:) = norm_lt_spikesAligned(i,:)/max(norm_lt_spikesAligned(i,:));
    norm_exd_ft_spikesAligned(i,:) = sum(exd_ft_spikesAligned{i})/size(exd_ft_spikesAligned{i},1);
    norm_exd_ft_spikesAligned(i,:) = norm_exd_ft_spikesAligned(i,:)/max(norm_exd_ft_spikesAligned(i,50:90));
end

% sort by weighted response time
weightedResponse=[];
for i = [1:32 38:51]
    alignby = norm_ft_spikesAligned(i,:)
    if ~isempty(find(alignby(56:75)>.75,1,'first'))
        tmp = find(alignby(56:75)>.75,1,'first');
    else [~, tmp] = max(alignby(50:75))
    end
    weightedResponse(i) = tmp;
end
farCells = [22:30 32];
[~,near_idx] = sort(weightedResponse([1:17 38:51]))
near_idx = nearCells(near_idx)
[~,far_idx] = sort(weightedResponse([farCells]))
far_idx = farCells(far_idx);

image(224*cat(1,norm_ft_spikesAligned(near_idx,25:100),zeros(1,76),norm_ft_spikesAligned(far_idx,25:100)))
title('Normalized first touch response')
set(gca,'xlim', [0 75],'xtick',[ 25 50 75 100],'xticklabel',[0 25 50 75])
ylabel('Cells')
xlabel('Time from touch (ms)')
box off

h_cb = colorbar
set(h_cb,'ytick',[0 .5 1])
% imagesc(cat(1,norm_lt_spikesAligned(1:17,:),zeros(1,201),norm_lt_spikesAligned_far))
% title('Normalized late touch response')
%
% set(gca,'xtick',[1 51 101 151],'xticklabel',[-50 0 50 100])
% ylabel('Cells')
% xlabel('Time from touch (ms)')
% box off
%print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\StackedPCTHHeatmapNearFar'])


%% L5C2 Normalized PCTH heatmap

baseRate = zeros(21,1);
sampleRate = zeros(21,1);
peakRate = zeros(21,1);
spikes = {};
ft_spikesAligned={};

for rec = [1:11]
    
    ftIdx = find(L5{rec}.S_ctk(9,:,:)==1);
    ltIdx = find(L5{rec}.S_ctk(12,:,:)==1);
    ltIdx = ltIdx(ltIdx < L5{rec}.t*L5{rec}.k-150);
    
    % clean up ftIdx to exclude second touch within 40ms
    
    diffMat = repmat(ftIdx,1,size(ltIdx,1))-repmat(ltIdx',size(ftIdx,1),1);
    
    doubleTouch = find(sum(diffMat>=-40 & diffMat <0,2));
    exd_ftIdx = ftIdx(setdiff(1:length(ftIdx),doubleTouch));
    
    
    ft_spikesAligned{rec} = zeros(numel(ftIdx),201);
    lt_spikesAligned{rec} = zeros(numel(ltIdx),201);
    exd_ft_spikesAligned{rec} = zeros(numel(exd_ftIdx),201);
    spikes{rec} = squeeze(L5{rec}.R_ntk);
    
    
    % sampleSpikes = []
    % for i = 1:L5{rec}.k
    %     sampleSpikes = cat(2,sampleSpikes, spikes(find(L5{rec}.S_ctk(15,:,i))));
    %
    %
    % end
    
    for i = 1:size(ft_spikesAligned{rec},1)
        ft_spikesAligned{rec}(i,:) = spikes{rec}(ftIdx(i)+[-50:150]);
    end
    
    for i = 1:size(lt_spikesAligned{rec},1)
        
        lt_spikesAligned{rec}(i,:) = spikes{rec}(ltIdx(i)+[-50:150]);
    end
    
    for i = 1:size(exd_ft_spikesAligned{rec},1)
        exd_ft_spikesAligned{rec}(i,:) = spikes{rec}(exd_ftIdx(i)+[-50:150]);
    end
    
    %baseRate(rec) = sum(sum(spikes{rec}(200:450,:)))/numel(spikes{rec}(200:450,:))*1000
    %sampleRate(rec)  = sum(resL5(rec).spkCountExploration)/sum(resL5(rec).timeExploration)*1000;
    
    peakRate(rec) = max(smooth(sum(ft_spikesAligned{rec}),3)/length(ft_spikesAligned{rec}))*1000;
end

figure(2);clf;set(gcf,'paperposition',[0 0 3 2])
norm_ft_spikesAligned = zeros(21,201);
norm_lt_spikesAligned = [];
norm_exd_ft_spikesAligned = [];
cmap = jet(256);
colormap(cmap(1:224,:));
L5Cells = [1:11];
for i = [1:11]
    norm_ft_spikesAligned(i,:) = sum(ft_spikesAligned{i})/size(ft_spikesAligned{i},1);
    norm_ft_spikesAligned(i,:) = norm_ft_spikesAligned(i,:)/max(norm_ft_spikesAligned(i,50:98));
    norm_lt_spikesAligned(i,:) = sum(lt_spikesAligned{i})/size(lt_spikesAligned{i},1);
    norm_lt_spikesAligned(i,:) = norm_lt_spikesAligned(i,:)/max(norm_lt_spikesAligned(i,:));
    norm_exd_ft_spikesAligned(i,:) = sum(exd_ft_spikesAligned{i})/size(exd_ft_spikesAligned{i},1);
    norm_exd_ft_spikesAligned(i,:) = norm_exd_ft_spikesAligned(i,:)/max(norm_exd_ft_spikesAligned(i,50:90));
end

% sort by weighted response time
weightedResponse=[];
for i = [1:11]
    alignby = norm_ft_spikesAligned(i,:)
    if ~isempty(find(alignby(56:75)>.75,1,'first'))
        tmp = find(alignby(56:75)>.75,1,'first');
    else [~, tmp] = max(alignby(50:75))
    end
    weightedResponse(i) = tmp;
end
[~,near_idx] = sort(weightedResponse([1:11]))


image(224*cat(1,norm_ft_spikesAligned(near_idx,25:100),zeros(1,76)))
title('Normalized first touch response')
set(gca,'xlim', [0 75],'xtick',[ 25 50 75 100],'xticklabel',[0 25 50 75])
ylabel('Cells')
xlabel('Time from touch (ms)')
box off

h_cb = colorbar
set(h_cb,'ytick',[0 .5 1])
% imagesc(cat(1,norm_lt_spikesAligned(1:17,:),zeros(1,201),norm_lt_spikesAligned_far))
% title('Normalized late touch response')
%
% set(gca,'xtick',[1 51 101 151],'xticklabel',[-50 0 50 100])
% ylabel('Cells')
% xlabel('Time from touch (ms)')
% box off
print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\StackedPCTHHeatmapNearFarL5'])

%% Fig 1 Grand PCTH mean

figure(1);clf;set(gcf,'paperposition',[0 0 2 1]);
gm={};
gmm = [];
pav = [];
for i =1:length(nearCells)
    tmpalgn = zeros(L{nearCells(i)}.k,3400);
    for j = 1:L{nearCells(i)}.k
        tmpalgn(j,:) = L{nearCells(i)}.R_ntk(1,[-400:2999]+find(L{nearCells(i)}.S_ctk(15,:,j)==1,1),j);
    end
    pav(i,:) = squeeze(mean(L{nearCells(i)}.S_ctk(15,[-400:2999]+find(L{nearCells(i)}.S_ctk(15,:,j)==1,1),:),3));
    gm{i} = (mean(tmpalgn,1));
    gmm(i,:) = gm{i};
    
end
smgmm= smooth(mean(gmm),10);
bar([0:10:3399],smgmm(5:10:3400)*1000,'k')
box off
set(gca,'xlim',[0 3400],'ylim',[0 6],'ytick',[0 2 4 6],'xtick',[400 1400 2400 3400],'xticklabel',[0 1 2 3])

print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\GrandMeanTrialL4PCTHC2'])
image((1-mean(pav))*255)
colormap(gray(256))
box off
print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\GrandMeanTrialL4PCTHC2_pole'])

%%
figure(1);clf;set(gcf,'paperposition',[0 0 2 1]);
gm={};
gmm = [];
pav = [];
for i =1:11
    tmpalgn =  zeros(L5{i}.k,3400);
    for j = 1:L5{i}.k
        tmpalgn(j,:) = L5{i}.R_ntk(1,[-400:2999]+find(L5{i}.S_ctk(15,:,j)==1,1),j);
     %   tmpalgn(j,:) = L5{i}.R_ntk(1,1:3400,j);
    end
    pav(i,:) = squeeze(mean(L5{i}.S_ctk(15,[-400:2999]+find(L5{i}.S_ctk(15,:,j)==1,1),:),3));
    gm{i} = (mean(tmpalgn,1));
    gmm(i,:) = gm{i};
    
end
smgmm= smooth(mean(gmm),10);
bar([0:10:3399],smgmm(5:10:3400)*1000,'k')
box off
set(gca,'xlim',[0 3400],'ylim',[0 50],'ytick',[0 2 4 6],'xtick',[400 1400 2400 3400],'xticklabel',[0 1 2 3])

print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\GrandMeanTrialL5PCTHC2'])
image((1-mean(pav))*255)
colormap(gray(256))
box off
print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\GrandMeanTrialL5PCTHC2_pole'])

%%
figure(1);clf;set(gcf,'paperposition',[0 0 2 1]);
gm={};
gmm = [];
pav = [];
for i =1:length(farCells)
    tmpalgn = zeros(L{farCells(i)}.k,3375);
    for j = 1:L{farCells(i)}.k
        tmpalgn(j,:) = L{farCells(i)}.R_ntk(1,[-400:2974]+find(L{farCells(i)}.S_ctk(15,:,j)==1,1),j);
    end
    pav(i,:) = squeeze(mean(L{farCells(i)}.S_ctk(15,[-400:2974]+find(L{farCells(i)}.S_ctk(15,:,j)==1,1),:),3));
    gm{i} = (mean(tmpalgn,1));
    gmm(i,:) = gm{i};
    
end
smgmm= smooth(mean(gmm),10);
bar([0:10:3374],smgmm(5:10:3375)*1000,'k')
box off
set(gca,'xlim',[0 3400],'ylim',[0 6],'ytick',[0 2 4 6],'xtick',[400 1400 2400 3400],'xticklabel',[0 1 2 3])

print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\GrandMeanTrialL4PCTHoutside'])
image((1-mean(pav))*255)
colormap(gray(256))
box off
print(gcf,'-depsc', ['C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\Fig2\Materials\GrandMeanTrialL4PCTHoutside_pole'])


