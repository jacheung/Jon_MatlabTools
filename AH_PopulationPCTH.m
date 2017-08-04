L4L5X = [22:30 32 52:62];
ca_su = [3   4   5  6  7 10 28 54 55 67 73 82  85 88 89 90 91 100 101  103 104 39 41 42 48 56 60 71 74 81 84 86]
sfig2dir = 'C:\Users\shires.DTS\Dropbox\NoiseProject\NoiseManuscript\sFigures\SFig3_psthsETC\Materials\'
L4C2 = [38 48 39 41 47 43 45 51  9 10 11 12 49  4 46  1  7 14 15 16 17 44  3 42 50  5  6  8 40  2 13];

lt_spikesAligned = {};
ft_spikesAligned = {};
%%
for rec = 52:84
    ft_crop = L{rec}.S_ctk(9,:,:);
    lt_crop = L{rec}.S_ctk(12,:,:);
    
    
    
    ft_crop(1,end-150:end,:) = 0;
    lt_crop(1,end-150:end,:) = 0;
    ft_touchIdx = [find(ft_crop==1)];
    lt_touchIdx = [find(lt_crop==1)];

    ft_spikesAligned{rec} = zeros(numel(ft_touchIdx),251);
    lt_spikesAligned{rec} = zeros(numel(lt_touchIdx),251);
    spikes = squeeze(L{rec}.R_ntk);
 
    for i = 1:size(ft_spikesAligned{rec},1)
        ft_spikesAligned{rec}(i,:) = spikes(ft_touchIdx(i)+[-100:150]);
    end
    
      for i = 1:size(lt_spikesAligned{rec},1)
        lt_spikesAligned{rec}(i,:) = spikes(lt_touchIdx(i)+[-100:150]);
      end
    
    figure(1);clf;set(gcf,'paperposition',[0 0 2 .75])
    bar(-100:150,sum(ft_spikesAligned{rec})/numel(ft_touchIdx),'k');hold on
    plot([0 0 ],[0 .2],'r')
    set(gca,'xlim',[-50 100],'ylim',[0 .4])
    axis off
    box off
    % print(gcf,'-depsc', [sfig2dir 'FT_PCTH_' num2str(rec)])
    
     figure(2);clf;set(gcf,'paperposition',[0 0 2 .75])
    bar(-100:150,sum(lt_spikesAligned{rec})/numel(lt_touchIdx),'k')

    set(gca,'xlim',[-50 100],'ylim',[0 .2])
    axis off
    box off
    % print(gcf,'-depsc', [sfig2dir 'LT_PCTH_' num2str(rec)])
end

%% Latency
for i = 1:82
    spthresh1(i) = mean(mean(ft_spikesAligned{i}(:,1:100)))+3*std(mean(ft_spikesAligned{i}(:,1:50),2));
    spthresh2(i) = mean(mean(lt_spikesAligned{i}(:,75:100)))+3*std(mean(lt_spikesAligned{i}(:,1:50),2));
    if isempty([find(mean(lt_spikesAligned{i}(:,102:151))>spthresh2(i),1) find(mean(ft_spikesAligned{i}(:,102:151))>spthresh1(i),1)])
       onsetLatency(i) = NaN;
    elseif isempty([find(mean(ft_spikesAligned{i}(:,102:151))>spthresh1(i),1)])
            onsetLatency(i) = find(mean(lt_spikesAligned{i}(:,102:151))>spthresh2(i),1)
    elseif isempty([find(mean(lt_spikesAligned{i}(:,102:151))>spthresh2(i),1)])
        onsetLatency(i) = find(mean(ft_spikesAligned{i}(:,102:151))>spthresh1(i),1)
    else
                        onsetLatency(i) = min([find(mean(ft_spikesAligned{i}(:,102:151))>spthresh1(i),1) find(mean(lt_spikesAligned{i}(:,102:151))>spthresh2(i),1)])

    end
end
onsetLatency(4)= 9;
onsetLatency(17)= 9;
onsetLatency(28)= NaN;
onsetLatency(39)= 5;
onsetLatency(41)= 6;
onsetLatency(42)= 6;
onsetLatency(49)= 7;
onsetLatency(59)= 6;
onsetLatency(61)= 6;
onsetLatency(72)= 4;
onsetLatency(74)= 4;
onsetLatency(78)= 4;
onsetLatency(79)= 16;
onsetLatency(80)= 4;
onsetLatency(81)= 11;

onsetLatencyOut = [18 20 21 26];

figure(1);clf;set(gcf,'paperposition',[0 0 2 1]);hold on
edges=(0.5:30.5);
hist_lat = histc(onsetLatency([1:17 38:51]),edges);
hist_latout = histc(onsetLatency([22:30 32]),edges);
hist_latL5 = histc(onsetLatency([52:62]),edges);

bar(1:31,hist_lat,'r','edgecolor','none')
bar(1:31,hist_latout,'facecolor',[.5 .5 .5],'edgecolor','none')
bar(1:31,hist_latL5,'facecolor','none','edgecolor',[0 .5 0])
set(gca,'xlim',[0 30],'ylim',[0 6],'ytick',[0 3 6])

print(gcf,'-depsc', [sfig2dir 'Latency_Hist'])


%% Phase

for rec = ca
    
 figure(1);clf;set(gcf,'paperposition',[0 0 2 .6]);hold on
 plot([-3.5 -3.5],[0 .5],'k','linewidth',2)

bar(r.phaseAlig(rec).meanBin.phase, r.phaseAlig(rec).meanBin.spk,'k')
 % plot(linspace(0,10,120),r.phase.fit(rec,:)/size(r.phaseAlig(rec).spk,1)*10000,'color',[.5 .5 .5]);
plot(linspace(-pi,pi,120),r.phase.fit(rec,:),'color',[.5 .5 .5]);
set(gca,'xlim',[-3.6 pi])
 axis off
 print(gcf,'-depsc', [sfig2dir 'Phase_' num2str(rec)])

end

%%
for rec = ca
    
 figure(1);clf;set(gcf,'paperposition',[0 0 2 .6]);hold on
 bar(.5:1:9.5,sum(reshape(sum(r.phaseAlig(rec).spk),12,10))/(12*size(r.phaseAlig(rec).spk,1))*1000,'k')
% plot(linspace(0,10,120),r.phase.fit(rec,:)/size(r.phaseAlig(rec).spk,1)*10000,'color',[.5 .5 .5]);
% plot(linspace(0,10,120),r.phase.fit(rec,:)/(size(r.phaseAlig(rec).spk,1))*1000,'color',[.5 .5 .5]);
set(gca,'xlim',[0 10])
 
 print(gcf,'-depsc', [sfig2dir 'Phase_' num2str(rec)])

end


%% Phase

for rec = 1
 figure(1);clf;set(gcf,'paperposition',[0 0 2 .6]);hold on
 bar(),12,10)))%/size(r.phaseAlig(rec).spk,1)*1000,'k')
% plot(linspace(0,10,120),r.phase.fit(rec,:)/size(r.phaseAlig(rec).spk,1)*10000,'color',[.5 .5 .5]);
 plot(linspace(0,10,120),r.phase.fit(rec,:),'color',[.5 .5 .5]);
set(gca,'xlim',[0 10])
 axis off
 print(gcf,'-depsc', [sfig2dir 'Phase_' num2str(rec)])

end

%%

figure(1);clf;set(gcf,'paperposition',[0 0 2 40])
for i = 1:31
    subplot(31,1,i)
        bar(-100:150,sum(ft_spikesAligned{L4C2(i)})/size(ft_spikesAligned{L4C2(i)},1),'k','edgecolor','none');hold on
        plot([0 0 ],[0 .2],'r')
        set(gca,'xlim',[-50 100])
 
%        subplot(31,2,i+31)
%         bar(-100:150,sum(lt_spikesAligned{L4C2(i)})/size(lt_spikesAligned{L4C2(i)},1),'k','edgecolor','none');hold on
%     plot([0 0 ],[0 .2],'r')
    axis off
end
 print(gcf,'-depsc', [sfig2dir '1stTouchL4C2sortedPCTH' num2str(rec)])
 %%
 figure(1);clf;set(gcf,'paperposition',[0 0 2 40])
 for i = 1:31
    subplot(31,1,i)
%         plot([onsetLatency(L4C2(i))-.5 resL4(L4C2(i)).idxT+.5], [1 1]*max(sum(lt_spikesAligned{L4C2(i)})/size(lt_spikesAligned{L4C2(i)},1)), '-','color',[221 204 119]/255)
        plot([6.5 resL4(L4C2(i)).idxT+6.5], [1 1]*max(sum(lt_spikesAligned{L4C2(i)})/size(lt_spikesAligned{L4C2(i)},1)), '-','color',[221 204 119]/255)

        hold on

        bar(-100:150,sum(lt_spikesAligned{L4C2(i)})/size(lt_spikesAligned{L4C2(i)},1),'k','edgecolor','none');
        plot([0 0 ],[0 .2],'r')
        set(gca,'xlim',[-50 100])

        %        subplot(31,2,i+31)
%         bar(-100:150,sum(lt_spikesAligned{L4C2(i)})/size(lt_spikesAligned{L4C2(i)},1),'k','edgecolor','none');hold on
%     plot([0 0 ],[0 .2],'r')
    axis off
end
 print(gcf,'-depsc', [sfig2dir 'LateTouchL4C2sortedPCTH' num2str(rec)])
