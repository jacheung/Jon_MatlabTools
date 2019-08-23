%Finding 'start of whisk'

close all
whiskAmpval = 3; %what amplitude value threshold to count as whisk start
thetaThresh = 5; %what value above rest to count as whisk start
window = [-25 25];
% tomask = 'none';
var2use = 'touchOnset';
figure(290);clf
figure(68);clf
PopMods = zeros(length(U),2);


for rec = 1:length(U)
    
    amps = squeeze(U{rec}.S_ctk(3,:,:));
    thetas = squeeze(U{rec}.S_ctk(1,:,:));
    %in case whiskOnset val is right at beginning of trial, this will block
    %those out as we can't build preWhisk histogram.
    mask = ones(size(amps));
    mask(1:abs(window(1)),:)=NaN;
    tOnsetmask = ones(size(amps));
    tOnsetmask(U{rec}.t-window(2):U{rec}.t,:) = NaN; %not allowing touches that dont have at least xms after onset to be defined.
    firstOnsets = squeeze(U{rec}.S_ctk(9,:,:));
    lateOnsets = squeeze(U{rec}.S_ctk(12,:,:));
    
    firstOnsets = firstOnsets.*tOnsetmask;
    lateOnsets = lateOnsets.*tOnsetmask;
    thetas = thetas.*mask;
    amps = amps.*mask;
    
    
    %     if strcmp(tomask,'yes')
    %         %use this to choose mask:
    %         %availtoend_mask, avail_mask, touchEx_mask, firsttouchEx_mask
    %         [objmask]= assist_touchmasks(U{rec});
    %         mask = objmask.samplingp; %only look within sampling p
    %         thetas = thetas.*mask;
    %         amps = amps.*mask;
    %     end
    
    
    %finding first index of whiskAmp crossing using theta or amp
    ampList = intersect(find(amps>whiskAmpval-.1),find(amps<whiskAmpval+.1));
    crossList = ampList+1;
    [~,x,~] = intersect(ampList,crossList);
    ampList(x) = [];
    restingWhisk = nanmedian(thetas(:));
    thetaList = find(thetas>restingWhisk+5);
    touchOnsetList = [find(firstOnsets==1);find(lateOnsets==1)];
    
    
    if strcmp(var2use,'theta')
        whiskDeterm = thetaList;    
        %ADDITIONAL: finding first instance of >whiskAmpVal in each trial if you
        %dont want to use ALL amps
        [~,p]=unique(floor(whiskDeterm/4000));
        whiskDeterm=whiskDeterm(p);
    elseif strcmp(var2use,'amp')
        whiskDeterm = ampList;
        %ADDITIONAL: finding first instance of >whiskAmpVal in each trial if you
        %dont want to use ALL amps
        [~,p]=unique(floor(whiskDeterm/4000));
        whiskDeterm=whiskDeterm(p);
    elseif strcmp(var2use,'touchOnset')
        whiskDeterm = touchOnsetList;        
    else
        error('select theta or amp or touchOnset for which to determine whisk onset')
    end
    
    
    %building spikes around variable index
    spikes = squeeze(U{rec}.R_ntk);
    spikesAligned = zeros(numel(whiskDeterm),numel(window(1):window(2)));
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(whiskDeterm(i)+[window(1):window(2)]);
    end
    
    %ANOVA to compare whether means are sig different from one another
    alignedPre = sum(spikesAligned)./numel(whiskDeterm);
    preWhisk = alignedPre(1:0-(window(1)));
    postWhisk = alignedPre(2-(window(1)):end);
    [p,tbl,stats] = anova1([preWhisk' postWhisk']); %anova1
    comp = multcompare(stats);
    
    modIdx=(mean(postWhisk) - mean(preWhisk))/(mean(postWhisk)+mean(preWhisk));
    
    
    %use this to plot first trial of each rec and see which amplitude points you're mapping
    figure(290);hold on;subplot(5,4,rec)
    plot(U{rec}.S_ctk(3,:,1),'k')
    hold on; plot(U{rec}.S_ctk(1,:,1),'b')
    if ~isempty(whiskDeterm(whiskDeterm<4000))
    hold on; plot([whiskDeterm(whiskDeterm<4000) whiskDeterm(whiskDeterm<4000)],[-10 60],'r')
    end
    xlabel('time from trial start (ms)')
    if rec == 4
        legend('amp','theta','OnsetIdx')
    end
    
    %Plotting histograms
    figure(68);hold on;subplot(5,4,rec)
    bar(window(1):window(2),sum(spikesAligned)/numel(whiskDeterm),'k');
    if comp(end)<.05
        if modIdx<0
            if strcmp(var2use,'touchOnset')
                text(window(2)*.5,max(sum(spikesAligned)/numel(whiskDeterm))*.75,'TouchOFF','color','red')
            else
                text(window(2)*.5,max(sum(spikesAligned)/numel(whiskDeterm))*.75,'WhOFF','color','red')
            end
        elseif modIdx>0
            if strcmp(var2use,'touchOnset')
                text(window(2)*.5,max(sum(spikesAligned)/numel(whiskDeterm))*.75,'TouchON','color','blue')
            else
                text(window(2)*.5,max(sum(spikesAligned)/numel(whiskDeterm))*.75,'WhON','color','blue')
            end
        end
    end
    
    PopMods(rec,:) = [comp(end) modIdx];
    set(gca,'xlim',window,'xtick',[window(1):25:window(2)]);
%     xlabel('Time from whisk onset (ms)')
%     ylabel('spks / ms')
    
    
end