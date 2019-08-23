[rc] = numSubplots(length(selectedCells));


plotrow = rc(1);
plotcolumn = rc(2);

for rec = 1:length(U)
    
    %% plot raster go v no go
    figure(40);clf
    subplot(4,2,[1 2]);hold on
    
    for i = 1:U{rec}.k
        if sum(U{rec}.R_ntk(1,:,i))>0 & U{rec}.meta.trialType(i)==1%plot go Trials
            plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
        else
        end
    end
    set(gca,'ylim',[0 U{rec}.k]+1,'xlim',[0 4000],'xticklabel',[])
    title('Go')
    
    subplot(4,2,[3 4]);hold on
    for i = 1:U{rec}.k
        if sum(U{rec}.R_ntk(1,:,i))>0 & U{rec}.meta.trialType(i)==0%plot NOgo Trials
            plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
        else
        end
    end
    set(gca,'ylim',[0 U{rec}.k]+1,'xlim',[0 4000])
    title('Nogo')
    
    % plot average spike rate for go v nogo
    goidx=find(U{rec}.meta.trialType==1);
    nogoidx=find(U{rec}.meta.trialType==0);
    gospks=zeros(length(goidx),4000);
    nogospks=zeros(length(nogoidx),4000);
    swindow=50;
    
    for i=1:length(goidx)
        gospks(i,1:length(U{rec}.R_ntk(:,:,goidx(i))))=U{rec}.R_ntk(:,:,goidx(i));
    end
    %(i,1:length(U{rec}.R_ntk(:,:,goidx(i))))... because sometimes spikes
    %arent full 4000 in length
    for i=1:length(nogoidx)
        nogospks(i,1:length(U{rec}.R_ntk(:,:,nogoidx(i))))=U{rec}.R_ntk(:,:,nogoidx(i));
    end
    
    tmp=find(U{rec}.S_ctk(9,:,:)==1)/4000;
    d = tmp-floor(tmp); %keep only decimals
    avail=min(d)*4000; %pole available time in ms
    onset=mean(U{rec}.meta.poleOnset)*1000;
    avggospks=smooth(mean(gospks)*1000,swindow);
    avgnogospks=smooth(mean(nogospks)*1000,swindow);
    
    subplot(4,2,[5 8]);hold on
    plot(1:length(gospks),avggospks); %smoothing with 50ms window
    ylabel('spks/s')
    xlabel('time from trial start')
    
    hold on
    plot(1:length(gospks),avgnogospks,'r');
    hold on
    plot([avail avail],[0 max(avggospks)],'k:')
    plot([onset onset],[0 max(avggospks)],'k:')
    
    
    %% PLOT raster sorted by motor position 
    
    [~,motorsIdx] = sort(U{rec}.meta.motorPosition);
    touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
    
    %sorted motor position raster
    figure(888);hold on;
    for k = 1:length(motorsIdx)
        
         if sum(U{rec}.R_ntk(:,:,motorsIdx(k)))>0
            plot(find(U{rec}.R_ntk(:,:,motorsIdx(k))==1),k,'k.')
        else
        end
    end
    
    %sorted motor position heat map 
   spikes=squeeze(U{rec}.R_ntk(:,:,motorsIdx))';
   figure(20830);clf;imagesc(imgaussfilt(flipud(spikes),[1 10],'padding','replicate'))
   set(gca,'ytick',[],'xtick',0:1000:4000,'xlim',[0 4000])

   
   %sorted motor position TOUCH heat map 
   elements = [-25:50];

       allspikes = squeeze(U{rec}.R_ntk(:,:,:));
       fullMat = nan(length(motorsIdx),length(elements));
    for k = 1:length(motorsIdx)
         if sum(U{rec}.R_ntk(:,:,motorsIdx(k)))>0
             currTouchidx = find(ceil(touchIdx./4000)==motorsIdx(k));
               if ~isempty(currTouchidx)
                   spikeMat = repmat(touchIdx(currTouchidx),1,numel(elements)) + repmat(elements,numel(currTouchidx),1);
                   fullMat(k,:) = mean(allspikes(spikeMat),1);
               end
         end
    end
                   
    figure(889);clf            
    imagesc(fullMat)
    hold on; plot([find(elements ==0) find(elements ==0)],[0 length(motorsIdx)],'-.w')
    set(gca,'xtick', 1:25:length(elements),'xticklabel',-25:25:max(elements) ,'ytick',[])



end