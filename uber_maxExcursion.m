
 figure(17);clf
 figure(70);clf
 figure(18);clf
 distmean=[];
 diststd=[];
for rec = 1:length(U)
    array=U{rec};
    [~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(array);
    varx=[1:6];
    for f=1:length(varx)
        [~, ~,V{rec}.hit{f}, V{rec}.miss{f}, V{rec}.FA{f}, V{rec}.CR{f}, V{rec}.go{f} , V{rec}.nogo{f}] = assist_vardistribution(U{rec},varx(f),prelixGo,prelixNoGo,[-25:0],[5:25]);
    end
    
    gomotor = U{rec}.meta.motorPosition(U{rec}.meta.trialType==1);
    nogomotor = U{rec}.meta.motorPosition(U{rec}.meta.trialType==0);
    motortheta = [];
    motorthetamean=[];
    
    % make raw values of theta instead of mean.
    %     for i = 1:length(gomotor)
    %     selectedgo = V{rec}.go{5}{i}<0; %selecting only protraction touches
    %     gotmp=[repmat(gomotor(i),sum(selectedgo),1) V{rec}.go{1}{i}(selectedgo)'];
    %     motortheta = [motortheta ; gotmp];
    %     end
    %     for k = 1:length(nogomotor)
    %     selectednogo = V{rec}.nogo{5}{k}<0;
    %     nogotmp=[repmat(nogomotor(k),sum(selectednogo),1) V{rec}.nogo{1}{k}(selectednogo)'];
    %     motortheta = [motortheta ; nogotmp]    ;
    %     end
    %
    
    for i = 1:length(gomotor)
        selectedgo = V{rec}.go{5}{i}<0; %selecting only protraction touches
        gotmp=[gomotor(i) mean(V{rec}.go{1}{i}(selectedgo))];
        motorthetamean = [motorthetamean ; gotmp];
    end
    
    for k = 1:length(nogomotor)
        selectednogo = V{rec}.nogo{5}{k}<0;
        nogotmp=[nogomotor(k) mean(V{rec}.nogo{1}{k}(selectednogo))];
        motorthetamean = [motorthetamean ; nogotmp]    ;
    end
    
    motorthetamean(isnan(motorthetamean(:,2)),:)=[];
    motorthetamean=sortrows(motorthetamean);
        if rec ==2 
            motorthetamean(1,:)=[];
        end
   

    [V{rec}.coeff, ~ , V{rec}.mu] = polyfit(motorthetamean(:,1),motorthetamean(:,2),2);
    f = polyval(V{rec}.coeff,motorthetamean(:,1),[],V{rec}.mu);
     
    figure(70); subplot(3,4,rec);
    scatter(motorthetamean(:,1),motorthetamean(:,2));
    hold on; plot(motorthetamean(:,1),f)
    xlabel('Motor Position (Distal - Proximal)'); ylabel('Theta at base')
    
    mp=U{rec}.meta.motorPosition;
    thetareqDB=polyval(V{rec}.coeff,mp,[],V{rec}.mu);
    dbtheta = max(thetareqDB(U{rec}.meta.trialType==1));
    hold on; plot([0 180000],[dbtheta dbtheta],'r')
    set(gca,'xlim',[min(motorthetamean(:,1))-10000 180000])
    
    
    
    %Find max theta excursion on each whisk cycle 
    [objmask]= assist_touchmasks(U{rec});  
    masktmp = objmask.samplingp; %set all pole avail time to NaN
    mask=ones(size(masktmp));
    mask(~isnan(masktmp))=NaN;
    
    amp_mask = ones(size(squeeze(U{rec}.S_ctk(3,:,:))));
    amp_mask(U{rec}.S_ctk(3,:,:)<2.5) = NaN; %amplitude mask used to filter out only periods of high whisking
    phase = squeeze(U{rec}.S_ctk(5,:,:));
    amp = squeeze(U{rec}.S_ctk(3,:,:));
    theta = squeeze(U{rec}.S_ctk(1,:,:));
    selectedPhase = mask.*amp_mask.*phase; %can add a mask here to adjust what time period to look at (i.e. pole avail to end trila mask)
    phases = find(selectedPhase>-.1 & selectedPhase<.1); %find all phases within this window
    
    samp_r=[-20:20];%look -20:20ms around those found phases
    phases=phases(phases+max(samp_r)<numel(theta));
    phases=phases(phases+min(samp_r)>0);
    
    pidx=repmat(phases,1,41)+repmat(samp_r,length(phases),1); %build indices
    maxidx=zeros(1,size(pidx,1));
    for f = 1:size(pidx,1)
        [~,maxtmp]=max(theta(pidx(f,:)));
        maxidx(f)=samp_r(maxtmp);
    end
    
    peakidx=unique(phases+maxidx'); %this is the idx for max excursion.
%     
%     tmpDB(rec)=dbtheta;
%     rest(rec)=nanmean(theta(amp<=0.05));
    
     %% test a trial to make sure you are looking at the right one
%     trial = 69; %shifted by 1 (ex. trial=0 ..> trial =1)
%     figure(trial+1);clf;plot(theta(:,trial+1));
%     validx=peakidx(floor(peakidx/4000)==trial);
%     excurx=round(((validx/4000)-trial)*4000);
%     for i = 1:length(excurx)
%         hold on; plot([excurx(i) excurx(i)],[-15 60],'r')
%     end
%     


    
    %% max excursion
    
    mp=U{rec}.meta.motorPosition;
    thetareqDB=polyval(V{rec}.coeff,mp,[],V{rec}.mu);
    dbtheta = max(thetareqDB(U{rec}.meta.trialType==1));
    
          filt= mp>min(mp)+10000; %farthest nogo filter
           mp(filt)=NaN;
%            filt = mp>=min(mp(U{rec}.meta.trialType)) +10000; %elim greater than min go motor+10k
%            mp(filt)=NaN;
%            filt2=mp<min(mp(U{rec}.meta.trialType)); %elim less than min go
%            mp(filt2)=NaN;
%            
%            filt = U{rec}.meta.trialType==0; %elim all gotrials
%            mp(~filt)=NaN;
           
%            filt = U{rec}.meta.trialType==0 .* U{rec}.meta.trialCorrect==1; %all CR trials
%            mp(~filt)=NaN;
           
%            filt = U{rec}.meta.trialType==0 .* U{rec}.meta.trialCorrect==0; %all FA trials
%            mp(~filt)=NaN;
           
    tidx = find(~isnan(mp)==1)-1;
    xplor=[];
    for j=1:length(tidx)
        normidx=peakidx(floor(peakidx/4000)==tidx(j));
        %xplor=[xplor;repmat(thetareq(j+1),length(normidx),1)-theta(normidx)];
        %xplor with normalizing to thetareq for each trial
        xplor=[xplor;repmat(dbtheta,length(normidx),1)-theta(normidx)]; %normalization to farthest go trial
    end
    ranges=[round(min(xplor)-2):round(max(xplor)+2)];
    vals=histc(xplor,ranges);
    figure(17);subplot(2,4,rec)
    bar(ranges,vals/(sum(vals)),'b');
    hold on; plot([0 0],[0 1],'r','LineWidth',1)
    xlabel('ThetaREQ - MAX Excursion');ylabel('Proportion of Trials')
    set(gca,'xlim',[-50 50],'ylim',[0 .1]);
    alpha(.8)
    
    figure(18);subplot(2,4,rec)
    bar(ranges,vals/(sum(vals)),'b');
    hold on; plot([0 0],[0 1],'r','LineWidth',1)
    xlabel('ThetaREQ - MAX Excursion');ylabel('Proportion of Trials')
    set(gca,'xlim',[-50 0],'ylim',[0 .1]);
    alpha(.8)
    
    distmean(rec)=mean(xplor);diststd(rec)=std(xplor);distmed(rec)=median(xplor);
end
