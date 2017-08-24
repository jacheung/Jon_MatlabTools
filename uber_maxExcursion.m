%**************************************************************************
% Purpose of this code is to quantify search strategy. 
%
% This is done by looking at the max excursion of the whisker at each
% protraction cycle and identify its distance traveled relative to the
% discrimination boundary. 
%
% Input: Uber Array 
%
%
%
%**************************************************************************

close all; clearvars -except U D BV
%% Parameters
% selectedMotorT filter to choose which type of trials to use to look at
% max excursion. Options are:
    % 'nogo10k' 'allgo' 'allnogo' 'DB10k' 'FA' 'CR'    
selectedMotorT = 'nogo10k'  ;
for rec = 1:length(U)
    array=U{rec};
    %% Find Theta Required to Contact Pole
    [~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(array);
    
    varx=[1:6];
    for f=1:length(varx)
        [~, ~,~, ~, ~,~, E{rec}.go{f} , E{rec}.nogo{f}] = assist_vardistribution(U{rec},varx(f),prelixGo,prelixNoGo,[-25:0],[5:25]);
    end
    
    if strcmp(array.meta.layer,'D') %DISCRETE 
        goprot=cell2mat(E{rec}.go{5})<0; nogoprot=cell2mat(E{rec}.nogo{5})<0;
        gothetas=cell2mat(E{rec}.go{1}); nogothetas=cell2mat(E{rec}.nogo{1});
        gothetamean = mean(gothetas(goprot));  nogothetamean = mean(nogothetas(nogoprot)); 
    else    
     
    %FOR SM AND BV
    gomotor = U{rec}.meta.motorPosition(U{rec}.meta.trialType==1);
    nogomotor = U{rec}.meta.motorPosition(U{rec}.meta.trialType==0);
    motortheta = [];
    motorthetamean=[];
    
    
    for i = 1:length(gomotor)
        selectedgo = E{rec}.go{5}{i}<0; %selecting only protraction touches
        gotmp=[gomotor(i) mean(E{rec}.go{1}{i}(selectedgo))];
        motorthetamean = [motorthetamean ; gotmp];
    end
    
    for k = 1:length(nogomotor)
        selectednogo = E{rec}.nogo{5}{k}<0;
        nogotmp=[nogomotor(k) mean(E{rec}.nogo{1}{k}(selectednogo))];
        motorthetamean = [motorthetamean ; nogotmp]    ;
    end
    
    motorthetamean(isnan(motorthetamean(:,2)),:)=[];
    motorthetamean=sortrows(motorthetamean);
    if rec ==2
        motorthetamean(1,:)=[];
    end
    
    
    [E{rec}.coeff, ~ , E{rec}.mu] = polyfit(motorthetamean(:,1),motorthetamean(:,2),2);
    f = polyval(E{rec}.coeff,motorthetamean(:,1),[],E{rec}.mu);
    
    figure(70); subplot(3,4,rec);
    scatter(motorthetamean(:,1),motorthetamean(:,2));
    hold on; plot(motorthetamean(:,1),f)
    xlabel('Motor Position (Distal - Proximal)'); ylabel('Theta at base')
    
    mp=U{rec}.meta.motorPosition;
    thetareqDB=polyval(E{rec}.coeff,mp,[],E{rec}.mu);
    dbtheta = max(thetareqDB(U{rec}.meta.trialType==1));
    hold on; plot([0 180000],[dbtheta dbtheta],'r')
    set(gca,'xlim',[min(motorthetamean(:,1))-10000 180000])
    end
    
    
    %% Find max theta excursion on each whisk cycle
    [objmask]= assist_touchmasks(array);
    masktmp = objmask.samplingp; %
%     mask=ones(size(masktmp));
    %mask(~isnan(masktmp))=NaN;
    mask = objmask.samplingp;
    
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
    
    %% test a trial to make sure that theta is aligned right
        trial = 1; %shifted by 1 (ex. trial=0 ..> trial =1)
        figure(trial+1);clf;plot(theta(:,trial+1));
        xlabel('Time from Trial Start (ms)');ylabel('Whisker Position')
        
        validx=peakidx(floor(peakidx/U{rec}.t)==trial);
        excurx=round(((validx/U{rec}.t)-trial)*U{rec}.t);
        for i = 1:length(excurx)
            hold on; scatter(excurx(i),theta(excurx(i),trial+1),'ro')
        end
%     %
%     set(gca,'xlim',[550 1300])
%     gos=U{1}.meta.trialType==1;
%     nogos = U{1}.meta.trialType==0;
    
    %% max excursion
    if strcmp(array.meta.layer,'D')
        dbtheta = mean([gothetamean nogothetamean]);
        tidx = find(array.meta.trialType==0);
    else
        mp=U{rec}.meta.motorPosition;
        thetareqDB=polyval(E{rec}.coeff,mp,[],E{rec}.mu);
        
        %changing DB based on whether SM or BV
        if strcmp(array.meta.layer,'SM')
            dbthetatmp = max(thetareqDB(U{rec}.meta.trialType==1));
            dbnogotheta = min(thetareqDB(U{rec}.meta.trialType==0));
            dbtheta = mean([dbthetatmp dbnogotheta]);
        elseif strcmp(array.meta.layer,'BV')
            dbtheta = max(thetareqDB(U{rec}.meta.trialType==1));
        end
        
        
        switch selectedMotorT
            case 'allgo'
                filt2=mp<min(mp(U{rec}.meta.trialType)); %elim less than min go
                mp(filt2)=NaN;
            case 'allnogo'
                filt = U{rec}.meta.trialType==0; %elim all gotrials
                mp(~filt)=NaN;
            case 'nogo10k'
                filt= mp>min(mp)+10000; %farthest nogo filter
                mp(filt)=NaN;
            case 'DB10k'
                filt = mp>=min(mp(U{rec}.meta.trialType)) +10000; %elim greater than min go motor+10k
                mp(filt)=NaN;
                
            case 'CR'
                filt = U{rec}.meta.trialType==0 .* U{rec}.meta.trialCorrect==1; %all CR trials
                mp(~filt)=NaN;
            case 'FA'
                filt = U{rec}.meta.trialType==0 .* U{rec}.meta.trialCorrect==0; %all FA trials
                mp(~filt)=NaN;
        end
        tidx = find(~isnan(mp)==1)-1;
    end
    xplor=[];
    for j=1:length(tidx)
        normidx=peakidx(floor(peakidx/array.t)==tidx(j));
        %xplor=[xplor;repmat(thetareq(j+1),length(normidx),1)-theta(normidx)];
        %xplor with normalizing to thetareq for each trial
        xplor=[xplor;repmat(dbtheta,length(normidx),1)-theta(normidx)]; %normalization to farthest go trial
    end
    distmean(rec)=mean(xplor);diststd(rec)=std(xplor);distmed(rec)=nanmedian(xplor);
    gxplor{rec}=xplor;
    
    %% Plotting Features for Search
    ranges=[round(min(xplor)-2):round(max(xplor)+2)];
    vals=histc(xplor,ranges);
    figure(17);subplot(2,5,rec)
    bar(ranges,vals/(sum(vals)),'k');
    hold on; plot([0 0],[0 1],'r','LineWidth',1)
    xlabel('Theta from Discrimination Boundary');ylabel('Proportion of Trial')
    set(gca,'xlim',[-50 50],'ylim',[0 .1]);
    alpha(.8)
    
    figure(18);subplot(2,5,rec)
    bar(ranges,vals/(sum(vals)),'k');
    hold on; plot([0 0],[0 1],'r','LineWidth',1)
    xlabel('Theta from Discrimination Boundary');ylabel('Proportion of Trial')
    set(gca,'xlim',[-50 0],'ylim',[0 .1]);
    alpha(.8)
    
end
%% PLOTTING MEDIANS
figure(376);clf;scatter(1:length(U),distmed,'k','linewidth',8)
hold on; plot([0 length(U)],[0 0], '-.ok')
set(gca,'ylim',[-20 20],'yaxislocation','right')
xlabel('Mouse Number'); ylabel('Search Bias')
[~,newIdx] = sort(abs(distmed));
for k = 1:length(U) %make new J array with U sorted by lowest median of search 
    J{k} = U{newIdx(k)};
end
 
camroll(-90)
%% PLOTTING UBER HISTO
ranges = [-50:50];
ghisto= histc(cell2mat(gxplor'),ranges);
    figure(58);
    bar(ranges,ghisto/(sum(ghisto)),'k');
    hold on; plot([0 0],[0 1],'r','LineWidth',1)
    xlabel('Theta from Discrimination Boundary');ylabel('Proportion of Trial')
    set(gca,'xlim',[-50 50],'ylim',[0 .1]);
    alpha(.8)
    nanmedian(cell2mat(gxplor'));
    nanstd(cell2mat(gxplor'));
%% K means clustering to identify number of search strategies
%probably want to do more initializations and see which cluster gives lowest cost
numClusters = 2;
clusters = kmeans(distmed',numClusters);
numMouse = 1:length(clusters);
[~,newIdx] = sort(distmed); %sort indices based on lowest median 
for k = 1:length(U) %make new J array with U sorted by lowest median of search 
    J{k} = U{newIdx(k)};
end


figure(10);clf;
%scatter(1:length(clusters),distmed);
colors = [{'b'},{'r'},{'g'},{'m'}];
for k = 1:numClusters
    hold on; scatter(numMouse(clusters==k),distmed(clusters==k),colors{k})
end

xlabel('Mouse Number'); ylabel('Median of Search')
legend('High Search Bias','Low Search Bias')

