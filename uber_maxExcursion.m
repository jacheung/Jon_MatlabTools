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
% SM=SM([1 3:10]);
clearvars -except U D BV SM N 
%% Parameters
% selectedMotorT filter to choose which type of trials to use to look at
% max excursion. Options are:
% 'nogo10k' 'allgo' 'allnogo' 'DB10k' 'FA' 'CR'

selectedMotorT = 'nogo10k'  ;
distmed=[];
popGo = [];
popNogo = [];
for rec = 1:length(U)
    array=U{rec};
    % Find Theta Required to Contact Pole
    [~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(array);
    
    
    if strcmp(array.meta.layer,'D') %DISCRETE
        varx=[1:6]; %not 1 to 7 b/c no follicleXY
        for f=1:length(varx)
            [~, ~,~, ~, ~,~, E{rec}.go{f} , E{rec}.nogo{f}] = assist_vardistribution(U{rec},varx(f),prelixGo,prelixNoGo,[-25:0],[5:25]);
        end
        
        goprot=cell2mat(E{rec}.go{5})<0; nogoprot=cell2mat(E{rec}.nogo{5})<0;
        gothetas=cell2mat(E{rec}.go{1}); nogothetas=cell2mat(E{rec}.nogo{1});
        gothetamean = mean(gothetas(goprot));  nogothetamean = mean(nogothetas(nogoprot));
        
        
        figure(70); subplot(2,5,rec);
            scatter(zeros(length(nogothetas),1)-1,nogothetas,'k');
            hold on; scatter(ones(length(gothetas),1),gothetas,'k');
            xlabel('Normalized Motor Positions'); ylabel('Theta at base')
            set(gca,'xlim',[-1.2 1.2])
            
            
    
    else
        varx=[1:7];
        for f=1:length(varx)
            [~, ~,~, ~, ~,~, E{rec}.go{f} , E{rec}.nogo{f}] = assist_vardistribution(U{rec},varx(f),prelixGo,prelixNoGo,[-25:0],[5:25]);
        end
        
        %FOR SM AND BV
        gomotor = U{rec}.meta.motorPosition(U{rec}.meta.trialType==1);
        nogomotor = U{rec}.meta.motorPosition(U{rec}.meta.trialType==0);
        motorthetamean=[];
        motorFollicleMean = [];
        
        for i = 1:length(gomotor)
            selectedgo = E{rec}.go{5}{i}<0; %selecting only protraction touches
            gotmp=[gomotor(i) mean(E{rec}.go{1}{i}(selectedgo))];
            if ~isempty(E{rec}.go{7}{i}(selectedgo,:))
                goFoltmp = [gomotor(i) mean(E{rec}.go{7}{i}(selectedgo,:),1)];
                motorthetamean = [motorthetamean ; gotmp];
                motorFollicleMean = [motorFollicleMean ; goFoltmp];
            end
        end
        
        for k = 1:length(nogomotor)
            selectednogo = E{rec}.nogo{5}{k}<0;
            nogotmp=[nogomotor(k) mean(E{rec}.nogo{1}{k}(selectednogo))];
            if ~isempty(E{rec}.nogo{7}{k}(selectednogo,:))
                nogoFoltmp = [nogomotor(k) mean(E{rec}.nogo{7}{k}(selectednogo,:),1)];
                motorthetamean = [motorthetamean ; nogotmp];
                motorFollicleMean = [motorFollicleMean ; nogoFoltmp];
            end
        end
        
        motorthetamean(isnan(motorthetamean(:,2)),:)=[];
        motorFollicleMean(isnan(motorFollicleMean(:,2)),:)=[];
        motorthetamean=sortrows(motorthetamean);
        motorFollicleMean=sortrows(motorFollicleMean);
        
        %         if rec ==2 %error in SM2
        %             motorthetamean(1,:)=[];
        %         end
        %check to make sure there are at least some touches in nogo
        %positions
        minGo = min(U{rec}.meta.motorPosition(U{rec}.meta.trialType==1));
        nogoMotors = find(motorthetamean(:,1)<minGo);
        
        if ~isempty(nogoMotors)
            
            [E{rec}.coeff, ~ , E{rec}.mu] = polyfit(motorthetamean(:,1),motorthetamean(:,2),2);
            f = polyval(E{rec}.coeff,motorthetamean(:,1),[],E{rec}.mu);
            [E{rec}.FolXCoeff, ~ , E{rec}.FolXMu] = polyfit(motorFollicleMean(:,1),motorFollicleMean(:,2),2);
            g = polyval(E{rec}.FolXCoeff,motorFollicleMean(:,1),[],E{rec}.FolXMu);
            [E{rec}.FolYCoeff, ~ , E{rec}.FolYMu] = polyfit(motorFollicleMean(:,1),motorFollicleMean(:,3),2);
            h = polyval(E{rec}.FolYCoeff,motorFollicleMean(:,1),[],E{rec}.FolYMu);
            
            
            figure(70); subplot(2,5,rec);
            scatter(motorthetamean(:,1),motorthetamean(:,2),'k');
            hold on; plot(motorthetamean(:,1),f,'k')
            xlabel('Normalized Motor Positions'); ylabel('Theta at base')
            
            mp=U{rec}.meta.motorPosition;
            thetareqDB=polyval(E{rec}.coeff,mp,[],E{rec}.mu);
            dbtheta = max(thetareqDB(U{rec}.meta.trialType==1));
            hold on; plot([0 180000],[dbtheta dbtheta],'-.k')
%             set(gca,'xlim',[U{rec}.meta.ranges],'xtick',[U{rec}.meta.ranges(1) mean(U{rec}.meta.ranges) U{rec}.meta.ranges(2)],'xticklabel',[-1 0 1])
            set(gca,'xtick',[U{rec}.meta.ranges(1) mean(U{rec}.meta.ranges) U{rec}.meta.ranges(end)],'xticklabel',[-1 0 1])
            
            
            
            figure(80);subplot(2,5,rec)
            
            scatter(motorFollicleMean(:,1),motorFollicleMean(:,2)./33,'k');
            %         hold on; plot(motorFollicleMean(:,1),g,'k'); %plot poly fit for follicle X
            ylabel('Follicle X')
            yyaxis right
            scatter(motorFollicleMean(:,1),motorFollicleMean(:,3)./33,'r');
            %         plot(motorFollicleMean(:,1),h,'r-'); %plot poly fit for follicle Y
            ylabel('Follicle Y')
            xlabel('Normalized Motor Positions');
%             set(gca,'xlim',[U{rec}.meta.ranges],'xtick',[U{rec}.meta.ranges(1) mean(U{rec}.meta.ranges) U{rec}.meta.ranges(end)],'xticklabel',[-1 0 1])
            set(gca,'xtick',[U{rec}.meta.ranges(1) mean(U{rec}.meta.ranges) U{rec}.meta.ranges(end)],'xticklabel',[-1 0 1])
            
        end
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
    
    if strcmp(U{1}.meta.layer,'BV')
        BV{rec}.peakIdx = peakidx;
    elseif strcmp(U{1}.meta.layer,'SM')
        SM{rec}.peakIdx = peakidx;
    elseif strcmp(U{1}.meta.layer,'D')
        D{rec}.peakIdx = peakidx;
    end
    %% test a trial to make sure that theta is aligned right
    %     trial = 1; %shifted by 1 (ex. trial=0 ..> trial =1)
    %     figure(trial+1);clf;plot(theta(:,trial+1));
    %     xlabel('Time from Trial Start (ms)');ylabel('Whisker Position')
    %
    %     validx=peakidx(floor(peakidx/U{rec}.t)==trial);
    %     excurx=round(((validx/U{rec}.t)-trial)*U{rec}.t);
    %     for i = 1:length(excurx)
    %         hold on; scatter(excurx(i),theta(excurx(i),trial+1),'ro')
    %     end
    %     %     %
    %     set(gca,'xlim',[550 1300])
    %     gos=U{1}.meta.trialType==1;
    %     nogos = U{1}.meta.trialType==0;
    
    %% max excursion
    if strcmp(array.meta.layer,'D')
        dbtheta = mean([gothetamean nogothetamean]); %theta = angle in b/t farthest go and closest nogo
%         dbtheta = max(gothetas(goprot));
        tidx = find(array.meta.trialType==0);
    else
        mp=U{rec}.meta.motorPosition;
        thetareqDB=polyval(E{rec}.coeff,mp,[],E{rec}.mu);
        
        %changing DB based on whether SM or BV
        if strcmp(array.meta.layer,'SM') || strcmp(array.meta.layer,'N')
            dbthetatmp = max(thetareqDB(U{rec}.meta.trialType==1));
            dbnogotheta = min(thetareqDB(U{rec}.meta.trialType==0));
            dbtheta = mean([dbthetatmp dbnogotheta]);
%                 dbtheta = max(thetareqDB(U{rec}.meta.trialType==1));
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
    %% Plotting Pole Availability Locations
    
    if strcmp(U{rec}.meta.layer,'BV')
        
        goMotor = [max(U{rec}.meta.ranges)-50000 max(U{rec}.meta.ranges)];
        goTheta = dbtheta - (polyval(E{rec}.coeff,goMotor,[],E{rec}.mu));
        nogoMotor = [min(U{rec}.meta.ranges) min(U{rec}.meta.ranges)+50000];
        nogoTheta = dbtheta - (polyval(E{rec}.coeff,nogoMotor,[],E{rec}.mu));
    elseif strcmp(U{rec}.meta.layer,'SM') || strcmp(U{rec}.meta.layer,'N')
        motorP = U{rec}.meta.motorPosition;
        go = U{rec}.meta.trialType == 1;
        nogo = ~go;
        goMotor = round([min(motorP(go)) max(motorP(go))],-4);
        goTheta = dbtheta - (polyval(E{rec}.coeff,goMotor,[],E{rec}.mu));
        nogoMotor = round([min(motorP(nogo)) max(motorP(nogo))],-4);
        nogoTheta = dbtheta - (polyval(E{rec}.coeff,nogoMotor,[],E{rec}.mu));
    elseif strcmp(U{rec}.meta.layer,'D')
        goTheta = dbtheta - [gothetamean gothetamean];
        nogoTheta = dbtheta - [nogothetamean nogothetamean];
    end
    
    popGo = [popGo ; goTheta];
    popNogo = [popNogo ; nogoTheta];
    
    
    %% Plotting Features for Search
    ranges=[round(min(xplor)-2):round(max(xplor)+2)];
    vals=histc(xplor,ranges);
    figure(17);subplot(2,5,rec)
    bar(ranges,vals/(sum(vals)),'k');
    hold on; plot([0 0],[0 1],'-.k','LineWidth',1)
    hold on; scatter(distmed(rec),.07,20,'rx','linewidth',2)
    hold on; plot([goTheta(1)-.5 goTheta(2)+.5],[.08 .08],'b','linewidth',5)
    hold on; plot([nogoTheta(1)-.5 nogoTheta(2)+.5],[.08 .08],'r','linewidth',5)
    
    set(gca,'xlim',[-50 50],'ylim',[0 .1]);
    alpha(.8)
    if rec == 3 || rec == 8
        xlabel('Theta from Discrimination Boundary');
    elseif rec == 1 || rec == 6
        ylabel('Proportion of Trial')
    end
    
end
% set(figure(70), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
% print(figure(70),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_thetaXmotor' ])
% set(figure(80), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
% print(figure(80),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_follicleXmotor' ])
% set(figure(17), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
% print(figure(17),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_searchDistribution'])

%% PLOTTING MEDIANSx

figure(376);

if strcmp(U{rec}.meta.layer,'N')
    hold on;scatter(ones(length(U),1).*0,distmed,40,'k')
    hold on;scatter(0,mean(distmed),100,'rx','linewidth',2)
    camroll(-90)
    hold on; plot([.5 3.5],[0 0], '-.k')
    set(gca,'ylim',[-20 20],'yaxislocation','right','xtick',[],'xlim',[.5 3.5])
    ylabel('Search Bias')
    
elseif strcmp(U{rec}.meta.layer,'D')
    hold on;scatter(ones(length(U),1).*1,distmed,'k')
    hold on;scatter(1,mean(distmed),100,'rx','linewidth',2)

elseif strcmp(U{rec}.meta.layer,'SM')
    hold on;scatter(ones(length(U),1).*2,distmed,'k')
    hold on;scatter(2,mean(distmed),100,'rx','linewidth',2)

elseif strcmp(U{rec}.meta.layer,'BV')
    hold on;scatter(ones(length(U),1).*3,distmed,'k')
    hold on;scatter(3,mean(distmed),100,'rx','linewidth',2)

end



% [~,newIdx] = sort(abs(distmed));
% for k = 1:length(U) %make new J array with U sorted by lowest median of search
%     J{k} = U{newIdx(k)};
% end


%% PLOTTING UBER HISTO
ranges = [-50:50];
ghisto= histc(cell2mat(gxplor'),ranges);
figure(58);
bar(ranges,ghisto/(sum(ghisto)),'k');
meanPopGo = mean(popGo);
meanPopNogo = mean(popNogo);
hold on; plot([meanPopGo(1)-.5 meanPopGo(2)+.5],[.045 .045],'b','linewidth',5)
hold on; plot([meanPopNogo(1)-.5 meanPopNogo(2)+.5],[.045 .045],'r','linewidth',5)

hold on; plot([0 0],[0 1],'-.k','LineWidth',1)
hold on;scatter(mean(distmed),.04,100,'rx','linewidth',2)
xlabel('Theta from Discrimination Boundary');ylabel('Proportion of Trial')
set(gca,'xlim',[-50 50],'ylim',[0 .05],'ytick',linspace(0,.05,3),'yticklabel',linspace(0,.05,3));
alpha(.8)
nanmedian(cell2mat(gxplor'));
nanstd(cell2mat(gxplor'));
print(figure(58),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_POPsearchDistribution'])

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

