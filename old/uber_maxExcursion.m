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

% type = {D,SM,BV};
type = {BV};
clearvars -except U D BV SM N type
close all
for p=1:length(type)
    close all
    U = type{p};
    %% Parameters
    % selectedMotorT filter to choose which type of trials to use to look at
    % max excursion. Options are:
    % 'nogo10k' 'allgo' 'allnogo' 'DB10k' 'FA' 'CR'
    
    selectedMotorT = 'nogo10k'  ;
    distmed=[];
    popGo = [];
    popNogo = [];
    mmdiff = zeros(1,length(U));
    totalFol_length = zeros(1,length(U));
    xyMaxChange = zeros(length(U),2);
    thetadiff = zeros(1,length(U));
    
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
            motorFollicleMean = [motorFollicleMean(:,1) motorFollicleMean(:,2:3)./33]; %conversion from pixels to mm  
            motorFollicleMean = [motorFollicleMean(:,1)    abs(motorFollicleMean(:,2)-max(motorFollicleMean(:,2)))      motorFollicleMean(:,3)-min(motorFollicleMean(:,3))]; %placing top right view of video as 0,0 
            
            
            
            %         if rec ==2 %error in SM2
            %             motorthetamean(1,:)=[];
            %         end
            %check to make sure there are at least some touches in nogo
            %positions
            minGo = min(U{rec}.meta.motorPosition(U{rec}.meta.trialType==1));
            nogoMotors = find(motorthetamean(:,1)<minGo);
            
            if ~isempty(nogoMotors)
                
                %THETAxMOTOR PLOTS 
                
                figure(70); subplot(2,5,rec);
                scatter(motorthetamean(:,1),motorthetamean(:,2),'k');
                hold on; plot(motorthetamean(:,1),f,'k')
                xlabel('Normalized Motor Positions'); ylabel('Theta at base')
                [E{rec}.coeff, ~ , E{rec}.mu] = polyfit(motorthetamean(:,1),motorthetamean(:,2),2);
                f = polyval(E{rec}.coeff,motorthetamean(:,1),[],E{rec}.mu);
               
                mp=U{rec}.meta.motorPosition;
                thetareqDB=polyval(E{rec}.coeff,mp,[],E{rec}.mu);
                dbtheta = max(thetareqDB(U{rec}.meta.trialType==1));
                hold on; plot([0 180000],[dbtheta dbtheta],'-.k')
                %             set(gca,'xlim',[U{rec}.meta.ranges],'xtick',[U{rec}.meta.ranges(1) mean(U{rec}.meta.ranges) U{rec}.meta.ranges(2)],'xticklabel',[-1 0 1])
                set(gca,'xtick',[U{rec}.meta.ranges(1) mean(U{rec}.meta.ranges) U{rec}.meta.ranges(end)],'xticklabel',[-1 0 1])
                
                %this bit here used to calculate the theta diff between
                %DB+/-.5mm
%                 tmp=[find(motorthetamean(:,1)>=mean(U{rec}.meta.ranges)+5000,1) find(motorthetamean(:,1)>=mean(U{rec}.meta.ranges)-5000,1)];
%                 thetadiff(rec)=diff(f(tmp));
                
                %FOLLICLExMOTOR PLOTS
               
                figure(80);subplot(2,5,rec)
                colormap(fliplr(redbluecmap))
                goFols = motorFollicleMean(:,1)>=mean(U{rec}.meta.ranges);
                nogoFols = motorFollicleMean(:,1)<mean(U{rec}.meta.ranges);
%                 scatter(motorFollicleMean(:,2),motorFollicleMean(:,3),[],motorFollicleMean(:,1)) %plot xy follicle coordinates, colors = motorPosition
                scatter(motorFollicleMean(goFols,2),motorFollicleMean(goFols,3),'b','filled') 
                hold on; scatter(motorFollicleMean(nogoFols,2),motorFollicleMean(nogoFols,3),'r','filled') 
                xlabel('AP Translation (mm)');ylabel('ML Translation(mm)');
                set(gca, 'YDir','reverse','XDir','reverse','ytick',0:1:2,'ylim',[0 2],'xlim',[0 2],'xtick',0:1:2)
                legend('Go','NoGo')
                
                
               
                if strcmp(U{1}.meta.layer,'BV')
                     %find angle difference between 1mmgo and 1mmnogo
                    gomm=find(motorthetamean(:,1)>=mean(U{rec}.meta.ranges)+10000,1);
                    nogomm=find(motorthetamean(:,1)<=mean(U{rec}.meta.ranges)-10000);
                    nogomm = nogomm(end);
                    mmdiff(rec) = motorthetamean(nogomm,2)-motorthetamean(gomm,2);
                    
                    %find follicle travel and xy change across distances
                    xyFollicle = sortrows(motorFollicleMean(:,2:3));
                    [FolCoeff, ~ , FolMu] = polyfit(xyFollicle(:,1),xyFollicle(:,2),2);
                    g = polyval(FolCoeff,xyFollicle(:,1),[],FolMu);
                    d = diff([xyFollicle(:,1) g(:)]);
                    totalFol_length(rec) = sum(sqrt(sum(d.*d,2)));
                    xyMaxChange(rec,:) = max(xyFollicle);
                end
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
    
    %% PLOTTING UBER HISTO
    ranges = [-50:50];
    ghisto= histc(cell2mat(gxplor'),ranges);
%     ghisto= histc(ganovacomparison{1},ranges);
    figure(58+p);
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
    
    % print(figure(58),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_POPsearchDistribution'])
    % set(figure(70), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
    % print(figure(70),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_thetaXmotor' ])
%     set(figure(80), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
%     print(figure(80),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_follicleXmotor' ])
%     set(figure(17), 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
%     print(figure(17),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_searchDistribution'])
    ganovacomparison{p} = cell2mat(gxplor') ;
end
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


%% ANOVA FOR DISTRIBUTION COMPARISON
% alpha = .01; %significance level of .01
% [p,tbl,stats] = anovan(group); %anova1
% comp = multcompare(stats); %comparison between all groups.
% bonfcorr = alpha/10; %post hoc bonferroni correction of pval alpha/n
%
% [comp(:,1:2) comp(:,end)<bonfcorr ]

% OR USE KSTEST2
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
figure(56);clf
for d = 1:length(colors)
    [f,x_values] = ecdf(ganovacomparison{d});
    hold on;
    F = plot(x_values,f,'color',rgb(colors{d}));
    set(F,'LineWidth',2);
end
xlabel('Theta from Discrimination Boundary');ylabel('Cumulative Distribution');
% title('Exploration Strategy Empirical CDF')
set(gca,'ytick',[0 .5 1],'xtick',[-50 0 50])
set(figure(56), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);
legend('Discrete','Semi-Continuous','Continuous','location','northwest')

[h,p] = kstest2(ganovacomparison{1},ganovacomparison{2});
[h2,p2] = kstest2(ganovacomparison{1},ganovacomparison{3});
[h3,p3] = kstest2(ganovacomparison{2},ganovacomparison{3});
ksnums.names = {'g1', 'g2' ,'reject?', 'pval'};
ksnums.vals = [1 2 h p; 1 3 h2 p2 ; 2 3 h3 p3];


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

