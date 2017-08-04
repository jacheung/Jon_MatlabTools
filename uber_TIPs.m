avgtouch=zeros(1,length(U));
avgpc=zeros(1,length(U));
avgduration=zeros(1,length(U));
Gvtnorm=cell(1,length(U));

for rec=1:length(U) %for each cell
    
    [touchespredecision,GprelixtimesGO,GprelixtimesNOGO,Gprelixpcth,Gprelixduration,Gprelixtimes] = assist_predecisionVar(U,rec);

    %plotting the peri-contact time histogram of predecision touches(calculated as diff between
    %previous touch). Only viewing touches with PCTH of <500ms
    totalpcth=sort(cell2mat(Gprelixpcth));
    totalpcth500=totalpcth(totalpcth<500);
    figure(25);subplot(2,1,1);hist(totalpcth500,50);
    title('Time between pre-decision touches')
    %xlabel('time(ms)')
    ylabel('# touches')
    text(max(totalpcth500)/2, max(hist(totalpcth500,50))*.9,['Average # of touches predecision:' num2str(touchespredecision)],'FontSize',8,'Color','black')
    
    %plotting the duration in histogram (calculated from onset:offset) of all
    %predecision touches.
    totalduration=sort(cell2mat(Gprelixduration));
    totalduration500=totalduration(totalduration<500);
    subplot(2,1,2);hist(totalduration500,50);
    title('Duration of pre-decision touches')
    xlabel('time(ms)');ylabel('# touches')
    text(max(totalduration500)/2, max(hist(totalduration500,50))*.9,['Average # of touches predecision:' num2str(touchespredecision)],'FontSize',8,'Color','black')
   
    pc=mean(U{rec}.meta.trialCorrect);
    avgtouch(1,rec)=touchespredecision;
    avgpc(1,rec)=pc;
    avgduration(1,rec)=mean(totalduration);
    
    Gprelix=Gprelixtimes;
    %% ASIDE: plotting 1 cell, # of touches vs percent correct 
    %This section can be uprooted and replaced by uber_motorposVSaccuracy
    
%     tmpGO=[cellfun(@numel,GprelixtimesGO);U{rec}.meta.trialCorrect(U{rec}.meta.trialType)]';
%     tmpNOGO=[cellfun(@numel,GprelixtimesNOGO);U{rec}.meta.trialCorrect(logical(1-U{rec}.meta.trialType))]';
%     tmp2GO=sortrows(tmpGO,1);tmp2GO(:,3)=tmp2GO(:,2).*tmp2GO(:,1);
%     tmp2NOGO=sortrows(tmpNOGO,1);tmp2NOGO(:,3)=tmp2NOGO(:,2).*tmp2NOGO(:,1);
%     % figure;plot(1:length(tmp2),tmp2(:,1)); hold on; scatter(1:length(tmp2),tmp2(:,3))
%     
%     [sorted sortedBy binBounds]=binslin(tmp2GO(:,1),tmp2GO(:,1:2),'equalE',14,0,12);
%     touchpcGO=zeros(length(sorted),2);
%     for i=1:length(sorted)
%         touchpcGO(i,:)=mean(sorted{i},1);
%     end
%     [sortedNO sortedByNO binBoundsNO]=binslin(tmp2NOGO(:,1),tmp2NOGO(:,1:2),'equalE',14,0,12);
%     touchpcNO=zeros(length(sortedNO),2);
%     for i=1:length(sortedNO)
%         touchpcNO(i,:)=mean(sortedNO{i},1);
%     end
%     figure(28);plot(touchpcGO(:,1),touchpcGO(:,2))
%     hold on; plot(touchpcNO(:,1),touchpcNO(:,2),'r');
%     xlabel('# of touches'),ylabel('% trials correct');
%     set(gca,'ylim',[0 1.2])
    
    %% ASIDE: plotting variable of 1st touch vs subsequent touches PRE DECISION
    %2= velocity 3 = amp, 4 = setpoint, 5=phase, 6=curvature
   
    varnames={'vel','amp','setpoint','phase'};
    for var = 2:5;
        vtouch=cell(1,numel(Gprelix));
        vtnorm=zeros(numel(Gprelix),20);vtnorm(:)=NaN;
        
        for j=1:length(Gprelix)
            if var == 2;
                for m=1:numel(Gprelix{j})
                    vtouch{j}(m)=mean(U{rec}.S_ctk(var,Gprelix{j}(m)-5:Gprelix{j}(m),j));
                end
            else
                vtouch{j}=U{rec}.S_ctk(var,Gprelix{j},j); %find associated var at touch
            end
        end
        
        for k = 1:length(Gprelix) %normalize across touches within that trial
            if numel(vtouch{k})>1
                vtnorm(k,1:numel(vtouch{k}))=(vtouch{k}-min(vtouch{k}))/(max(vtouch{k})-min(vtouch{k}));
            end
            if numel(vtouch{k})==1
                vtnorm(k)=1;
            end;
        end
        figure(41);subplot(2,2,var-1);
        %for k =1:length(vtnorm)
            %scatter(1:sum(~isnan(vtnorm(k,:))),vtnorm(k,1:sum(~isnan(vtnorm(k,:)))),'k'); hold on
        %end
        
    plot(nanmean(vtnorm,1),'k');xlim([0 15]);ylim([0 1])
    Gvtnorm{rec}=vtnorm;
    ylabel(['Normalized' varnames(var-1) 'within Trial']);xlabel('Touch Number');
    end
    
    
    %% ASIDE: plotting setpoint for hits,CR,FA,misses NEED TO FIX THIS
    hits=vtouch(find(U{rec}.meta.trialType==1 & U{rec}.meta.trialCorrect==1));
    CR=vtouch(find(U{rec}.meta.trialType==0 & U{rec}.meta.trialCorrect==1));
    FA=vtouch(find(U{rec}.meta.trialType==0 & U{rec}.meta.trialCorrect==0));
    miss=vtouch(find(U{rec}.meta.trialType==1 & U{rec}.meta.trialCorrect==0));
    
    figure(44);clf
    
    [sorted sortedBy binBounds]=binslin(cell2mat(hits),cell2mat(hits),'equalE',50,-50,50);
    
    hist(cell2mat(hits));hist(cell2mat(FA),'r');
    histfit(cell2mat(hits));hold on;histfit(cell2mat(CR));histfit(cell2mat(FA));histfit(cell2mat(miss))
    
    
    %% ASIDE: Plotting whisker occupancy 
    GprelixtimesGO
  
end

%want to find the mean of the distribution but can't decide which prob dist
%to use.
% fit=fitdist(totalpcth500','lognormal');
% fitdf=pdf(fit,1:max(totalpcth500));
% hold on; plot(1:max(totalpcth500),fitdf*1000);

%% Population level for L5b cells, # of touches correlated with percent correct
%look at population level and see if there is correlation between % correct and duration of predecision spikes

figure(27);scatter(avgtouch,avgpc);lsline
xlabel('Mean # of touches pre lick'), ylabel('Percent Trials (w/ predecision touches) Correct');

[rho,pval]=corr([avgtouch; avgpc],'type','pearson');

figure(28);scatter(avgduration,avgpc);lsline
xlabel('Mean duration of touches pre lick'), ylabel('Percent Trials (w/ predecision touches) Correct');

figure(29);scatter3(avgduration,avgtouch,avgpc);
xlabel('Average Touch Duration');ylabel('Average # of Touches'),zlabel('Percent Trials (w/ predecision touches) Correct')
% on the population level, plotting time between pre decision touches and
% percent correct.




