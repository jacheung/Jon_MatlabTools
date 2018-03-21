clear all;close all;

mouselist = {'AH0667','AH0668','AH0712','AH0716','AH0717'};

for mouse = 1:length(mouselist)
    close all
    mouseNum = mouselist{mouse};
    
    d = (['Z:\Users\Jon\DATA\Behavior\OnAxes']);
    cd(d);
    filelist=dir([d filesep '*.mat']);
    list = [];
    for f = 1:length(filelist)
        if isequal(filelist(f).name(33:38),mouseNum) == 1
            list = [list f];
        end
    end
    
    gpsycho = [];
    B{mouse}.name = mouseNum;
    B{mouse}.acc = zeros(length(list),3);
    
    
    for k = 1:length(list)
        load([ filelist(list(k)).name]);
        
        %%
        hor=saved.MotorsSection_previous_pole_positions;
        lat=saved.MotorsSection_previous_pole_positions_lat;
        arc=saved.SidesSection_previous_trial_types==97;
        rad=saved.SidesSection_previous_trial_types==114;
        norm=saved.SidesSection_previous_trial_types==110;
        go = saved.SidesSection_previous_sides==114;
        nogo = saved.SidesSection_previous_sides==108;
        
        %trial correct
        hist = saved_history;
        tmp = hist.AnalysisSection_PercentCorrect;
        tmp2=[0;tmp(1:end-1)];
        rewtmp=cell2mat(tmp)-cell2mat(tmp2);
        tcorr=(rewtmp>=0)';
        lastT=find((go(1:end-1).*tcorr),1,'last');
        %lastT = 520 %for AH0667 170718
        
        tcorr = tcorr(1:lastT);
        atrim = arc(1:lastT);
        ntrim = norm(1:lastT);
        rtrim = rad(1:lastT);
        gotrim = go(1:lastT);
        nogotrim = nogo(1:lastT);
        %% Sorting Trials based on Type
        atrials=sum(atrim);
        ahit=sum(tcorr.*atrim.*gotrim);
        aCR=sum(tcorr.*atrim.*nogotrim);
        aFA=sum(~tcorr.*atrim.*nogotrim);
        amiss=sum(~tcorr.*atrim.*gotrim);
        a = [ahit aCR aFA amiss]./atrials;
        
        aacc = (ahit +aCR)/atrials;
        
        ntrials=sum(ntrim);
        nhit=sum(tcorr.*ntrim.*gotrim);
        nCR=sum(tcorr.*ntrim.*nogotrim);
        nFA=sum(~tcorr.*ntrim.*nogotrim);
        nmiss=sum(~tcorr.*ntrim.*gotrim);
        n = [nhit nCR nFA nmiss]./ntrials;
        nacc = (nhit +nCR)/ntrials;
        
        rtrials=sum(rtrim);
        rhit=sum(tcorr.*rtrim.*gotrim);
        rCR=sum(tcorr.*rtrim.*nogotrim);
        rFA=sum(~tcorr.*rtrim.*nogotrim);
        rmiss=sum(~tcorr.*rtrim.*gotrim);
        r=[rhit rCR rFA rmiss]./rtrials;
        racc = (rhit +rCR)/rtrials;
        
        B{mouse}.acc(k,:) = [aacc nacc racc];
        B{mouse}.name = {'arc','normal','radial'};
        %% psychometric curves
        psys{1} = [tcorr(atrim)' gotrim(find(atrim==1))' hor(atrim)']; %arc task
        psys{2} = [tcorr(rtrim)' gotrim(find(rtrim==1))' lat(rtrim)']; %radial task
        psys{3} = [tcorr(ntrim)' gotrim(find(ntrim==1))' hor(ntrim)']; %continuous
        psys{3} = psys{3}(2:end,:); %error in 1st trial sometimes setting everything to 0
        
        for j = 1:length(psys)
            agoT = find(psys{j}(:,2)==1);anogoT = find(psys{j}(:,2)==0);tmp=psys{j}(:,3);
            
            discrimB = round(mean([min(tmp(agoT)) max(tmp(anogoT))]),-3);
            goR = round(max(tmp(agoT)),-3); nogoR = round(min(tmp(anogoT)),-3);
            [gosorted, sortedBy, binBounds] = binslin(psys{j}(:,3),psys{j},'equalE',6,discrimB,goR);
            [nogosorted, sortedBy, binBounds] = binslin(psys{j}(:,3),psys{j},'equalE',6,nogoR,discrimB);
            means=cell2mat(cellfun(@mean,[nogosorted ;gosorted],'uniformoutput',0));
            psycho=[1-means(1:5,1) ;means(6:10,1)];
            gpsycho = [gpsycho psycho]
            
            
            if j == 3
                figure(54);hold on; plot(1:10,psycho,'k','linewidth',.1);
                set(gca,'xtick',[1 5.5 10],'xticklabel',[-1 0 1],'ytick',[0 .25 .5 .75 1],'ylim', [ 0 1])
            else
                figure(53);subplot(2,1,j);hold on; plot(1:10,psycho,'k','linewidth',.1);
                set(gca,'xtick',[1 5.5 10],'xticklabel',[-1 0 1],'ytick',[0 .25 .5 .75 1],'ylim', [ 0 1])
                if j == 1
                    title('Arc (TOP) , Radial (BOTTOM) Psychometric Curves')
                end
                
            end
        end
        ylabel('P(lick)');xlabel('Discrimination Boundary');
        
        %% Are Mice Improving in the Task Overtime?
        avec=tcorr(atrim);
        rvec=tcorr(rtrim);
        nvec=tcorr(ntrim);
        figure(8);clf;
        plot(1:length(avec),smooth(avec,25),'b')
        hold on; plot(1:length(rvec),smooth(rvec,25),'c')
        set(gca,'xlim',[25 length(rvec)-25],'ylim',[.25 1])
        xlabel('Number of Trials'); ylabel('% Accuracy')
        legend('ARC','RADIAL')
        %hold on; plot(1:length(nvec),smooth(nvec,20),'m')
        
        %% Basic Bar Charts
        types= [n;a;r];
        
        myC=[0 0 1; 1 0 0 ; 0 1 0 ; 0 0 0]; % make a colors list
        figure(82);clf;H=bar(types,'stacked');
        for k=1:4
            set(H(k),'facecolor',myC(k,:))
        end
        ylabel('Proportion of Trials')
        set(gca,'xticklabel',{['norm (n = ' num2str(ntrials) ')' ] ;...
            ['arc (n = ' num2str(atrials) ')']  ;['radial (n = ' num2str(rtrials) ')']})
        
        figure(84);clf;scatter(hor(ntrim)/10000,lat(ntrim)/10000,'b')
        hold on; scatter(hor(atrim)/10000,lat(atrim)/10000,'bx');
        hold on; scatter(hor(rtrim)/10000,lat(rtrim)/10000,'b+');
        xlabel('Horizontal Motor Position (mm)')
        ylabel('Lateral Motor Position (mm)')
        
        figure(86);clf;scatter(hor(go(1:end-1)),lat(go(1:end-1)),'b')
        hold on;scatter(hor(nogo(1:end-1)),lat(nogo(1:end-1)),'r')
        xlabel('Horizontal Motor Position (mm)')
        ylabel('Lateral Motor Position (mm)')
        
        
    end
    figure(53);
    subplot(2,1,1);plot([1:10],mean(gpsycho(:,[1:3:end]),2),'r','linewidth',4);
    subplot(2,1,2);plot([1:10],mean(gpsycho(:,[2:3:end]),2),'r','linewidth',4);
    figure(54);
    title ('Continuous Psychometric Curves')
    plot([1:10],mean(gpsycho(:,[3:3:end]),2),'r','linewidth',4);
    
    figure(45);
    plot([1:10],mean(gpsycho(:,[3:3:end]),2),'k','linewidth',4); %continuous
    hold on; plot([1:10],mean(gpsycho(:,[1:3:end]),2),'b','linewidth',4); %arc
    plot([1:10],mean(gpsycho(:,[2:3:end]),2),'c','linewidth',4);%radial
    legend('Horizontal','Arc','Radial','location','southeast');
    set(gca,'xtick',[1 5.5 10],'xticklabel',[-1 0 1],'ytick',[0 .25 .5 .75 1],'ylim', [ 0 1])
    xlabel('Normalized Motor Positions');ylabel('Lick Probability')
    %     print(figure(45),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\' mouseNum '_onAxesPsycho'])
    
    
    mouseTotal{mouse} = [mean(gpsycho(:,[3:3:end]),2); mean(gpsycho(:,[1:3:end]),2); mean(gpsycho(:,[2:3:end]),2)];
    
    ALLgpsycho{mouse} = gpsycho;
end

popAVG = [];
figure(10);clf
for d = 1:length(B)
    hold on; errorbar(1:3,mean(B{d}.acc),std(B{d}.acc),'-o','color',[.5 .5 .5])
    popAVG = [popAVG ; B{d}.acc];
end
errorbar(1:3,mean(popAVG),std(popAVG),'k-o','linewidth',4)
set(gca,'xlim',[.5 3.5],'xtick',[1 2 3],'xticklabel',{'Arc','Horizontal','Radial'},...
    'ylim',[0 1],'ytick',[0:.25:1],'yticklabel',[0:25:100])
xtickangle(45)
ylabel('Percent Accuracy')

% set(gcf, 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
print(figure(10),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\' mouseNum '_onAxes_POPaccuracy'])

alpha = .01; %significance level of .01 
[p,tbl,stats] = anova1(popAVG); %anova1 
comp = multcompare(stats); %comparison between all groups. 
bonfcorr = alpha/numel(mouselist); %post hoc bonferroni correction of pval alpha/n

[comp(:,1:2) comp(:,end)<bonfcorr ]

% Plotting group psychometric curves with ALL shaded psychos
popAll = mean(cell2mat(mouseTotal),2);
semAll = std(cell2mat(mouseTotal),0,2)./sqrt(numel(mouselist));
tn = tinv([0.025 .975],numel(mouselist)-1);
% ci = repmat(popAll,1,2)+(tn.*semAll);
ci = (tn.*semAll);

contAll = popAll(1:10); contci = ci(1:10,:); contSEM = semAll(1:10);
arcAll = popAll(11:20); arcci = ci(11:20,:); arcSEM = semAll(11:20);
radAll = popAll(21:30); radci = ci(21:30,:); radSEM = semAll(21:30);


    
    figure(438);clf
    errorbar([1:10],radAll,radSEM,'-c','linewidth',2)
    hold on; errorbar([1:10],contAll,contSEM,'-k','linewidth',2)
    errorbar([1:10],arcAll,arcSEM,'-b','linewidth',2)
    
    figure(439);clf
    boundedline([1:10],radAll,radci(:,2),'c');
    hold on; boundedline([1:10],contAll,contci(:,2),'k')
    boundedline([1:10],arcAll,arcci(:,2),'b')

    
% for i = 1:length(ALLgpsycho)
%    figure(438);
%     h=plot([1:10],ALLgpsycho{i}(:,[3:3:end]),'color',[0 0 0]+alpha,'linewidth',1);  %continuous
%         for d = 1:length(h) 
%             h(d).Color(4)=.2;
%         end
%     hold on; h=plot([1:10],ALLgpsycho{i}(:,[1:3:end]),'color',[0 0 1],'linewidth',1);  %arc
%     for d = 1:length(h) 
%             h(d).Color(4)=.2;
%         end
%     h=plot([1:10],ALLgpsycho{i}(:,[2:3:end]),'color',[0 1 1],'linewidth',1); %radial
%     for d = 1:length(h) 
%             h(d).Color(4)=.2;
%         end
% end

%     plot([1:10],contAll,'k','linewidth',4); %continuous
%     hold on; plot([1:10],arcAll,'b','linewidth',4); %arc
%     plot([1:10],radAll,'c','linewidth',4);%radial
    legend('Radial','Horizontal','Arc','location','southeast');
    set(gca,'xtick',[1 5.5 10],'xticklabel',[-1 0 1],'ytick',[0 .25 .5 .75 1],'ylim', [ 0 1])
    xlabel('Normalized Motor Positions');ylabel('Lick Probability')
    text(8.5,.20,['n = ' num2str(numel(mouselist))],'FontSize',14)
    print(figure(438),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\onAxes_POPpsycho'])
    
    
    
    
