clear all;close all; 

mouseNum = 'AH0668';
d = (['Z:\Users\Jon\DATA\Behavior\OnAxes']);
cd(d);
filelist=dir([d filesep '*.mat']);
list = [];
for f = 1:length(filelist)
    if isequal(filelist(f).name(33:38),mouseNum) == 1
        list = [list f];
    end
end
samp=4
gpsycho = []
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
    
    ntrials=sum(ntrim);
    nhit=sum(tcorr.*ntrim.*gotrim);
    nCR=sum(tcorr.*ntrim.*nogotrim);
    nFA=sum(~tcorr.*ntrim.*nogotrim);
    nmiss=sum(~tcorr.*ntrim.*gotrim);
    n = [nhit nCR nFA nmiss]./ntrials;
    
    rtrials=sum(rtrim);
    rhit=sum(tcorr.*rtrim.*gotrim);
    rCR=sum(tcorr.*rtrim.*nogotrim);
    rFA=sum(~tcorr.*rtrim.*nogotrim);
    rmiss=sum(~tcorr.*rtrim.*gotrim);
    r=[rhit rCR rFA rmiss]./rtrials;
    
    B.rAcc(samp)=(rhit+rCR)/rtrials;
    B.nAcc(samp)=(nhit+nCR)/ntrials;
    B.aAcc(samp)=(ahit+aCR)/atrials;
    
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
    
    
    %% Plot all Psycho
    
    
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


%%
ravg = mean(B.rAcc);rstd = std(B.rAcc);
aavg = mean(B.aAcc);astd = std(B.aAcc);
navg = mean(B.nAcc);nstd = std(B.nAcc);

m = [navg aavg ravg];
s = [nstd astd rstd];

figure(10);errorbar([1:3],m,s,'o')
title('Population Average n=2 (2 sessions ea)');ylabel('% Accuracy')
set(gca,'xlim',[0 4],'ylim',[.5 1],'xtick',[1:3],'xticklabel',{'normal';'arc';'radial'},'ytick',...
    [.5:.1:1])