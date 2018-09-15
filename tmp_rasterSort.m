if strcmp(U{1}.meta.layer,'L5b')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','I23:J50');
    disp(U{1}.meta.layer);
    layer = 'L5b';
elseif strcmp(U{1}.meta.layer,'L3')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','B25:C44');
    disp(U{1}.meta.layer);
    layer = 'L3';
elseif strcmp(U{1}.meta.layer,'L4')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','O23:P28');
    disp(U{1}.meta.layer);
    layer = 'L4';
elseif strcmp(U{1}.meta.layer,'L3Out')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','B77:C81');
    disp(U{1}.meta.layer);
    layer = 'L3Out';
else strcmp(U{1}.meta.layer,'L5bOut')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','I75:J80');
    disp(U{1}.meta.layer);
    layer = 'L5bOut';
end
txt=txt(~isnan(ndata));
ndata=ndata(~isnan(ndata));

 
for p=1:length(ndata)
    cellcode = txt{p};
    rec = p;
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
    title('NoGo')
    
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
%     print(gcf,'-dpng',['Z:\Users\Jon\Projects\Characterization\' layer '\Figures\' cellcode '_' 'Behavioral_Sorted_Raster'])
end
    %% plot raster go v no go v CR v FA v HIT v Miss
    figure(41);clf
    subplot(6,2,[1 2]);hold on
    for i = 1:U{rec}.k
        if sum(U{rec}.R_ntk(1,:,i))>0 & U{rec}.meta.trialType(i)==1 & U{rec}.meta.trialCorrect(i)==1%plot HIT
            plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
        else
        end 
    end
    %set(gca,'ylim',[0 U{rec}.k]+1,'xticklabel',[])
    set(gca,'xticklabel',[])
    title('Hit')
    
    
    subplot(6,2,[3 4]);hold on
    for i = 1:U{rec}.k
        if sum(U{rec}.R_ntk(1,:,i))>0 & U{rec}.meta.trialType(i)==1 & U{rec}.meta.trialCorrect(i)==0%plot MISS
            plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
        else
        end 
    end
    set(gca,'xticklabel',[])
    title('Miss')
    
    subplot(6,2,[5 6]);hold on
    for i = 1:U{rec}.k
        if sum(U{rec}.R_ntk(1,:,i))>0 & U{rec}.meta.trialType(i)==0 & U{rec}.meta.trialCorrect(i)==1%plot CR
            plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
        else
        end 
    end
    set(gca,'xticklabel',[])
    title('Correct Rejection')
     
    subplot(6,2,[7 8]);hold on
    for i = 1:U{rec}.k
        if sum(U{rec}.R_ntk(1,:,i))>0 & U{rec}.meta.trialType(i)==0 & U{rec}.meta.trialCorrect(i)==0%plot FA
            plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
        else
        end 
    end
    title('False Alarm')
    
    goidx=find(U{rec}.meta.trialType==1);
    nogoidx=find(U{rec}.meta.trialType==0);
    coridx=find(U{rec}.meta.trialCorrect==1);
    erroridx=find(U{rec}.meta.trialCorrect==0);
    
    hitidx=intersect(goidx,coridx);
    missidx=intersect(goidx,erroridx);
    cridx=intersect(nogoidx,coridx);
    faidx=intersect(nogoidx,erroridx);
    
    hitspks=zeros(length(hitidx),4000);
    missspks=zeros(length(missidx),4000);
    crspks=zeros(length(cridx),4000);
    faspks=zeros(length(faidx),4000);
    swindow=50;
    
    for i=1:length(hitidx)
        hitspks(i,:)=U{rec}.R_ntk(:,:,hitidx(i));
    end
    
    for i=1:length(missidx)
        missspks(i,:)=U{rec}.R_ntk(:,:,missidx(i));
    end
    
    for i=1:length(cridx)
        crspks(i,:)=U{rec}.R_ntk(:,:,cridx(i));
    end
    
    for i=1:length(faidx)
        faspks(i,:)=U{rec}.R_ntk(:,:,faidx(i));
    end
    
    avghitspks=smooth(mean(hitspks)*1000,swindow);
    avgmissspks=smooth(mean(missspks)*1000,swindow);
    avgcrspks=smooth(mean(crspks)*1000,swindow);
    avgfaspks=smooth(mean(faspks)*1000,swindow);
    
    subplot(6,2,[9 12]);hold on
    plot(1:length(gospks),avghitspks,'g'); hold on
    plot(1:length(gospks),avgmissspks,'r'); hold on
    plot(1:length(gospks),avgcrspks,'b'); hold on
    plot(1:length(gospks),avgfaspks,'c'); hold on
    plot([avail avail],[0 max(avggospks)],'k:')
    plot([onset onset],[0 max(avggospks)],'k:')
    set(gca,'xtick',[])
    ylabel('spks/s')
    xlabel('time from trial start')
end