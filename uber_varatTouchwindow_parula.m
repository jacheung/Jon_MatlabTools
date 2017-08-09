if strcmp(U{1}.meta.layer,'L5b')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','I23:J65');
    disp(U{1}.meta.layer);
    layer = 'L5b';
elseif strcmp(U{1}.meta.layer,'L3')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','B25:C44');
    disp(U{1}.meta.layer);
    layer = 'L3';
elseif strcmp(U{1}.meta.layer,'L4')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','O23:P28');
    disp(U{1}.meta.layer);
    layer = 'L4';
elseif strcmp(U{1}.meta.layer,'L3Out')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','B77:C82');
    disp(U{1}.meta.layer);
    layer = 'L3Out';
elseif strcmp(U{1}.meta.layer,'L5bOut')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','I75:J82');
    disp(U{1}.meta.layer);
    layer = 'L5bOut';
else strcmp(U{1}.meta.layer,'L5bInt')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','V75:W81');
    disp(U{1}.meta.layer);
    layer = 'L5bInt';
end
txt=txt(~isnan(ndata));
thetaall=cell(1,length(U));
phaseall=cell(1,length(U));
ampall=cell(1,length(U));
spall=cell(1,length(U));

for p=34
    rec = p;
    normbinmin=10; %min in each bin to use for normalizing
    normbinwin=30; %window for normalizing from 0ms:binwin
    cellcode = txt{p};
    [varspikes, ~, ~] = assist_varAtTouch(U{rec},[-25:50]);
    % First 6 columns will be values for the variables
    % 1) THETA 2) AMP 3) SETPOINT 4) PHASE 5) MAX KAPPA 6) PRE TOUCH VELOCITY
    % Last columns will be the spikes around your given window
    %%
    %Plot theta at touch
    
    [sorted sortedBy binBounds]=binslin(varspikes(:,1),varspikes(:,7:82),'equalE',41,-100,100);
    thetarange=[-97.5:5:97.5];
    
    thetatouchspks=zeros(length(sorted),79);
    for j=1:length(sorted)
        thetatouchspks(j,2:77)=mean(sorted{j},1);
        thetatouchspks(j,78)=size(sorted{j},1);
    end
    thetatouchspks(:,1)=thetarange;
    %thetatouchspks = thetatouchspks(all(~isnan(thetatouchspks),2),:); %remove
    %any row with a NaN value in it.. not good so use one below this
    thetatouchspks = thetatouchspks(~any(isnan(thetatouchspks(:,2)),2),:);%remove NaN rows USING column for theta
    for k=1:size(thetatouchspks,1)
        thetatouchspks(k,2:77)=smooth(thetatouchspks(k,2:77));
        if thetatouchspks(k,78)>normbinmin
            thetatouchspks(k,79)=max(thetatouchspks(k,26:26+normbinwin));
        end
    end
    
    %thetatouchspks(8,79)=.05 %tmp for AH since he wanted .05spks/ms
    figure(30+g);clf;h1 = subplot(2,3,1);
    imagesc(thetatouchspks(:,2:77))
    hold on
    colormap(gca,parula)
    caxis([0 max(thetatouchspks(:,79))])
    plot([25 25],[25 0],'w:')
    set(gca,'Ydir','normal','ytick',(1:length(thetatouchspks(:,1))),'yticklabel',[thetatouchspks(:,1)],'xtick',(0:25:75),'xticklabel',(-25:25:50));
    for k=1:size(thetatouchspks,1)
        text(20,k,num2str(thetatouchspks(k,78)),'FontSize',8,'Color','white')
    end
    axis('square')
    xlabel('time from touch onset (ms)')
    title('Theta at touch')
    
    %% Plot amplitude
    [ampsorted ampsortedBy ampbinBounds]=binslin(varspikes(:,2),varspikes(:,7:82),'equalX',13);
    %amprange=[1:2:99];
    
    for y=1:length(ampbinBounds)-1
        amprange(y)=mean([ampbinBounds(y),ampbinBounds(y+1)]);
    end
    
    amptouchspks=zeros(length(ampsorted),79);
    for j=1:length(ampsorted)
        amptouchspks(j,2:77)=mean(ampsorted{j},1);
        %amptouchspks(j,1)=mean(ampsortedBy{j});
        amptouchspks(j,78)=size(ampsorted{j},1);
    end
    
    amptouchspks(:,1)=amprange;
    amptouchspks = amptouchspks(all(~isnan(amptouchspks),2),:); %remove NaN rows
    
    for k=1:size(amptouchspks,1)
        amptouchspks(k,2:77)=smooth(amptouchspks(k,2:77));
        if amptouchspks(k,78)>normbinmin
            amptouchspks(k,79)=max(amptouchspks(k,26:26+normbinwin));
        end
    end
    
    yticksamp = cell(length(amptouchspks(:,1)));
    for i = 1:length(yticksamp)
        yticksamp{i} = sprintf('%0.2g',amptouchspks(i,1));
    end
    
    subplot(2,3,2);
    imagesc(amptouchspks(:,2:77))
    hold on
    plot([25 25],[25 0],'w:')
    set(gca,'Ydir','normal','xtick',(0:25:75),'xticklabel',(-25:25:50),'ytick',[1:size(amptouchspks,1)],'yticklabel',yticksamp);
    %hCBar=colorbar;
    colormap(gca,parula)
    caxis([0 max(thetatouchspks(:,79))])
    axis('square')
    
    for k=1:size(amptouchspks,1)
        text(20,k+.1,num2str(amptouchspks(k,78)),'FontSize',8,'Color','white')
    end
    title('Amplitude at Touch')
    %% setpoint
    % [spsorted spsortedBy spbinBounds]=binslin(varspikes(:,3),varspikes(:,7:82),'equalE',101,-100,100);
    % sprange=[-99:2:99];
    [spsorted spsortedBy spbinBounds]=binslin(varspikes(:,3),varspikes(:,7:82),'equalX',13);
    
    for y =1:length(spbinBounds)-1
        sprange(y)=mean([spbinBounds(y), spbinBounds(y+1)]);
    end
    
    
    sptouchspks=zeros(length(spsorted),79);
    for j=1:length(spsorted)
        sptouchspks(j,2:77)=mean(spsorted{j},1);
        sptouchspks(j,1)=mean(spsortedBy{j});
        sptouchspks(j,78)=size(spsorted{j},1);
    end
    sptouchspks(:,1)=sprange;
    sptouchspks = sptouchspks(all(~isnan(sptouchspks),2),:); %remove NaN rows
    
    for k=1:size(sptouchspks,1)
        sptouchspks(k,2:77)=smooth(sptouchspks(k,2:77));
        if sptouchspks(k,78)>normbinmin
            sptouchspks(k,79)=max(sptouchspks(k,26:26+normbinwin));
        end
    end
    subplot(2,3,3);
    imagesc(sptouchspks(:,2:77))
    hold on
    caxis([0 max(thetatouchspks(:,79))]) %use this to set max spk rate scale of all plots
    plot([25 25],[25 0],'w:')
    colormap(gca,parula)
    ytickssp = cell(1,length(sptouchspks(:,1)));
    for i = 1:length(ytickssp)
        ytickssp{i} = sprintf('%0.2g',sptouchspks(i,1));
    end
    
    set(gca,'Ydir','normal','xtick',(0:25:75),'xticklabel',(-25:25:50),'ytick',[1:size(sptouchspks,1)],'yticklabel',ytickssp);
    for k=1:size(sptouchspks,1)
        text(20,k+.1,num2str(sptouchspks(k,78)),'FontSize',8,'Color','white')
    end
    axis('square')
    title('Setpoint at Touch')
    %% phase
    [phasesorted phasesortedBy phasebinBounds]=binslin(varspikes(:,4),varspikes(:,7:82),'equalX',13);
    
    phasetouchspks=zeros(length(phasesorted),79);
    for j=1:length(phasesorted)
        phasetouchspks(j,2:77)=mean(phasesorted{j},1);
        phasetouchspks(j,1)=mean(phasesortedBy{j});
        phasetouchspks(j,78)=size(phasesorted{j},1);
    end
    
    phasetouchspks = phasetouchspks(all(~isnan(phasetouchspks),2),:); %remove NaN rows
    
    for k=1:size(phasetouchspks,1)
        phasetouchspks(k,2:77)=smooth(phasetouchspks(k,2:77));
        if phasetouchspks(k,78)>normbinmin %if bin has more than 10 trials
            phasetouchspks(k,79)=max(phasetouchspks(k,26:26+normbinwin)); %find max spiking in each bin
        end
    end
    
    subplot(2,3,4);
    imagesc(phasetouchspks(:,2:77))
    hold on
    plot([25 25],[25 0],'w:')
    set(gca,'xtick',[0:5:75]);
    set(gca,'xticklabel',[-25:5:50]);
    title('Phase at Touch')
    %hCBar=colorbar;
    colormap(gca,parula)
    caxis([0 max(thetatouchspks(:,79))])
    set(gca,'Ydir','normal','ytick',[1:2.5:12],'yticklabel',{'-pi','-pi/2',0,'pi/2','pi'},'xtick',(0:25:75),'xticklabel',(-25:25:50))
    for k=1:size(phasetouchspks,1)
        text(20,k+.1,num2str(phasetouchspks(k,78)),'FontSize',8,'Color','white')
    end
    axis('square')
    %print(gcf,'-dpng',['Z:\Users\Jon\Projects\Characterization\' layer '\Figures\' cellcode '_' 'Variables_at_touch_HM'])
    %% max curvature
    [kappasorted kappasortedBy kappabinBounds]=binslin(varspikes(:,5),varspikes(:,7:82),'equalX',13);
    
    kappatouchspks=zeros(length(kappasorted),79);
    for j=1:length(kappasorted)
        kappatouchspks(j,2:77)=mean(kappasorted{j},1);
        kappatouchspks(j,1)=mean(kappasortedBy{j});
        kappatouchspks(j,78)=size(kappasorted{j},1);
    end
    
    kappatouchspks = kappatouchspks(all(~isnan(kappatouchspks),2),:); %remove NaN rows
    
    for k=1:size(kappatouchspks,1)
        kappatouchspks(k,2:77)=smooth(kappatouchspks(k,2:77));
        if kappatouchspks(k,78)>normbinmin
            kappatouchspks(k,79)=max(kappatouchspks(k,26:26+normbinwin));
        end
    end
    
    subplot(2,3,5);
    imagesc(kappatouchspks(:,2:77))
    hold on
    plot([25 25],[25 0],'w:')
    set(gca,'xtick',[0:5:75]);
    set(gca,'xticklabel',[-25:5:50]);
    title('Max Kappa at Touch')
    %hCBar=colorbar;
    colormap(gca,parula)
    caxis([0 max(thetatouchspks(:,79))])
    axis('square')
    
    ytickskap = cell(1,length(kappatouchspks(:,1)));
    for i = 1:length(ytickskap);
        ytickskap{i} = sprintf('%0.2g',kappatouchspks(i,1));
    end
    
    set(gca,'Ydir','normal','xtick',(0:25:75),'xticklabel',(-25:25:50),'ytick',[1:length(kappasortedBy)],...
        'yticklabel',ytickskap)
    for k=1:size(kappatouchspks,1)
        text(20,k+.1,num2str(kappatouchspks(k,78)),'FontSize',8,'Color','white')
    end
    
    %% pretouch velocity (average 3ms before touch)
    [pvelsorted pvelsortedBy pvelbinBounds]=binslin(varspikes(:,6),varspikes(:,7:82),'equalE',41,-10000,10000);
    pvelrange=[-9750:500:9750];
    
    pveltouchspks=zeros(length(pvelsorted),79);
    for j=1:length(pvelsorted)
        pveltouchspks(j,2:77)=mean(pvelsorted{j},1);
        pveltouchspks(j,1)=mean(pvelsortedBy{j});
        pveltouchspks(j,78)=size(pvelsorted{j},1);
    end
    
    pveltouchspks(:,1)=pvelrange;
    pveltouchspks = pveltouchspks(all(~isnan(pveltouchspks),2),:); %remove NaN rows
    
    
    for k=1:size(pveltouchspks,1)
        pveltouchspks(k,2:77)=smooth(pveltouchspks(k,2:77));
        if pveltouchspks(k,78)>normbinmin
            pveltouchspks(k,79)=max(pveltouchspks(k,26:26+normbinwin));
        end
    end
    
    subplot(2,3,6);
    imagesc(pveltouchspks(:,2:77))
    hold on
    plot([25 25],[25 0],'w:')
    set(gca,'xtick',[0:5:75]);
    set(gca,'xticklabel',[-25:5:50]);
    title('Vel Pre-Touch')
    axis('square')
    hCBar=colorbar;
    colormap(gca,parula)
    caxis([0 max(thetatouchspks(:,79))])
    
    
    yticksvel = cell(1,length(pveltouchspks(:,1)));
    for i = 1:length(yticksvel)
        yticksvel{i} = sprintf('%0.2g',pveltouchspks(i,1)/1000);
    end
    
    set(gca,'Ydir','normal','xtick',(0:25:75),'xticklabel',(-25:25:50),'ytick',[1:size(pveltouchspks,1)],'yticklabel',yticksvel);
    for k=1:size(pveltouchspks,1)
        text(20,k+.1,num2str(pveltouchspks(k,78)),'FontSize',8,'Color','white')
    end
    
    thetaall{p}=thetatouchspks;
    ampall{p}=amptouchspks;
    spall{p}=sptouchspks;
    phaseall{p}=phasetouchspks;
    
    %%
    set(gcf, 'Units', 'pixels', 'Position', [10, 100, 1000, 400]);
    %set(h1,'Units','pixels','Position',[
    print(gcf,'-dpng',['Z:\Users\Jon\Projects\Characterization\' layer '\Figures\' cellcode '_' 'Variables_at_touch_NORMALIZED'])
    
    t.varNames = {'touch theta','touch amplitude', 'touch setpoint', ...
        'touch phase', 'pre touch velocity', 'max kappa during touch'}
    
    %% FWHM and MODIDX
    thetasum=(mean(thetatouchspks(:,32:46),2));
    binthresh=sum(thetatouchspks(:,78))/numel(thetatouchspks(:,78))*.25;%min number of touches required in each bin
    thetasum=thetasum(thetatouchspks(:,78)>binthresh); %filter out bins with touches less than binthresh
    thetamod=(max(thetasum)-min(thetasum))/(max(thetasum)+min(thetasum));
    
    
end
