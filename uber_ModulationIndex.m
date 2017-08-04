
%Touch On Mod Idx
%% Touch v no touch FR
L3TotalPre=zeros(1,length(L3));
L3TotalPost=zeros(1,length(L3));
L3TotalPreOff=zeros(1,length(L3));
L3TotalPostOff=zeros(1,length(L3));
L3ONmodidx=zeros(1,length(L3));
L3OFFmodidx=zeros(1,length(L3));
for j = 1:length(L3)
    %ID Touch Onset FR
    touchIdx = [find(L3{j}.S_ctk(9,:,:)==1);find(L3{j}.S_ctk(12,:,:)==1)]; %all periods of touch
    touchIdx = touchIdx;
    spikesAligned = zeros(numel(touchIdx),201); %empty matrix for size
    spikes = squeeze(L3{j}.R_ntk);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-100:100]);
    end
    
    PreTouchSpikes = spikesAligned(:,1:100); %first 100ms BEFORE touch
    AvgPreTouchSpikes = mean(mean(PreTouchSpikes)); %population average pretouch
    L3TotalPre(j)=AvgPreTouchSpikes;    
    
    PostTouchSpikes=spikesAligned(:,101:150); %first 50ms AFTER touch
    PostTouchSpikesMean=mean(PostTouchSpikes);
    touchBin=zeros(1,size(PostTouchSpikes,2)/2);
    touchBin=(PostTouchSpikesMean(:,1:2:end-1)+PostTouchSpikesMean(:,2:2:end))/2;
    MaxPostTouchSpikes = (max(touchBin));
    L3TotalPost(j)=MaxPostTouchSpikes;
    
    %figure;
    %bar(-100:100,sum(spikesAligned)/numel(touchIdx),'k')    
    
    %ID Touch Offset FR  
    offtouchIdx = [find(L3{j}.S_ctk(10,:,:)==1);find(L3{j}.S_ctk(13,:,:)==1)]; %all periods of touch
    offtouchIdx = offtouchIdx;
    spikesAlignedoff = zeros(numel(offtouchIdx),201); %empty matrix for size
    spikes = squeeze(L3{j}.R_ntk);
    
    for i = 1:size(spikesAlignedoff,1)
        spikesAlignedoff(i,:) = spikes(offtouchIdx(i)+[-100:100]);
    end
       
    PreTouchOffSpikes = spikesAlignedoff(:,1:100); %first 100ms BEFORE touch offset
    AvgPreTouchOffSpikes = mean(mean(PreTouchOffSpikes)); %population average pretouch
    L3TotalPreOff(j)=AvgPreTouchOffSpikes;    
    
    PostTouchOffSpikes=spikesAlignedoff(:,101:150); %first 50ms AFTER touch offset
    PostTouchOffSpikesMean=mean(PostTouchOffSpikes);
    touchBin=zeros(1,size(PostTouchOffSpikes,2)/2);
    touchBin=(PostTouchOffSpikesMean(:,1:2:end-1)+PostTouchOffSpikesMean(:,2:2:end))/2;
    MaxPostTouchOffSpikes = (max(touchBin));
    L3TotalPostOff(j)=MaxPostTouchOffSpikes;
end

    L3TotalPost=L3TotalPost*1000; %convert all to spks/s
    L3TotalPre=L3TotalPre*1000;
    L3TotalPreOff=L3TotalPreOff*1000;
    L3TotalPostOff=L3TotalPostOff*1000;
%Modulation Indx    
for k=1:length(L3)
    L3ONmodidx(k) = (L3TotalPost(k)-L3TotalPre(k))/(L3TotalPost(k)+L3TotalPre(k));
    L3OFFmodidx(k) = (L3TotalPostOff(k)-L3TotalPreOff(k))/(L3TotalPostOff(k)+L3TotalPreOff(k));
end

%% Whisk v no Whisk FR (excluding periods of touch) 
L3TotalHighWhisk=zeros(1,length(L3));
L3TotalLowWhisk=zeros(1,length(L3));
L3WhiskModIdx=zeros(1,length(L3));

%Theta
ThetaHigh=zeros(1,length(L3));
ThetaLow=zeros(1,length(L3));
L3Thetamodix=zeros(1,length(L3));

%SetPoint
SPHigh=zeros(1,length(L3));
SPLow=zeros(1,length(L3));
L3SPmodix=zeros(1,length(L3));
for j=1:length(L3)
        
        touchOnIdx = [find(L3{j}.S_ctk(9,:,:)==1); find(L3{j}.S_ctk(12,:,:)==1)];
        %touchOnIdx = touchOnIdx+70 %do I need to add this to mask out spikes for
        %periods post touch?
        touchOffIdx = [find(L3{j}.S_ctk(10,:,:)==1); find(L3{j}.S_ctk(13,:,:)==1)];
        touchOffIdx = touchOffIdx+70;

        touchEx_mask = ones(size(squeeze(L3{j}.S_ctk(1,:,:))));

        for i = 1:length(touchOnIdx)
            touchEx_mask(touchOnIdx(i):touchOffIdx(i)) = NaN;
        end
        touchEx_mask(1:100,:) = 1;
    %Whisk    
    highamp_mask=squeeze(L3{j}.S_ctk(3,:,:)>2.5); %ID amplitudes of high whisking 
    selectedHighamp = touchEx_mask.*highamp_mask; %mask out touches for periods of contact during high whisk
    highamp_idx=find(~isnan(selectedHighamp));
    highamp_idxone = find(selectedHighamp(highamp_idx) ~= 0);
    highampSpikes = selectedHighamp(highamp_idxone);
    highampSpikes(:,2) = L3{j}.R_ntk(highamp_idxone);

    lowamp_mask=squeeze(L3{j}.S_ctk(3,:,:)<1.5); %ID amplitudes of no whisking 
    selectedLowamp = touchEx_mask.*lowamp_mask; %mask out touches
    lowamp_idx=find(~isnan(selectedLowamp));
    lowamp_idxone = find(selectedLowamp(lowamp_idx) ~= 0);
    lowampSpikes = selectedLowamp(lowamp_idxone);
    lowampSpikes(:,2) = L3{j}.R_ntk(lowamp_idxone);

    L3TotalHighWhisk(j)=mean(highampSpikes(:,2));
    L3TotalLowWhisk(j)=mean(lowampSpikes(:,2));
    
    %Theta
    theta = squeeze(L3{j}.S_ctk(1,:,:));
    selectedtheta = touchEx_mask.*theta;
    theta_Idx = find(~isnan(selectedtheta));
    thetaSpikes = selectedtheta(theta_Idx);
    thetaSpikes(:,2) = L3{j}.R_ntk(theta_Idx);
        
    [sorted sortedBy binBounds]=binslin(thetaSpikes(:,1),thetaSpikes(:,2),'equalX',13);

    theta_sp = zeros(12,1);
    for i = 1:12;
        theta_sp(i) = sum(sorted{i})/numel(sorted{i})*1000;% 160125 changed back to 1000 from 100(which was stupid to even change to...)since we want it to be spks/s
    end
    [~,ic]=min(abs(theta_sp(theta_sp>0))); 
    ThetaLow(j)=theta_sp(ic);
    [~,ib]=max(abs(theta_sp)); 
    ThetaHigh(j)=theta_sp(ib);
    
    %Setpoint
    setpoint = squeeze (L3{j}.S_ctk(4,:,:));
    selectedSP = touchEx_mask.*setpoint; 
    SP_Idx = find(~isnan(selectedSP));
    SPSpikes = selectedSP(SP_Idx);
    SPSpikes(:,2) = L3{j}.R_ntk(SP_Idx);

    [sorted sortedBy binBounds]=binslin(SPSpikes(:,1),SPSpikes(:,2),'equalX',13);
    FR_sp = zeros(12,1);

    for i = 1:12;
        FR_sp(i) = sum(sorted{i})/numel(sorted{i})*1000; %1000 for ocnversoin to spks/s
    end

    [~,ic]=min(abs(FR_sp(theta_sp>0))); %finds index for where max abs value
    SPLow(j)=FR_sp(ic);%sets ThetaLow to that value based on index (this ensures you get that value +/-)
    [~,ib]=max(abs(FR_sp)); 
    SPHigh(j)=FR_sp(ib); 
end
   L3TotalHighWhisk=L3TotalHighWhisk*1000;%conversion to spks/s
   L3TotalLowWhisk=L3TotalLowWhisk*1000;
   
%Modulation Index
for k=1:length(L3)
    L3WhiskModIdx(k)=(L3TotalHighWhisk(k)-L3TotalLowWhisk(k))/(L3TotalHighWhisk(k)+L3TotalLowWhisk(k));
    L3Thetamodidx(k)=(ThetaHigh(k)-ThetaLow(k))/(ThetaHigh(k)+ThetaLow(k));
    L3SPmodidx(k)=(SPHigh(k)-SPLow(k))/(SPHigh(k)+SPLow(k));
end

%% Concatenate all
m= 255;
primemap = zeros(m , 3);
t= rgb('black'); %default midline color

color1=rgb('red');
color2=rgb('blue');

%%%%
extremes = linspace(0,1,127)';
midpoint=1;


%%%color 1%%%%
primemap(1:127, 1) = extremes;

R = linspace(color2(1),t(1),127);
G = linspace(color2(2),t(2),127);
B = linspace(color2(3),t(3),127);

T = [R', G', B'];
primemap(1:127,:) = T;

%%%%color 2%%%%

%%%%
primemap(129:255,1)=extremes;

R = linspace(t(1),color1(1),127);
G = linspace(t(2),color1(2),127);
B = linspace(t(3),color1(3),127);

T = [R', G', B'];

primemap(129:255,:) = T;
%%%%

T = [0,   0,   1;   %// color1
     0,   0,   0;   %// midline color
     1,   0,   0]; %// color2  
 

 x = [0
     128
     255];
 
 
primemap = interp1(x/255,T,linspace(0,1,255));

%L3donuts=cat(1,L3ONmodidx,L3OFFmodidx,L3WhiskModIdx,L3Thetamodidx,L3SPmodidx);
L3donuts=cat(1,L3ONmodidx,L3OFFmodidx,L3WhiskModIdx);
%colormap('parula')
colormap(primemap);
imagesc(L3donuts);
hCBar=colorbar;
set(hCBar,'yaxislocation','right');
set(hCBar,'YTick',[1:1:10]);
caxis([-1 1]);
title('L3 Modulation Index');
set(gca,'YTick',(1:1:5));
set(gca,'YTickLabel',{'Touch Onset', 'Touch Offset','Whisking Amplitude','Whisking Position','Whisking Setpoint'});
set(gca,'XTickLabel',{});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Repeat Above but for L5b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Touch v no touch FR
L5bTotalPre=zeros(1,length(L5b));
L5bTotalPost=zeros(1,length(L5b));
L5bTotalPreOff=zeros(1,length(L5b));
L5bTotalPostOff=zeros(1,length(L5b));
L5bONmodidx=zeros(1,length(L5b));
L5bOFFmodidx=zeros(1,length(L5b));
for j = 1:length(L5b)
    %ID Touch Onset FR
    touchIdx = [find(L5b{j}.S_ctk(9,:,:)==1);find(L5b{j}.S_ctk(12,:,:)==1)]; %all periods of touch
    touchIdx = touchIdx;
    spikesAligned = zeros(numel(touchIdx),201); %empty matrix for size
    spikes = squeeze(L5b{j}.R_ntk);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-100:100]);
    end
    
    PreTouchSpikes = spikesAligned(:,1:100); %first 100ms BEFORE touch
    AvgPreTouchSpikes = mean(mean(PreTouchSpikes)); %population average pretouch
    L5bTotalPre(j)=AvgPreTouchSpikes;    
    
    PostTouchSpikes=spikesAligned(:,101:150); %first 50ms AFTER touch
    PostTouchSpikesMean=mean(PostTouchSpikes);
    touchBin=zeros(1,size(PostTouchSpikes,2)/2);
    touchBin=(PostTouchSpikesMean(:,1:2:end-1)+PostTouchSpikesMean(:,2:2:end))/2;
    MaxPostTouchSpikes = (max(touchBin));
    L5bTotalPost(j)=MaxPostTouchSpikes;
    
    %figure;
    %bar(-100:100,sum(spikesAligned)/numel(touchIdx),'k')    
    
    %ID Touch Offset FR  
    offtouchIdx = [find(L5b{j}.S_ctk(10,:,:)==1);find(L5b{j}.S_ctk(13,:,:)==1)]; %all periods of touch
    offtouchIdx = offtouchIdx;
    spikesAlignedoff = zeros(numel(offtouchIdx),201); %empty matrix for size
    spikes = squeeze(L5b{j}.R_ntk);
    
    for i = 1:size(spikesAlignedoff,1)
        spikesAlignedoff(i,:) = spikes(offtouchIdx(i)+[-100:100]);
    end
       
    PreTouchOffSpikes = spikesAlignedoff(:,1:100); %first 100ms BEFORE touch offset
    AvgPreTouchOffSpikes = mean(mean(PreTouchOffSpikes)); %population average pretouch
    L5bTotalPreOff(j)=AvgPreTouchOffSpikes;    
    
    PostTouchOffSpikes=spikesAlignedoff(:,101:150); %first 50ms AFTER touch offset
    PostTouchOffSpikesMean=mean(PostTouchOffSpikes);
    touchBin=zeros(1,size(PostTouchOffSpikes,2)/2);
    touchBin=(PostTouchOffSpikesMean(:,1:2:end-1)+PostTouchOffSpikesMean(:,2:2:end))/2;
    MaxPostTouchOffSpikes = (max(touchBin));
    L5bTotalPostOff(j)=MaxPostTouchOffSpikes;
end

    L5bTotalPost=L5bTotalPost*1000; %convert all to spks/s
    L5bTotalPre=L5bTotalPre*1000;
    L5bTotalPreOff=L5bTotalPreOff*1000;
    L5bTotalPostOff=L5bTotalPostOff*1000;
%Modulation Indx    
for k=1:length(L5b)
    L5bONmodidx(k) = (L5bTotalPost(k)-L5bTotalPre(k))/(L5bTotalPost(k)+L5bTotalPre(k));
    L5bOFFmodidx(k) = (L5bTotalPostOff(k)-L5bTotalPreOff(k))/(L5bTotalPostOff(k)+L5bTotalPreOff(k));
end

%% Whisk v no Whisk FR (excluding periods of touch) 
L5bTotalHighWhisk=zeros(1,length(L5b));
L5bTotalLowWhisk=zeros(1,length(L5b));
L5bWhiskModIdx=zeros(1,length(L5b));

%Theta
ThetaHigh=zeros(1,length(L5b));
ThetaLow=zeros(1,length(L5b));
L5bThetamodix=zeros(1,length(L5b));

%SetPoint
SPHigh=zeros(1,length(L5b));
SPLow=zeros(1,length(L5b));
L5bSPmodix=zeros(1,length(L5b));
for j=1:length(L5b)
        
        touchOnIdx = [find(L5b{j}.S_ctk(9,:,:)==1); find(L5b{j}.S_ctk(12,:,:)==1)];
        %touchOnIdx = touchOnIdx+70 %do I need to add this to mask out spikes for
        %periods post touch?
        touchOffIdx = [find(L5b{j}.S_ctk(10,:,:)==1); find(L5b{j}.S_ctk(13,:,:)==1)];
        touchOffIdx = touchOffIdx+70;

        touchEx_mask = ones(size(squeeze(L5b{j}.S_ctk(1,:,:))));

        for i = 1:length(touchOnIdx)
            touchEx_mask(touchOnIdx(i):touchOffIdx(i)) = NaN;
        end
        touchEx_mask(1:100,:) = 1;
    %Whisk    
    highamp_mask=squeeze(L5b{j}.S_ctk(3,:,:)>2.5); %ID amplitudes of high whisking 
    selectedHighamp = touchEx_mask.*highamp_mask; %mask out touches for periods of contact during high whisk
    highamp_idx=find(~isnan(selectedHighamp));
    highamp_idxone = find(selectedHighamp(highamp_idx) ~= 0);
    highampSpikes = selectedHighamp(highamp_idxone);
    highampSpikes(:,2) = L5b{j}.R_ntk(highamp_idxone);

    lowamp_mask=squeeze(L5b{j}.S_ctk(3,:,:)<1.5); %ID amplitudes of no whisking 
    selectedLowamp = touchEx_mask.*lowamp_mask; %mask out touches
    lowamp_idx=find(~isnan(selectedLowamp));
    lowamp_idxone = find(selectedLowamp(lowamp_idx) ~= 0);
    lowampSpikes = selectedLowamp(lowamp_idxone);
    lowampSpikes(:,2) = L5b{j}.R_ntk(lowamp_idxone);

    L5bTotalHighWhisk(j)=mean(highampSpikes(:,2));
    L5bTotalLowWhisk(j)=mean(lowampSpikes(:,2));
    
    %Theta
    theta = squeeze(L5b{j}.S_ctk(1,:,:));
    selectedtheta = touchEx_mask.*theta;
    theta_Idx = find(~isnan(selectedtheta));
    thetaSpikes = selectedtheta(theta_Idx);
    thetaSpikes(:,2) = L5b{j}.R_ntk(theta_Idx);
        
    [sorted sortedBy binBounds]=binslin(thetaSpikes(:,1),thetaSpikes(:,2),'equalX',13);

    theta_sp = zeros(12,1);
    for i = 1:12;
        theta_sp(i) = sum(sorted{i})/numel(sorted{i})*1000;% 160125 changed back to 1000 from 100(which was stupid to even change to...)since we want it to be spks/s
    end
    [~,ic]=min(abs(theta_sp(theta_sp>0))); 
    ThetaLow(j)=theta_sp(ic);
    [~,ib]=max(abs(theta_sp)); 
    ThetaHigh(j)=theta_sp(ib);
    
    %Setpoint
    setpoint = squeeze (L5b{j}.S_ctk(4,:,:));
    selectedSP = touchEx_mask.*setpoint; 
    SP_Idx = find(~isnan(selectedSP));
    SPSpikes = selectedSP(SP_Idx);
    SPSpikes(:,2) = L5b{j}.R_ntk(SP_Idx);

    [sorted sortedBy binBounds]=binslin(SPSpikes(:,1),SPSpikes(:,2),'equalX',13);
    FR_sp = zeros(12,1);

    for i = 1:12;
        FR_sp(i) = sum(sorted{i})/numel(sorted{i})*1000; %1000 for ocnversoin to spks/s
    end

    [~,ic]=min(abs(FR_sp(theta_sp>0))); %finds index for where max abs value
    SPLow(j)=FR_sp(ic);%sets ThetaLow to that value based on index (this ensures you get that value +/-)
    [~,ib]=max(abs(FR_sp)); 
    SPHigh(j)=FR_sp(ib); 
end
   L5bTotalHighWhisk=L5bTotalHighWhisk*1000;%conversion to spks/s
   L5bTotalLowWhisk=L5bTotalLowWhisk*1000;
   
%Modulation Index
for k=1:length(L5b)
    L5bWhiskModIdx(k)=(L5bTotalHighWhisk(k)-L5bTotalLowWhisk(k))/(L5bTotalHighWhisk(k)+L5bTotalLowWhisk(k));
    L5bThetamodidx(k)=(ThetaHigh(k)-ThetaLow(k))/(ThetaHigh(k)+ThetaLow(k));
    L5bSPmodidx(k)=(SPHigh(k)-SPLow(k))/(SPHigh(k)+SPLow(k));
end

%% Concatenate all
figure;
%L5bdonuts=cat(1,L5bONmodidx,L5bOFFmodidx,L5bWhiskModIdx,L5bThetamodidx,L5bSPmodidx);
L5bdonuts=cat(1,L5bONmodidx,L5bOFFmodidx,L5bWhiskModIdx);
colormap(primemap)
imagesc(L5bdonuts)
hCBar=colorbar
get(hCBar)
set(hCBar,'yaxislocation','right')
set(hCBar,'YTick',[1:1:10])
caxis([-1 1])
title('L5b Modulation Index')
set(gca,'YTick',(1:1:5));
set(gca,'YTickLabel',{'Touch Onset', 'Touch Offset','Whisking Amplitude','Whisking Position','Whisking Setpoint'});
set(gca,'XTickLabel',{});


%Full concact
Alldonuts=cat(2,L3donuts,L5bdonuts);
figure;
colormap(primemap)
imagesc(Alldonuts)
hCBar=colorbar
get(hCBar)
set(hCBar,'yaxislocation','right')
set(hCBar,'YTick',[1:1:10])
caxis([-1 1])
title('Modulation Index')
set(gca,'YTick',(1:1:5));
set(gca,'YTickLabel',{'Touch Onset', 'Touch Offset','Whisking Amplitude','Whisking Position','Whisking Setpoint'});
set(gca,'XTickLabel',{});