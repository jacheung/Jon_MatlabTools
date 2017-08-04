function [mask] = assist_touchmasks(array)    
% Function used to pull out masks based on touches and pole availability

% OUTPUT: 
% touchEx_mask = mask out all periods of touch 
% firsttouchEx_mask = mask out first touches wthin trials
% availtotrialend_mask = mask out periods from when pole rises to end of
% trial
% avail_mask = mask out periods from when pole is available 

%% this is used to mask out touches and allow analysis for phase, setpoint, and amplitude
    tmp=find(array.S_ctk(9,:,:)==1)/4000;
    d = tmp-floor(tmp); %keep only decimals
    avail=round(min(d)*4000); %pole available time in ms
    offset=round(array.meta.poleOffset*1000);
    offmax=find(offset>4000);
    offset(offmax)=4000;
    spikes = squeeze(array.R_ntk);
    
    firsttouchIdx = [find(array.S_ctk(9,:,:)==1)];
    firsttouchOffIdx = [find(array.S_ctk(10,:,:)==1)];
    touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
    touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
    if touchOffIdx<(numel(spikes)-70) %ensures that offset is within index boundaries
        touchOffIdx = touchOffIdx+70;
    end
    
    firsttouchEx_mask = ones(size(squeeze(array.S_ctk(1,:,:))));
    for i = 1:length(firsttouchIdx)
        firsttouchEx_mask(firsttouchIdx(i):firsttouchOffIdx(i))=NaN;
    end
    
    touchEx_mask = ones(size(squeeze(array.S_ctk(1,:,:))));
    for i = 1:length(touchOnIdx)
        touchEx_mask(touchOnIdx(i):touchOffIdx(i)) = NaN;
    end
    touchEx_mask(1:100,:) = 1; %since bleedover from end of trials before, tihs ensure we keep end
    
    availtotrialend_mask=ones(size(squeeze(array.S_ctk(1,:,:))));
    availtotrialend_mask(avail:size(availtotrialend_mask,1),:) = NaN; %block out all periods after pole availability
    
    avail_mask=ones(size(squeeze(array.S_ctk(1,:,:))));
    for j=1:array.k
        avail_mask(avail:offset(j),j)=NaN; %block out all periods between pole availability
    end
    
    sampling_mask=NaN(size(squeeze(array.S_ctk(1,:,:))));
    for f = 1:array.k
        sampling_mask(round(array.meta.poleOnset(f)*1000):round((array.meta.poleOnset(f)+.75)*1000),f)=1;
    end
   

    mask.touch = touchEx_mask;
    mask.first = firsttouchEx_mask;
    mask.availend = availtotrialend_mask;
    mask.avail = avail_mask;
    mask.samplingp = sampling_mask;
    mask.samp_notouch = (touchEx_mask.*sampling_mask);