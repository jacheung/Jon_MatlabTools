%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will find the peak protraction of each whisk cycle within
% specific timepoints of a trial. I.e. Default is pole avail to first lick.
% 
% maskString vals = 'avail2lick' OR 'sampling'
%
% Will need to edit if you want to swap windows of where you want to look
% for peak protractions. 
%
%
%

function [P] = findMaxProtraction(array,maskString)

%% Find max theta excursion on each whisk cycle
    
    [objmask]= assist_touchmasks(array);
    
    if strcmp(maskString,'avail2lick')
        mask = objmask.availtolick;
    elseif strcmp(maskString,'sampling')
        mask = objmask.samplingp;
    end
    
    
    amp_mask = ones(size(squeeze(array.S_ctk(3,:,:))));
    amp_mask(array.S_ctk(3,:,:)<2.5) = NaN; %amplitude mask used to filter out only periods of high whisking
    phase = squeeze(array.S_ctk(5,:,:));
    amp = squeeze(array.S_ctk(3,:,:));
    setpoint = squeeze(array.S_ctk(4,:,:));
    theta = squeeze(array.S_ctk(1,:,:));
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
    
    P.peakidx=unique(phases+maxidx'); %this is the idx for peak protraction 
    P.trialNums = floor(P.peakidx/array.t)+1;
    P.theta = [theta(P.peakidx)];
    P.phase = [phase(P.peakidx)];
    P.amp= [amp(P.peakidx)];
    P.setpoint = [setpoint(P.peakidx)];
