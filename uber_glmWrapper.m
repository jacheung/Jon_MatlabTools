clear U
layer = 'L5b';
cellNum = [1:9 11:13 15 17:20 23:24 26 28:38 41:45];
trialCutoffs = repmat([1 200],numel(cellNum),1);
trialCutoffs(30,:)=[1 95];

%%
clear U
layer = 'L3';
cellNum = [3:8 11:14 16:22];
trialCutoffs = [1 200;1 200; 1 200; 1 200; 1 200; 1 200; 44 200; 1 200; 1 200; 1 200; 1 200; 1 200; 1 200; 1 200; 1 200; 1 200; 1 200];
%%
clear U
layer = 'L4';
cellNum = [1:6];
trialCutoffs = [1 200;1 200; 1 200; 1 200; 1 200; 1 200];
%%
clear U
layer = 'L3Out';
cellNum = [3:8];
trialCutoffs = [1 200;1 200; 1 200; 1 200; 1 200; 1 200];
%%
clear U
layer = 'L5bOut';
cellNum = [1:5 7];
trialCutoffs = repmat([1 200],numel(cellNum),1);
%%
clear U
layer = 'L5bInt';
cellNum = [1:7];
trialCutoffs = repmat([1 200],numel(cellNum),1);

%%
clear U
layer = 'BV';
cellNum = [1:8];
trialCutoffs = [148 348;190 390;214 414;89 289; 238 438 ; 233 433 ; 62 262 ; 262 462] ;

%%
clear U
layer = 'SM';
% cellNum = [2 11 16 23 27 34:38]; %used 13 insteda of 11 b/c no touches
% (same mouse though) 
cellNum = [2 13 16 23 27 34 36:38]; %removed 35 b/c no CR touches 
trialCutoffs = repmat([1 200],numel(cellNum),1);
%%
clear U
layer = 'N';
cellNum = [1:6];
trialCutoffs = repmat([1 800],numel(cellNum),1);
%% 

for cellStep = 1:length(cellNum)
    
    %loadSUDataFinalizedJC(cellNum(cellStep), SU)
     cd(['Z:\Users\Jon\Projects\Characterization\' layer '\TArrays'])
     load(['trial_array_' num2str(cellNum(cellStep)) '.mat'])
     cd(['Z:\Users\Jon\Projects\Characterization\' layer '\Contacts'])
     load(['ConTA_' num2str(cellNum(cellStep)) '.mat'])
        
    d.varNames = {'thetaAtBase', 'velocity', 'amplitude', 'setpoint', 'phase', ...
        'deltaKappa','M0Adj','FaxialAdj', 'firstTouchOnset', 'firstTouchOffset', ...
        'firstTouchAll', 'lateTouchOnset','lateTouchOffset','lateTouchAll','PoleAvailable','beamBreakTimes'}
    d.cellNum = cellNum(cellStep)
    
    d.t = max(cellfun(@(x)round(x.whiskerTrial.time{1}(end)*1000)+1,T.trials(T.whiskerTrialInds))); 
    
%     d.t = max([cellfun(@(x)numel(x.whiskerTrial.theta{1}),T.trials(T.whiskerTrialInds)) ...
%        cellfun(@(x)numel(x.whiskerTrial.time{1}),T.trials(T.whiskerTrialInds))] );
%    [~,useTrials] = intersect(T.trialNums,UgoodtrialnumsCrop{cellStep});
   [~,useTrials] = intersect(T.trialNums,T.whiskerTrialNums);
    useTrials = useTrials(useTrials >= trialCutoffs(cellStep,1) & useTrials <= trialCutoffs(cellStep,2));

    d.k = length(useTrials);
    d.u = 1;
    d.c = 16;
    d.S_ctk = nan(d.c, d.t, d.k);
    d.R_ntk = zeros(1, d.t, d.k);

    
useTrials = useTrials';
traj = 1;


    for i = 1:length(useTrials);
        display(i)
        timeIdx = round(T.trials{useTrials(i)}.whiskerTrial.time{traj}*1000)+1;
        
        
        theta = nan(1,d.t);
        theta(timeIdx) = T.trials{useTrials(i)}.whiskerTrial.thetaAtBase{traj};
        
        nanidx = find(isnan(theta));
        nanidx = nanidx(nanidx > 2 & nanidx < length(timeIdx)-2);
        for j = nanidx
            if j>=3999
                theta(j)=nanmean(theta(j+[-2:0]));
            else
                theta(j) = nanmean(theta(j+[-2:2]));
            end
        end
        
        firstTouchOn    = [];
        firstTouchOff   = [];
        firstTouchAll   = [];
        
        lateTouchOn     = [];
        lateTouchOff    = [];
        lateTouchAll    = [];
        
        if isfield(contacts{useTrials(i)},'segmentInds') 
            if ~isempty(contacts{useTrials(i)}.segmentInds{traj})
            firstTouchOn  = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(1,1)));
            firstTouchOff = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(1,2)));
            firstTouchAll = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(1,1):contacts{useTrials(i)}.segmentInds{traj}(1,2)));
            
            if size(contacts{useTrials(i)}.segmentInds{traj},1)>1
                lateTouchOn  = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(2:end,1)));
                lateTouchOff = round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(2:end,2)));
                
                for j = 2:size(contacts{useTrials(i)}.segmentInds{traj},1)
                lateTouchAll = cat(2,lateTouchAll,round(1000*T.trials{useTrials(i)}.whiskerTrial.time{traj}(contacts{useTrials(i)}.segmentInds{traj}(j,1):contacts{useTrials(i)}.segmentInds{traj}(j,2))));
                end
                
                end
            end 
        end
        
        [hh amplitude  filteredSignal setpoint amplitudeS setpointS phase phaseS] =  SAHWhiskerDecomposition(theta);
        
        travelIn = 0;
        travelOut = 250;
        
        vel = diff([0 theta])/.001;
        
        pinIn = round(1000*T.trials{useTrials(i)}.pinDescentOnsetTime + travelIn);
        pinOut = min([d.t round(1000*T.trials{useTrials(i)}.pinAscentOnsetTime + travelOut)]);
        
        d.S_ctk(1,:,i) = theta;
        d.S_ctk(2,:,i) = vel;
        d.S_ctk(3,:,i) = amplitude;
        d.S_ctk(4,:,i) = setpoint;
        d.S_ctk(5,:,i) = phase;
        d.S_ctk(6,timeIdx,i) = T.trials{useTrials(i)}.whiskerTrial.deltaKappa{traj};
        d.S_ctk(7,timeIdx,i) = contacts{useTrials(i)}.M0comboAdj{traj};
        d.S_ctk(8,timeIdx,i) = contacts{useTrials(i)}.FaxialAdj{traj};
        d.S_ctk(9,firstTouchOn,i) = 1;
        d.S_ctk(10,firstTouchOff,i) = 1;
        d.S_ctk(11,firstTouchAll,i) = 1;
             d.S_ctk(12,lateTouchOn,i) = 1;
        d.S_ctk(13,lateTouchOff,i) = 1;
        d.S_ctk(14,lateTouchAll,i) = 1;  
        d.S_ctk(15,:,i) = 0;
        d.S_ctk(15,pinIn:pinOut,i) = 1;
        d.S_ctk(16, ceil(1000*T.trials{useTrials(i)}.beamBreakTimes(T.trials{useTrials(i)}.beamBreakTimes > 0 & T.trials{useTrials(i)}.beamBreakTimes < d.t/1000)),i) = 1;

        d.S_ctk(6,setdiff(1:d.t,contacts{useTrials(i)}.contactInds{traj}),i) = 0;
        d.S_ctk(7,setdiff(1:d.t,contacts{useTrials(i)}.contactInds{traj}),i) = 0;
        d.S_ctk(8,setdiff(1:d.t,contacts{useTrials(i)}.contactInds{traj}),i) = 0;
        
        
        
        
        
    end
    
%end
    

    d.S_ctk(15,1:min(find(nansum(squeeze(d.S_ctk(9,:,:))')))-1,:)= 0;  % define pole availablitity onset
    for i =1:length(useTrials);
        
        spiketimes = round(T.trials{useTrials(i)}.spikesTrial.spikeTimes/10-T.whiskerTrialTimeOffset*1000);
        spiketimes = spiketimes(spiketimes>0 & spiketimes <=d.t);
        d.R_ntk(1,spiketimes,i) = 1;
    end
    
    U{cellStep} = d;
%     
%     
%     U{cellStep}.meta.SUnum             = cellNum(cellStep)
%     U{cellStep}.meta.cellName          = SU.cellName{cellNum(cellStep)};
%    % U{cellStep}.meta.recordingLocation = SU.recordingLocation{cellNum(cellStep)};
%     U{cellStep}.meta.trialArrayName    = SU.trialArrayName{cellNum(cellStep)}
%   %  U{cellStep}.meta.behaviorArrayName = SU.behaviorArrayName{cellNum(cellStep)};
%     U{cellStep}.meta.contactsArrayName = SU.contactsArrayName{cellNum(cellStep)};
    U{cellStep}.meta.motorPosition     = cellfun(@(x)x.behavTrial.motorPosition,T.trials(useTrials));
    U{cellStep}.meta.goPosition        = cellfun(@(x)x.behavTrial.goPosition,T.trials(useTrials));
    U{cellStep}.meta.nogoPosition      = cellfun(@(x)x.behavTrial.nogoPosition,T.trials(useTrials));
    U{cellStep}.meta.trialType         = cellfun(@(x)x.behavTrial.trialType,T.trials(useTrials));
    U{cellStep}.meta.trialCorrect      = cellfun(@(x)x.behavTrial.trialCorrect,T.trials(useTrials));
    U{cellStep}.meta.poleOnset         = cellfun(@(x)x.behavTrial.pinDescentOnsetTime,T.trials(useTrials));
    U{cellStep}.meta.poleOffset        = cellfun(@(x)x.behavTrial.pinAscentOnsetTime,T.trials(useTrials));
    U{cellStep}.meta.layer             = layer;
    U{cellStep}.meta.ranges            = unique([cellfun(@(x)x.behavTrial.nogoPosition,T.trials(useTrials)),cellfun(@(x)x.behavTrial.goPosition,T.trials(useTrials))]);
    U{cellStep}.meta.sweepArrayName    = {};%SU.sweepArrayName{cellNum(cellStep)};
    U{cellStep}.meta.harddrive         = {};
   % U{cellStep}.meta.C2distance        = sqrt(SU.recordingLocation{cellNum(cellStep)}(1)^2+SU.recordingLocation{cellNum(cellStep)}(2)^2)
  %  U{cellStep}.meta.isC2              = SU.distance{cellNum(cellStep)} < .16;
   % U{cellStep}.meta.isPV              = SU.isU{cellNum(cellStep)};
    U{cellStep}.meta.manipulation      = {};    
    U{cellStep}.meta.performingTrials  = [];
    U{cellStep}.meta.performance       = [];
    U{cellStep}.meta.stimTrials        = {};
    
    U{cellStep}.whisker.follicleX      = cellfun(@(x)x.whiskerTrial.follicleCoordsX,T.trials(useTrials));
    U{cellStep}.whisker.follicleY      = cellfun(@(x)x.whiskerTrial.follicleCoordsY,T.trials(useTrials));
    U{cellStep}.whisker.barPos         = cellfun(@(x)x.whiskerTrial.barPos,T.trials(useTrials),'UniformOutput',false);
end    

% end
% 
% 
