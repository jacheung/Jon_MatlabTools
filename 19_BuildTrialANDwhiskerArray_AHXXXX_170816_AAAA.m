%%
numOfCore = feature('numcores')
parpool(numOfCore);


%% DONT FORGET TO CLEAR clear all VARIBALES WHEN DOING MULTIPLE BUILDS DONT WANT RUN OVER FROM 
%PROJECT INFO 
cellNumberForProject = '19';
projectName = upper('S2_PHIL');
mouseName  =  'AHXXXX'; 
sessionName = '170816'; %VIDEO AND BEHAVIOR 
code        = 'AAAA'; %EPHUS
cellnum     = 'XXXXXX'; %EPHUS

%RECORDING LOCATION INFO 
depth = 847;  % in um from pia
recordingLocation = [0 0 0];  % % in mm


%GENERAL INFO 
category1 = lower('');
category2 = lower('');
category3 = lower('');
category4 = lower('');
category5 = lower('');
cellNotes = 'AAAA1 DOUBLE CELL!!!! 1/2';

%TRIAL INFO adjust if ephus files are bad 
sweepnums =  [2:685];
% sweepnums =  [174:202, 204:229,231:274, 276:285,287:288,290:380];

%PERFORMANCE INFO     
prePerformanceRegion = [];
goodPerformanceRegion = [];
postPerformanceRegion = [521:686]; 

% SPIKE INFO
spikesNormalRegion = [];
spikesFastBrokeInCellRegion = [];



%VIDEO INFO 
NASletter = 'Z' ; %%%letter for which nas to use either z or y capitalized (for video)
videoloc    = 'PHILLIP';%%
whiskerTrialTimeOffset = .15; % in seconds to account for frametrig lag


%% Build Behavior
fn =['Z:\Users\Phil\Data\Behavior\' mouseName filesep 'data_@pole_contdiscrim_obj_' mouseName '_' sessionName 'a.mat'];
b = Solo.BehavTrialArray(fn, sessionName,'loadingAutosaveFile'); 
%%%% NOTE YOU CAN ADD ANOTHER VARIABLE INPUT ('loadingAutosaveFile') AND THE PROGRAM 
%%%% WILL SKIP SKIP THE LAST TRIAL, I DONT KNOW HOW THIS WILL IMPACT LATER
%%%% BUT SHOULD BE FINE BECASUE WE USE 'INTERSECT' BELOW TO REMOVE NON
%%%% INTERSECTING TRIALS
%%%% ANALYSIS THOUGH !!!!!!!
b.trim = [1 0];  
b.performanceRegion = sweepnums; 
figure; 
b.plot_scored_trials
pc = performance(b)
save(['Z:\Users\Phil\Data\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'],'b')

%% Build Spike Array

cd(['Z:\Users\Phil\Data\EPHUS\' mouseName]); 
s = LCA.SweepArray(cellnum, code, sweepnums);
%% set threshold 

 s = s.set_primary_threshold_all(1); % in mV
% % % % % % % 
% % % % % % % threshVars =       [.2    .4      .8     1      .8      .5          .3         ];
% % % % % % % trialVars =     [2     48     100    300    500      622        642    686   ];%first number should be first trial and last should be last will be 1 larger than threshVars
% % % % % % % for k = 1:numel(threshVars)
% % % % % % % s = s.set_primary_threshold([trialVars(k):trialVars(k+1)],repmat(threshVars(k),size([trialVars(k):trialVars(k+1)]))); % trials and mV
% % % % % % % end




% % % % % s = s.set_primary_threshold(350:499,repmat(1.4,size(350:499)));
% % % % % s = s.set_primary_threshold(500:650,repmat(2,size(500:650)));
% % % % % s = s.set_primary_threshold(651:680,repmat(.8,size(651:680)));
% % % % % s = s.set_primary_threshold(681:713,repmat(.45,size(681:713)));

% s = s.set_artifact_threshold_all(-1);
%  s = s.set_artifact_threshold(78:150,repmat(-.8,size(78:150)));
%s.sweeps{15}.trialNum=231 %% used for when there is a glitch and one of the
%trial numbers are not assigned for the sweep
s.viewer


    
%% USE THIS TO SET UP GOOD SPIKES FOR BELOW, OLD CODE THAT USED TO REMOVED 
%%%% LICKS BUT NOW JUST PACKAGES SPIKES FOR LATER...
% by filtering
sweepTrialNums = cellfun(@(x)x.trialNum,s.sweeps);
sweepTrialNums(find(sweepTrialNums==0))= sweepTrialNums(find(sweepTrialNums==0)+1)-1;
[c, ia, ib] = intersect(b.trialNums,sweepTrialNums);
sweepBeamBreakTimes = cellfun(@(x)x.beamBreakTimes,b.trials(ia),'uniformoutput',0);
for i = 1:length(s.sweeps)
% % % % % % % % % % % % % % % % % % % % % %         [indx,indy] = find(abs(repmat(sweepBeamBreakTimes{i}*10000,1,length(s.sweeps{i}.spikeTimes))...
% % % % % % % % % % % % % % % % % % % % % %                               -repmat(s.sweeps{i}.spikeTimes',length(sweepBeamBreakTimes{i}),1))<100);
    badSpikes{i} = [];
    goodSpikes{i} = setdiff(find(s.sweeps{i}.spikeTimes),badSpikes{i});
    %   s.sweeps{i}.spikeTimes = s.sweeps{i}.spikeTimes(goodSpikes{i})
    %     s.sweeps{i}.spikeWaveforms{1} = s.sweeps{i}.spikeWaveforms{1}(goodSpikes{i});
    %     s.sweeps{i}.spikesWaveforms{2} = s.sweeps{i}.spikeWaveforms{2}(goodSpikes{i});
end


%% section for EVEN WORSE licking artifact where you will need to specify
% % % % % % % % % %% values for defining what licks are
% % % % % % % % % 
% % % % % % % % % sw = s.get_spike_waveforms_all;
% % % % % % % % % swall = cellfun(@(x)x.spikeWaveforms{2},sw,'uniformoutput',0);
% % % % % % % % % bad = cell(size(swall));
% % % % % % % % % good = cell(size(swall));
% % % % % % % % % rise=10; %AP rise time in ms
% % % % % % % % % 
% % % % % % % % % %use this to plot BEFORE
% % % % % % % % % swall = cellfun(@(x)x.spikeWaveforms{2},sw,'uniformoutput',0);
% % % % % % % % % swcat = [swall{:}];
% % % % % % % % % swscale = (swcat - repmat(mean(swcat(1:rise,:)),30,1))./repmat(max(swcat)-mean(swcat(1:rise,:)),30,1);
% % % % % % % % % %normalized to min of trial swscale = (swcat - repmat(min(swcat),30,1))./repmat( max(swcat)-min(swcat),30,1);
% % % % % % % % % figure;plot(swscale)
% % % % % % % % % title('unfiltered')
% % % % % % % % % 
% % % % % % % % % for i=1:length(swall) %built to filter out each trial individually so can cross validate w/ licks 
% % % % % % % % %     swcat = swall{i};
% % % % % % % % %     if ~isempty(swall{i});
% % % % % % % % %         swscale = (swcat - repmat(mean(swcat(1:rise,:)),30,1))./repmat(max(swcat)-mean(swcat(1:rise,:)),30,1);
% % % % % % % % %          %bad{i} = find(swscale(1,:)<-999999999);%no filter
% % % % % % % % %        bad{i} = find(swscale(1,:)<-.2|swscale(22,:)<-.8|swscale(19,:)>.25); %this is what you change to filter out bad spikes
% % % % % % % % %         good{i} = setdiff(1:size(swscale,2),bad{i});
% % % % % % % % %         newswall{i} = swcat(:,good{i}); %rewrites swall to keep only good spike columns
% % % % % % % % %     end
% % % % % % % % % end
% % % % % % % % % 
% % % % % % % % % % use this graph to compare to your swscale and see if these are the spike
% % % % % % % % % % you want. If not, change the filter settings through "bad" 
% % % % % % % % % swcatall = [newswall{:}];
% % % % % % % % % scatscale = (swcatall - repmat(mean(swcatall(1:rise,:)),30,1))./repmat(max(swcatall)-mean(swcatall(1:rise,:)),30,1);
% % % % % % % % % figure; plot(scatscale)
% % % % % % % % % title('Test Filter')

 
%%  use this to get the spike array setup baby
sw = s.get_spike_waveforms_all;
swall = cellfun(@(x)x.spikeWaveforms{2},sw,'uniformoutput',0);
bad = cell(size(swall));
good = cell(size(swall));
rise=10; % change to change where the spikes are aligned to (from 1 to rise mean...)
[trash, spikeCountPerTrial] = cellfun(@size, swall);
%% use this to make single array of spikes
swcat = [swall{:}];
swscale = (swcat - repmat(mean(swcat(1:rise,:)),30,1))./repmat(max(swcat)-mean(swcat(1:rise,:)),30,1);
% swaligned = swcatall - repmat(sum(swcatall(1:10,:))/10,30,1);
swaligned = swcat - repmat(sum(swcat(1:10,:))/10,30,1);
swaligned12 = swcat - repmat(swcat(12,:),30,1); %aligh to time = 12 ms 
swaligned13 = swcat - repmat(swcat(13,:),30,1); %aligh to time = 13 ms 
swaligned15 = swcat - repmat(swcat(15,:),30,1); %aligh to time = 15 ms 
%normalized to min of trial swscale = (swcat - repmat(min(swcat),30,1))./repmat(max(swcat)-min(swcat),30,1);
size(swscale)
%% view all spikes (dont really need will do below)
figure;
plot(swscale)
title('unfiltered scale (normalized)')
figure;
plot(swaligned)
title('unfiltered aligned (not normalized)')

%% use this to sort spikes 
% % BELOW IS JUST AN EXAMPLE OF THE PROGRAM
% % [aligned,spksToPlotGood,spksToPlotBad]= ...
% % spike_sorting_PSM(TrialIndCount,swall,scaleORaligned, preSpikeBaselineTime,varargin)

preSpikeBaselineTime = (10:12); % set to align the spikes from mean of ms X1 to ms X2

TrialIndCount = [(1:1000:2000)]

trialsPerPlot = 25; 
%%% DONT CHANGE THIS VVVVVV THE PROGRAM WILL AUTO DETECT END OF SPIKES
TrialIndCount =(1:trialsPerPlot:20000);%%program will autotrim this don't worry


 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%USE ON OF THE BELOW 2 FUNCTIONS%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% you can do this multiple times and cross reference the good
%%%%%%% vairable (saving it as a different variable like good2 good3 ect.
%%%%%%% so that you can take advantage of both scale and aligned setup just
%%%%%%% make sure to use set difference below 
%% THIS IS FOR NEW SETS 

setToScaleOrAligned = 'scale'; 
[indexingMat,spksToPlotGood,spksToPlotBad, good, bad]= ...
spike_sorting_PSM(TrialIndCount,swall,setToScaleOrAligned, preSpikeBaselineTime);

%% THIS IS FOR IF YOU HAVE AN 'indexingMat' ALSO LOADS TrialIndCount 
startOnPlotNum = 1;
setToScaleOrAligned = 'aligned'; 
[indexingMat,spksToPlotGood,spksToPlotBad, good, bad]= ...
spike_sorting_PSM(TrialIndCount,swall,setToScaleOrAligned, preSpikeBaselineTime,indexingMat,startOnPlotNum);
%%%%NOTE THIS WILL LOAD ALL THE POINTS FROM THE INDEXINGMAT SO IF FOR
%%%%EXAMPLE YOU HAVE 4 POINTS ON PLOT NUMBER 10, AND YOU WANT TO DO PLOT 10
%%%%OVER, YOU HAVE TO CLICK 4 TIMES TO OVERWRITE ALL THE PREVIOUS POINTS.
%%%%JUST USE 'DUMMY POINTS' IF NEEDED (POINTS OUTSIDE OF SPIKES THAT DONT
%%%%REMOVE ANY SPIKES) (ALSO NOTE ALL THE POINTS WILL BE VISABLE SO YOU
%%%%SHOULDNT ACCIDNETALLY MISS THIS IF YOU PAY ATTENTION TO THE PLOTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% for multiple good arrays 
%%%%use 1 for aligned and 2 for scale 
for k = 1
goodArray{k} = good;
badArray{k} = bad;
indexingMatArray{k} = indexingMat;
end
%% combine good and bad arrays 
%setdif(A, B) returns data in A that isn't in B
for k=1:length(s.sweeps)
        goodIntersect{k} = intersect(goodArray{1}{k}',goodArray{2}{k}');
        good1diff{k} = setdiff(goodArray{1}{k}', goodArray{2}{k}');
        good2diff{k} = setdiff(goodArray{2}{k}', goodArray{1}{k}');%good spikes that are in 2 but not 1
        bad1diff{k} = setdiff(badArray{1}{k}', badArray{2}{k}');
        bad2diff{k} = setdiff(badArray{2}{k}', badArray{1}{k}');%bad spikes in 2 that are not in 1
end


%% plot the above spikes to figure out which is better

selectedArray = goodSpikes;%set selcted array to look at  
clear newswall
for k = 1:numel(swall)
newswall{k} = swcat(:,selectedArray{k}); 
end
selectedArraySize = sum(cell2mat(cellfun(@numel,selectedArray, 'UniformOutput',false)))
allselectedArray = [newswall{:}];
allselectedArrayAligned =  allselectedArray - repmat(sum(allselectedArray(preSpikeBaselineTime,:))/numel(preSpikeBaselineTime),30,1);
allselectedArrayScale = (allselectedArray - repmat(mean(allselectedArray(1:rise,:)),30,1))./repmat(max(allselectedArray)-mean(allselectedArray(preSpikeBaselineTime,:)),30,1);

%% look at selected array spikes 
figure
selectedArraySize
numSpikesPerPlot = 1000;
iters = 1:numSpikesPerPlot:selectedArraySize;
if iters(end)~=selectedArraySize
    iters(end+1)=selectedArraySize;
end
numel(iters)-1
for i = 1:numel(iters)-1
    i  
plot(allselectedArrayAligned(:,iters(i):iters(i+1)-1))
pause()
end
close
%%
%%%%if i really wanted to i could make a new swall (or like swfraction) of
%%%%only the spikes of interest by indexing them from swall and pulling
%%%%those spikes out then use spike sorter again to get those spikes back 

 %% save vars with unique name based on date and time (and cell info)        Z:\Users\Phil\Data\SpikesData\indexingVars
 dateString = datestr(now,'yymmdd_HHMM');
 save(['Z:\Users\Phil\Data\SpikesData\indexingVars\indexingMat_' setToScaleOrAligned '_' mouseName '_' sessionName '_' cellnum '_' code '_' dateString '.mat'],...
    'indexingMat','spksToPlotGood' ,'spksToPlotBad', 'good', 'bad', 'TrialIndCount', 'preSpikeBaselineTime')
%%
badSize = sum(cell2mat(cellfun(@numel,bad, 'UniformOutput',false)))
goodSize = sum(cell2mat(cellfun(@numel,good, 'UniformOutput',false)))

allBadSpikes = [spksToPlotBad{:}];
allBadSpikesAligned =  allBadSpikes - repmat(sum(allBadSpikes(preSpikeBaselineTime,:))/numel(preSpikeBaselineTime),30,1);
allBadSpikesScale = (allBadSpikes - repmat(mean(allBadSpikes(1:rise,:)),30,1))./repmat(max(allBadSpikes)-mean(allBadSpikes(preSpikeBaselineTime,:)),30,1);

allGoodSpikes = [spksToPlotGood{:}];
allGoodSpikesAligned =  allGoodSpikes - repmat(sum(allGoodSpikes(preSpikeBaselineTime,:))/numel(preSpikeBaselineTime),30,1);
allGoodSpikesScale = (allGoodSpikes - repmat(mean(allGoodSpikes(preSpikeBaselineTime,:)),30,1))./repmat(max(allGoodSpikes)-mean(allGoodSpikes(preSpikeBaselineTime,:)),30,1);

%% look at bad spikes 
figure
badSize
numSpikesPerPlot = 100;
iters = 1:numSpikesPerPlot:badSize;
if iters(end)~=badSize
    iters(end+1)=badSize;
end
numel(iters)-1
for i = 1:numel(iters)-1
    i  
plot(allBadSpikesAligned(:,iters(i):iters(i+1)-1))
pause()
end
close
%% look at good spikes 
figure
goodSize
numSpikesPerPlot = 1000;
iters = 1:numSpikesPerPlot:goodSize;
if iters(end)~=goodSize
    iters(end+1)=goodSize;
end
numel(iters)-1
for i = 1:numel(iters)-1
    i  
plot(allGoodSpikesAligned(:,iters(i):iters(i+1)-1))
pause()
end
close
%% Finalize the sorting of good spikes also cross ref good vars here!!
%%%%% goodSpikes are literally all spikes use to matter before but dont
%%%%% need now, left in in case we want to cross correlate two sets of
%%%%% spikes for example if we do this for alinged spikes and scaled spikes
%%%%% and are careful to never leave spikes out but include (if needed)
%%%%% lik or noise that is recognized as a spike then we can use my program
%%%%% to output the good spikes at 'goodSpikes' ann then it spikes will be
%%%%% sorted really well. 


allgood = cell(size(swall)); %integrate both good(removed artifact fr om waveform) and goodSpikes (removed artifact based on beambreaktimes) 
for i=1:length(s.sweeps)
        allgood{i} = intersect(goodSpikes{i},good{i}');%goodSpikes is for lick removed and good is filtered by user
end
%%


spikes_trials = s.get_spike_times; % takes into account allgood to replace spikes
for i = 1:length(spikes_trials.spikesTrials)
    spikes_trials.spikesTrials{i}.spikeTimes = spikes_trials.spikesTrials{i}.spikeTimes(allgood{i});
    viewsp{i}=sw{i}.spikeWaveforms{2}(:,[allgood{i}]);
end

% % % % % %view all good waveforms
% % % % % spcat=[viewsp{:}];
% % % % % finalsw=(spcat - repmat(mean(spcat(1:rise,:)),30,1))./repmat(max(spcat)-mean(spcat(1:rise,:)),30,1);
% % % % % figure;plot(finalsw)
% % % % % title('Filtered')

%SAVE THE INDEXING MAT SO YOU DONT HAVE TO SORT SPIKES AGAIN JUST HAVE TO
%LOAD IT NOTE THIS DOES NOT SAVE MULTIPLE SO YOU HAVE TO DO THAT TO SAVE IF
%YOU USE SCALE AND ALIGNED SPIKES CROSS REFERENCE TO MAKE YOUR SPIKE ARRAY 
save(['Z:\Users\Phil\Data\SpikesData\indexingVars\indexingMat_final_'  mouseName '_' sessionName '_' cellnum '_' code '.mat'],...
    'indexingMat','spksToPlotGood' ,'spksToPlotBad', 'good', 'bad', 'preSpikeBaselineTime')

%%

%%########%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%%%%#################%%%%%
%%%%%%%%%%%%%%%%WHISKER ARRAY BUILDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%..........................................................%%%%%%%


d= (['Z:\Data\Video\' videoloc filesep mouseName filesep sessionName filesep])
load(['Z:\Users\Phil\Data\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'])

cd(d)
filelist=dir([d '*.measurements'])

dirTrialNums=[];
%trialNums=[100:241];  % enter which trial nums to process 
%includef=cell(length(trialNums),1);
includef=cell(length(filelist),1);

% Assign the trial numbers to existing .measurements files in the directory
% NOTE : This assumes that the .measurements files have a four digit number
% with leading zeros corresponding to trial number in string positions 29:32
% of the file name.  These index numbers may need to be changed to match up 
% to the numerical code of the trial number. 



for i=1:length(filelist);
    dirTrialNums(i)=str2double(strtok(filelist(i).name(15:18),'.')); % extract out the trial number from each measurements file present in directory
end
trialNums = dirTrialNums;
for i = 1: length(filelist)
    includef{i} = filelist(i).name(1:end-13);
end
%% Optional section for cross correlating behavior and video trials, if you didn't pay attention to trial numbers
vv = nan(max(dirTrialNums),1);

for i = 1:length(dirTrialNums) 
    
    
    if ~isempty(find(dirTrialNums == i,1)); %create matrix showing trial location
        fidx = find(dirTrialNums == i); %show index in matrix where trial is 

        tmp = load([filelist(fidx).name(1:end-13) '.bar']);
         vv(i) = tmp(1,2);
    else
    end
    
end

gngThreshold = nanmean(vv) % for continuous pole postion with equal width go/nogo ranges, mean vv = the transition between go/nogo.
vv2 = vv >= gngThreshold; % threshold for pole position on 1 side or the other
vv3 = vv < gngThreshold;
vvDiff = vv2-vv3;
figure
plot([vvDiff],'.')
hold on
plot(vv3,'ro')

%%  Build behavior number vector
bv = zeros(max(b.trialNums),1);

bv(b.trialNums) = b.trialTypes*2-1; %1=GO and -1=NOGO and 0=NAN
[c, lags] = xcorr(bv,vvDiff); %correlation between bv (behavioral GO and NOGO) and vv2 (video GO and NOGO trials)
[mc, mx] = max(c); %find max c value and max x value
lag_shift = lags(mx) %find lag in between both bv and vv2

for i = 1: length(filelist)
    includef{i} = filelist(i).name(1:end-13);
end

% self inputted values since tossed first 132 trials and know lag shift
% lag_shift = -1

trialNums = dirTrialNums+lag_shift;  % correct the trial numbers

hold on
plot(bv,'go')

% restrict processing to trials...
startTrial = 1;
endTrial = 99999999;
includef = includef(trialNums >= startTrial & trialNums <=endTrial);
trialNums =  trialNums(trialNums >= startTrial & trialNums <=endTrial);

%% Step 1 - Run without the mask for the first 8 trials to
% determine what MASK you need
% 
% Whisker.makeAllDirectory_WhiskerTrial(d,0,...'mask', [900 700 200; 50 62 0],...
%     'trial_nums',trialNums(1:8),'include_files',includef(1:8),...
%     'barRadius',10,'faceSideInImage', 'top', 'framePeriodInSec',.001,...
%     'imagePixelDimsXY',[460 270],'pxPerMm',33,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','leftward')
% 
% Whisker.makeAllDirectory_WhiskerSignalTrial(d,'polyRoiInPix',[132-33 132+33])%,'follicleExtrapDistInPix',80);
% Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'r_in_mm',2,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);
% wl = Whisker.WhiskerTrialLiteArray(d);
% 
% tid = 0; % Set trajectory ID to view
% Whisker.view_WhiskerTrialLiteArray(wl,tid)

%% Step 2 - make your mask

trialnum = '10';
vidName = ['Z:\Data\Video\' videoloc filesep mouseName filesep sessionName filesep mouseName 'x' sessionName '-' trialnum '.mp4'];
vidFile = VideoReader(vidName);
figure(500)
imagesc(readFrame(vidFile))
if exist('maskPoints')
   hold on 
   plot(maskPoints(1,:),maskPoints(2,:))
end
points =ginput; %many clicks makes a nice fit!!!!
points = round(points);

maskPointsA = rot90(points(:,1),3);
maskPointsB = rot90(points(:,2),3);
maskPoints = [maskPointsA; maskPointsB];

close
%% Step 3 - Run with mask! for the first 8 trials and see how it looks

Whisker.makeAllDirectory_WhiskerTrial(d,0,'mask', maskPoints,...
    'trial_nums',trialNums(1:8),'include_files',includef(1:8),...
    'barRadius',10,'faceSideInImage', 'top', 'framePeriodInSec',.001,...
    'imagePixelDimsXY',[560 347],'pxPerMm',33,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','leftward')

Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef(1:8),'polyRoiInPix',[99-33 99+33],'follicleExtrapDistInPix',33);
Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef(1:8),'r_in_mm',3,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);
wl = Whisker.WhiskerTrialLiteArray(d);

tid = 0; % Set trajectory ID to view
Whisker.view_WhiskerTrialLiteArray(wl,tid)  % Open the viewer



%%%%%%%%%% check mask - load .WST file first
%%% this trial number is what you will see for the mask  
trialnum ='105';
WSTName = ['Z:\Data\Video\' videoloc filesep mouseName filesep sessionName filesep mouseName 'x' sessionName '-' trialnum '_' 'WST' '.mat']
load(WSTName)%load file based on trial above to test mask 
tp = [1.0 3.5];
figure;ws.plot_fitted_whisker_time_projection(0,'k',tp)
hold on; ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
hold on; ws.plot_mask(0,'g',tp)



%% save all variable to run multiple later                      Z:\Users\Phil\Data\BUILDER\SavedBuilderVars

save(['Z:\Users\Phil\Data\BUILDER\SavedBuilderVars\buildingVars_' cellNumberForProject])
%MAKE UNIQUE TO PREVENT ACCIDENTAL OVERWRITE CAN DELETE WHEN DONE
dateString = datestr(now,'yymmdd_HHMM');
save(['C:\Users\maire\Desktop\savedVarsFromBuilder\buildingVars_' cellNumberForProject '_' dateString])



%% Step 3 - run everything
% select matching files
%tmp = cellfun(@(x)str2num(x(15:end)),includef);
%incf_idx = find(tmp>= 12 & tmp <=85);
tic 

Whisker.makeAllDirectory_WhiskerTrial(d,0,'mask',maskPoints,...
    'trial_nums',trialNums,'include_files',includef,...
    'barRadius',10,'faceSideInImage', 'top', 'framePeriodInSec',.001,...
    'imagePixelDimsXY',[560 346],'pxPerMm',33,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','leftward')

Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef,'polyRoiInPix',[99-33 99+33],'follicleExtrapDistInPix',33);
Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef,'r_in_mm',3,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);;
wl = Whisker.WhiskerTrialLiteArray(d);
save([d mouseName sessionName '-WTLIA.mat'],'wl');

tid = 0; % Set trajectory ID to view
Whisker.view_WhiskerTrialLiteArray(wl,tid)
timeSpent = toc
hours = timeSpent/(60*60)

try
%%%%%%%SAVE THE VARIABLES TO LOAD ON YOUR PC AND MAKE AND SAVE THE T
%%%%%%%VARIABLE (WITH ALL THE PROJECT DETAILS ONLY WORK ON MY PC) 
save(['Z:\Users\Phil\Data\BUILDER\savedFinishedVars\finishedVars_' cellNumberForProject])
%MAKE UNIQUE TO PREVENT ACCIDENTAL OVERWRITE CAN DELETE WHEN DONE

%using date from original save so we can match them. should load from
%'buildingVars'
save(['C:\Users\maire\Desktop\savedVarsFromBuilder\finishedVars_' cellNumberForProject '_' dateString])

catch 
    save(['C:\Users\Public\finishedVars_' cellNumberForProject '_' dateString])
    display('had to save to local')
end

%%########%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%%%%#################%%%%%
%%%%%%%%%%%%%%%%WHISKER ARRAY BUILDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%..........................................................%%%%%%%




%%

%%%%%%%%%%%%%%%%%%FIRST BUILD WHISKER ARRAY AND THEN Load whisker data for session:
load([NASletter ':\Data\Video\' videoloc filesep mouseName filesep sessionName filesep mouseName sessionName '-WTLIA.mat'],'wl')

%%%%%%%%%%%%%%%%%%%%% Save everything
spikes_trials = s.get_spike_times; % takes into account allgood to replace spikes
for i = 1:length(spikes_trials.spikesTrials)
    spikes_trials.spikesTrials{i}.spikeTimes = spikes_trials.spikesTrials{i}.spikeTimes(allgood{i});
end

save(['Z:\Users\Phil\Data\SpikesData\SweepArrays\sweepArray_' mouseName '_' sessionName '_' cellnum '_' code '_' cellNumberForProject '.mat'],'s')
save(['Z:\Users\Phil\Data\SpikesData\indexingVars\indexingMat_'  mouseName '_' sessionName '_' cellnum '_' code '_' cellNumberForProject '.mat'],...
    'indexingMat','spksToPlotGood' ,'spksToPlotBad', 'good', 'bad', 'TrialIndCount')
save(['Z:\Users\Phil\Data\SpikesData\SpikeArrays\spikes_trials_'  mouseName '_' sessionName '_' cellnum '_' code '_' cellNumberForProject '.mat'],'spikes_trials')


T = LCA.TrialArray(b,spikes_trials,wl);
% to make this work projectDetails you must make that a property in
% +LCA/@TrialArray/TrialArray.m file and then copy over the @projectDetails
% folder into +LCA folder. Won't affect anyone elses stuff-PSM 
T.projectDetails = LCA.projectDetails;


T.projectDetails.cellNumberForProject = cellNumberForProject;
T.projectDetails.projectName = projectName;

T.projectDetails.category1 = category1;
T.projectDetails.category2 = category2;
T.projectDetails.category3 = category3;
T.projectDetails.category4 = category4;
T.projectDetails.category5 = category5;
T.projectDetails.cellNotes = cellNotes;
        
T.projectDetails.prePerformanceRegion = prePerformanceRegion;
T.projectDetails.goodPerformanceRegion = goodPerformanceRegion;
T.projectDetails.postPerformanceRegion = postPerformanceRegion; 

T.projectDetails.spikesNormalRegion = spikesNormalRegion; 
T.projectDetails.spikesFastBrokeInCellRegion = spikesFastBrokeInCellRegion;



T.whiskerTrialTimeOffset = whiskerTrialTimeOffset;
T.depth = depth;
T.recordingLocation = recordingLocation;



for i = 1:length(T.trials)
    T.trials{i}.spikesTrial.spikeTimes = T.trials{i}.spikesTrial.spikeTimes(:);
end

save(['Z:\Users\Phil\Data\TrialArrayBuilders\trial_array_'  mouseName '_' sessionName '_' cellnum '_' code '_' cellNumberForProject '.mat'],'T')
cd('Z:\Users\Phil\Data\TrialArrayBuilders\');

figure;T.plot_spike_raster(0,'BehavTrialNum')
T


save(['Z:\Users\Phil\Data\Characterization\trial_array_'  cellNumberForProject '.mat'],'T')
%date and time to make unique files in case of accidental overwirte can
%delete these when not needed anymore
save(['Z:\Users\Phil\Data\Characterization\trial_array_'  cellNumberForProject '_' datestr(datetime) '.mat'],'T') 









