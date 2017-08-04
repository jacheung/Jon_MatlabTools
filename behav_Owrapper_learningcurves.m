clear all; close all; clc
%
groups = [{'ContLearningCurves'} {'SemiLearningCurves'}];

for i=1:length(groups)

%wrap all behavioral data files
[hist, info, trim] = behav_bdatawrapper(['Z:\Users\Jon\DATA\Behavior\' groups{i}]);
%info{2}.MotorsSection_no_pole_position_ant=60000;
clear O
for d = 1:length(hist)
    output = behav_processBehavior_v2(info{d},hist{d},trim{d});
    output.trialResults(:,8) = d;
    output.mouseName = info{d}.SavingSection_MouseName;
    output.trialResultsNames = {'trial number', 'bar at go', 'bar at nogo', 'trials correct', 'errors','Dprime','use trials','session','trialstart','motor pos','most ant go','most ant nogo','numtouches'};
    O(d) = output;
end


%Data sheet with all recorded mouse numbers

gname = strmatch(groups{i},'ContLearningCurves');
if gname == 1
        [~, txt, ~] = xlsread('RecordedMice','F14:F20');
else
        [~, txt, ~] = xlsread('RecordedMice','G14:G23');
end


mousel = cell(1,numel(txt));

%sort o structure array by mouse numbers and plot dprime/accuracy per mouse
for j=1:length(txt)
    tmp=strfind({O.mouseName},txt{j}); 
    logarray = ~cellfun(@isempty,tmp);
    mouseIdx=find(logarray);
    goodo=O(mouseIdx);
    mousel{j} =behav_smoothBehaviorWindow(goodo,100);
end

 %PLOTS
behav_learningcurves(mousel,4,100); %plot learning curves for each mouse
end