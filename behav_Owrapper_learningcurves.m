function [] = learningCurves(behavDataLocation)

clear all; close all; clc
%
 %groups = [{'ContLearningCurves'} {'SemiLearningCurves'}];
groups = [{'ContLearningCurves'}];
for i=1:length(groups)

%wrap all behavioral data files
[hist, info, trim] = behav_bdatawrapper('Z:\Users\Jon\DATA\Behavior\ContLearningCurves');
clear O
uniqueMouseNums = nan(length(hist),1); 
for d = 1:length(hist)
    output = behav_processBehavior_v2(info{d},hist{d},trim{d});
    output.trialResults(:,8) = d;
    output.mouseName = info{d}.SavingSection_MouseName;
    output.trialResultsNames = {'trial number', 'bar at go', 'bar at nogo', 'trials correct', 'errors','Dprime','use trials','session','trialstart','motor pos','most ant go','most ant nogo','numtouches'};
    O(d) = output;
    uniqueMouseNums(d) = info{d}.SavingSection_MouseName;
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
behav_learningcurves(mousel,4,200); %plot learning curves for each mouse
end

