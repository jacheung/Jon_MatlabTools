function [varAligned, prevarAligned] = assist_varAtTouchFORCE(array,range,variable)
% This function will find the indices of touch, the spikes around that
% touch given by the RANGE (ie -25:50ms) you provide. Lastly it'll find the
% theta, phase, amplitude, setpoint, max kappa, and pre touch velocity.

% Simple organization tool for visualizing parameters at touch onset

% First 6 columns will be values for the variables 
% 1) THETA 2) AMP 3) SETPOINT 4) PHASE 5) MAX KAPPA 6) PRE TOUCH VELOCITY 
% Last columns will be the spikes around your given window 

%INPUT: vector with time points you want to view (ie [-25:50])



window = numel(range);

%set touch ranges
touchOnIdx = [find(array.S_ctk(9,:,:)==1); find(array.S_ctk(12,:,:)==1)];
touchOffIdx = [find(array.S_ctk(10,:,:)==1); find(array.S_ctk(13,:,:)==1)];
spikes = squeeze(array.R_ntk);
touchOnIdx = touchOnIdx(touchOnIdx<(numel(spikes)-range(end)));
touchOffIdx = touchOffIdx(1:length(touchOnIdx));
touchOnIdx = sort(touchOnIdx,1);
touchOffIdx = sort(touchOffIdx,1);

%align variables all around ALL touches (range before/after touch)
varAligned=zeros(numel(touchOnIdx),window);
var= squeeze(array.S_ctk(variable,:,:));
for i = 1:size(varAligned,1)
    varAligned(i,:) = var(touchOnIdx(i)+range);
    if variable == 1%% this is made so we get a change in theta from touch
        varAligned(i,:)=varAligned(i,:)-var(touchOnIdx(i));%added negative sign at front (i think this is right)
   % elseif variable == 6
        %varAligned(i,:)=varAligned(i,:)-var(touchOnIdx(i));
    end
end

%align variables all around PREDECISION touches (range before/after touch)
[prelixtouches ,~, ~, ~ ,~ ,~] = assist_predecisionVar(array);
prevarAligned=cell(1,length(prelixtouches));
for d = 1:length(prelixtouches)
    theta=array.S_ctk(1,:,d);
    var=array.S_ctk(variable,:,d);
    prevarAligned{d}=zeros(numel(prelixtouches{d}),window+1);
    for f = 1:numel(prelixtouches{d})
        prevarAligned{d}(f,2:end)=var(prelixtouches{d}(f)+range);
        prevarAligned{d}(f,1)=theta(prelixtouches{d}(f));
        if variable==1
            prevarAligned{d}(f,2:end)=prevarAligned{d}(f,2:end)-var(prelixtouches{d}(f));
        end
    end
end



%% 
%this variable is used to build scatters like Pammer 2013 fig 8
% radialvar=cell(1,length(touchrange)); 
% 
% for i=1:length(touchrange)
%     radialvar{i}(1,:)= array.S_ctk(1,touchrange(i,1):touchrange(i,2));
%     radialvar{i}(2,:)= array.S_ctk(6,touchrange(i,1):touchrange(i,2)); 
%     radialvar{i}(3,:)= array.S_ctk(8,touchrange(i,1):touchrange(i,2));
% end
% 
% radialvar2=cell2mat(radialvar);
% scatter(radialvar2(1,:),radialvar2(2,:));%plotting all theta during touch against Moadj
%scatter(radialvar2(1,:),radialvar2(3,:)); %plotting all theta during touch against Faxialadj
%% THIS CAN BE IGNORED - USEFUL WHEN WE USED THIS TO CORRELATE WITH SPIKES
%variables to be measured 
vartouch = zeros(numel(touchOnIdx),6);%theta 1, amp 2, setpoint 3, phase 4, max curvature 5, vel 6
for i = 1:length(touchOnIdx)
    vartouch(i,1)=array.S_ctk(1,touchOnIdx(i)); %theta at touch
    vartouch(i,2)=array.S_ctk(3,touchOnIdx(i)); %amp at at touch
    vartouch(i,3)=array.S_ctk(4,touchOnIdx(i)); %setpoint at touch
    vartouch(i,4)=array.S_ctk(5,touchOnIdx(i)); %phase at touch
        kwin=array.S_ctk(6,touchOnIdx(i):touchOffIdx(i)); %get values in touch window
    vartouch(i,5)= max(abs(kwin)); %find idx of max kappa within each touch window, neg or pos
                     %use idx to pull max kappa
    vartouch(i,6)=mean(array.S_ctk(2,touchOnIdx(i)-4:touchOnIdx(i)-1)); %finds mean of velocity (-4:-1ms) before touch
end

%mash variables and spikes into a single matrix
varAligned=horzcat(vartouch(:,1),varAligned);
prevarAligned=prevarAligned(~cellfun('isempty',prevarAligned)); %remove empty cells  
prevarAligned=cell2mat(prevarAligned');
%%


