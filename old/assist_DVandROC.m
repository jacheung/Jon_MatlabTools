function [listDVclose,listDVfar,ROC,AUC] = assist_DVandROC(farPSTH,closePSTH)
%INPUTS: 
% far = vector listing all far trial numbers
% close = vector listing all close trial numbers 
% farPSTH = mean number of spikes at each sample
    % ie: average number of spikes generated within touch window 
    % ie: average number of spikes generated into predecision window
% closePSTH = mean number of spikes at each sample

%OUTPUTS: 
%this function will output a list of decision variables at close/far which
%you can plot in a histogram to see distribution. Will also output ROC
%which is a two column matrix that contains pairs of probabilities for
%deciphering one PSTH from another. Plot it to see. Lastly, will output the
%AUC for the ROC. 


%% DECISION VARIABLE/ROC
%equation for "Decision Variable" for hit = ti (mean (hit - evaluating hit) - mean (CR)) 
%equation for "Decision Variable" for CR = ti (mean (hit) - mean (CR-evaluating CR)) 
closePSTH=closePSTH(~isnan(closePSTH)); %elim all NaNs
farPSTH=farPSTH(~isnan(farPSTH));
DVclose=zeros(1,length(closePSTH));
DVfar=zeros(1,length(farPSTH));

%Build Decision Variables for both vars. Meaning aligning on axis from 0 as
%most sim to farthest, least sim, for selected PSTH to mean PSTH
for j=1:length(DVclose)
    tmpclose=closePSTH;
    tmpclose(j)=[];%all close trials - evaluating trial
    DVclose(j)=dot(closePSTH(j),(mean(tmpclose))-(mean(farPSTH)));
end
for j=1:length(DVfar)
    tmpfar=farPSTH;
    tmpfar(j)= [];%all far trials - evaluating trial
    DVfar(j)=dot(farPSTH(j),(mean(closePSTH))-(mean(tmpfar)));
end
DVclose=DVclose*1000;%converting to spks/s^2 instead of spks/ms^2
DVfar=DVfar*1000;%converting to spks/s^2 instead of spks/ms^2

listDVclose=zeros(2,length(DVclose));
listDVfar=zeros(2,length(DVfar));
listDVclose=unique(DVclose);
for j=1:length(listDVclose)
    listDVclose(2,j)=numel(find(DVclose==listDVclose(1,j)));
end
listDVfar=unique(DVfar);
for j=1:length(listDVfar)
    listDVfar(2,j)=numel(find(DVfar==listDVfar(1,j)));
end

%Run ROC analysis between both variables 
[ROCclose,ROCfar]=assist_ROCcurves(DVclose,DVfar);
ROC = [ROCfar' ROCclose'];
ROC = unique(ROC,'rows');%minimize the amount of pairs to shrink file size

 
 if sum(ROC)>0
    AUC = abs(trapz(ROCfar,ROCclose));
 else
     AUC = 0;
 end

 

