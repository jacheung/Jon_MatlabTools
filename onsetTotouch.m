function [ttimes] = onsetTotouch(array)
%Function used to find time of first touch after pole onset. Outputs ttimes
%which is a cell array organized by trial type (H,FA,CR,MISS). Two column
%matrix with column 1 as time to touch and column 2 as trial number. 

%Variables for classification
nogos =find(array.meta.trialType ==0);
gos = find(array.meta.trialType ==1);
corr = find(array.meta.trialCorrect ==1);
incorr = find(array.meta.trialCorrect == 0);
hits = intersect(gos,corr);
FA = intersect(nogos,incorr);
CR = intersect(nogos,corr);
miss = intersect(gos,incorr);
ponset = round((array.meta.poleOnset(1)).*1000,0);

%trial outcome filters
type = {[hits],[FA],[CR],[miss]};
ttimes = cell(1,4);
colors = {'b','g','r','k'};
for b = 1:length(type)
     if ~isempty(type{b}) %check to make sure there are trials
         trs = type{b};
         for d = 1:length(type{b})
             if ~isempty(find(squeeze(array.S_ctk(9,:,trs(d)))==1))
             ttimes{b}(d,:) = [find(squeeze(array.S_ctk(9,:,trs(d)))==1)- ponset' trs(d)];
             else 
             ttimes{b}(d,:) = [nan trs(d)];
             end
%              ttimes{b}(ttimes{b}(:,1) == 0,:) = []; %elim non touch trials
         end       
     end
end
       