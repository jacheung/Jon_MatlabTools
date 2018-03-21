%Function for finding the touch times from whisk onset to touch. Couple
%specifications is...
% 1) only looks within window of pole availability to lick
% 2) only looks at protraction touches
% 3) toss touch times >100ms
%Written by JC 180315 along with findMaxMinProtraction



function [timetotouch,trialnums] = uber_touchTiming(array)

%P struct for finding peak and trough protraction values
% [P] = findMaxMinProtraction(array,'avail2lick');
[P] = findMaxMinProtraction(array);
[objmask]= assist_touchmasks(array);
mask = objmask.availtolick;

%Variables for classification
nogos =find(array.meta.trialType ==0);
gos = find(array.meta.trialType ==1);
corr = find(array.meta.trialCorrect ==1);
incorr = find(array.meta.trialCorrect == 0);
hits = intersect(gos,corr);
FA = intersect(nogos,incorr);
CR = intersect(nogos,corr);
miss = intersect(gos,incorr);

%trial outcome filters
type = {[hits],[FA],[CR],[miss]};
colors = {'b','g','r','k'};
for b = 1:length(type)
    if ~isempty(type{b}) %check to make sure there are trials
        tmpmask=mask(:,type{b});
        phases = squeeze(array.S_ctk(5,:,type{b}));
        touchIdx = [find(array.S_ctk(9,:,type{b})==1);find(array.S_ctk(12,:,type{b})==1)];
        %finding only touches within masked window (avail2lick)
        touchIdx= intersect(touchIdx,find(tmpmask==1));
        
        %remove retraction touches
         touchIdx(phases(touchIdx)>0)=[];

        if isempty(touchIdx)
            timetotouch{b} = [];
            trialnums{b} = [];
        else
            toss=[]; %tossing out touches where we cant find trough
            neartrough = zeros(length(touchIdx),0);
            for k = 1:length(touchIdx)
                if isempty(find(P.troughIdx<touchIdx(k),1,'last'))
                    neartrough(k) = 1;
                    toss = [toss k];
                else
                    neartrough(k) = find(P.troughIdx<touchIdx(k),1,'last');
                end
            end
            
            tmpIdx = P.troughIdx(neartrough);
            TTidx = [tmpIdx touchIdx];
            TTidx(toss,:) = [];
            [~,keep] = unique(TTidx(:,1)); %removing doubling ups since sometimes we dont capture the troughIdx
           
            
            TTidx=TTidx(keep,:);
            touchtimes = TTidx(:,2)-TTidx(:,1);
            
            %eliminate any touches > 2 std from the median values
%           oob = find(touchtimes>median(touchtimes) + 2*std(touchtimes));
            oob = find(touchtimes>100);
            touchtimes(oob) = []; 
            tnumIdx = ceil(TTidx(:,2)/4000);
            tnumIdx(oob) = [];
            
            
            timetotouch{b} = touchtimes;
            trialnums{b} = type{b}(tnumIdx)';
        end
%                  figure(50);hold on;histogram(touchtimes,'binedges',[0:2: 100],'facecolor',colors{b},'normalization','probability')
    else
        timetotouch{b} = [];
        trialnums{b} = [];
    end
    
end

%% Test plots for identifying touch and protraction start
% thetas = squeeze(array.S_ctk(1,:,:));
% 
% xmin=1
% xmax=length(TTidx)
% n=1
% trial=round(xmin+rand(1,n)*(xmax-xmin))
% 
% % trial = 75
% ranges = (trial*4000)+1:(trial*4000)+4000;
% 
% figure(58);clf;
% plot(thetas(ranges(1):ranges(end)));
% 
% plotIdx = find(TTidx(:,2)<ranges(end) & TTidx(:,2)>ranges(1));
% 
% 
% plottouches=TTidx(plotIdx,2);
% plotonset=TTidx(plotIdx,1);
% hold on; scatter(plottouches-ranges(1),thetas(plottouches),'r')
% hold on; scatter(plotonset-ranges(1),thetas(plotonset),'g')