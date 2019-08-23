%Function for finding the touch times from whisk onset to touch. Couple
%specifications is...
% 1) only looks within window of pole availability to lick
% 2) only looks at protraction touches
% 3) toss touch times >100ms
%Written by JC 180315 along with findMaxMinProtraction



function [timetotouch,trialnums,velocity,trough] = uber_touchTiming(array)

%P struct for finding peak and trough protraction values
% [P] = findMaxMinProtraction(array,'avail2lick');
[P] = findMaxMinProtraction(array,1);
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
        thetas = squeeze(array.S_ctk(1,:,type{b}));
        phases = squeeze(array.S_ctk(5,:,type{b}));
        touchIdx = [find(array.S_ctk(9,:,type{b})==1);find(array.S_ctk(12,:,type{b})==1)];
        %finding only touches within masked window (avail2lick)
        touchIdx= intersect(touchIdx,find(tmpmask==1));
        
        %remove retraction touches
         touchIdx(phases(touchIdx)>0)=[];

        if isempty(touchIdx)
            timetotouch{b} = [];
            trialnums{b} = [];
            velocity{b} = [];
            trough{b} = [];
        else
            toss=[]; %tossing out touches where we cant find trough
            neartrough = zeros(length(touchIdx),0);
            for k = 1:length(touchIdx)
                if isempty(find(P.troughidx<touchIdx(k),1,'last'))
                    neartrough(k) = 1;
                    toss = [toss k];
                else
                    neartrough(k) = find(P.troughidx<touchIdx(k),1,'last');
                end
            end
            
            tmpIdx = P.troughidx(neartrough);
            TTidx = [tmpIdx touchIdx];
            TTidx(toss,:) = [];
            [~,keep] = unique(TTidx(:,1)); %removing doubling ups since sometimes we dont capture the troughidx
           
            
            TTidx=TTidx(keep,:);
            touchtimes = TTidx(:,2)-TTidx(:,1);
            velotmp = thetas(TTidx(:,2))-thetas(TTidx(:,1))./touchtimes;
            troughtmp = thetas(TTidx(:,1));
            
            %eliminate any touches > 2 std from the median values
%           oob = find(touchtimes>median(touchtimes) + 2*std(touchtimes));
            oob = find(touchtimes>100);
            touchtimes(oob) = []; 
            tnumIdx = ceil(TTidx(:,2)/array.t);
            tnumIdx(oob) = [];
            velotmp(oob) = [];
            troughtmp(oob) = [];
            
            trough{b} = troughtmp; 
            velocity{b} = velotmp;
            timetotouch{b} = touchtimes;
            trialnums{b} = type{b}(tnumIdx)';
        end
%                  figure(50);hold on;histogram(touchtimes,'binedges',[0:2: 100],'facecolor',colors{b},'normalization','probability')
    else
        velocity{b} = [];
        trough{b} = [];
        timetotouch{b} = [];
        trialnums{b} = [];
    end
    
end

%% % Test plots for identifying touch and protraction start
thetas = squeeze(array.S_ctk(1,:,:));
phases = squeeze(array.S_ctk(5,:,:));
% touchIdx = [find(array.S_ctk(9,:,type{b})==1);find(array.S_ctk(12,:,type{b})==1)];
touchIdx = [find(array.S_ctk(9,:,:)==1);find(array.S_ctk(12,:,:)==1)];
touchIdx= intersect(touchIdx,find(mask==1));

proTouchFilt =phases(touchIdx)<0;
touchIdx = touchIdx(proTouchFilt);


alltroughs = P.troughidx;      
troughsIdx = intersect(alltroughs,find(mask==1));

xmin=1;
xmax=array.k;
trialIdx=round(xmin+rand(1,1)*(xmax-xmin));

selranges = (trialIdx*array.t)+1:(trialIdx*array.t)+array.t;
plotTouch = mod(intersect(touchIdx,selranges),array.t);
plottroughs = mod(intersect(troughsIdx,selranges),array.t);
plotthetas = thetas(selranges);
figure(20);clf
hold on; plot(plotthetas,'k')
hold on; scatter(plotTouch,plotthetas(plotTouch),'ro')
hold on; scatter(plottroughs,plotthetas(plottroughs),'go')


%% 
xmin=1;
xmax=max(unique(ceil(TTidx./4000)));
trialIdx=round(xmin+rand(1,1)*(xmax-xmin));
tnum = type{b}(trialIdx)

thetaranges = (tnum*array.t)+1:(tnum*array.t)+array.t;
selTheta = thetas(thetaranges);
selTouchIdx = find(touchIdx<thetaranges(end) & touchIdx>thetaranges(1));

idxranges = (trialIdx*array.t)+1:(trialIdx*array.t)+array.t;
plotIdx = find(TTidx(:,2)<idxranges(end) & TTidx(:,2)>idxranges(1));

plottouches=TTidx(plotIdx,2);
plottouchesALL = touchIdx(plotIdxtall)
plotonset=TTidx(plotIdx,1);
touchtp = mod(plottouchesALL,array.t);
onsetp = mod(plotonset,array.t);

figure(58);clf;
plot(selTheta,'k');
hold on; scatter(touchtp,selTheta(touchtp),'r')
hold on; scatter(onsetp,selTheta(onsetp),'g')