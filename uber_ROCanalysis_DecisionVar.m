   

%% find ALL touches and get spikes from 6:50ms after touch
touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];    
    spikes = squeeze(U{rec}.R_ntk);
    touchIdx = touchIdx(touchIdx<(numel(spikes)-151)); %elim last touches
    spikesAligned = zeros(numel(touchIdx),45);
    
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[6:50]);
    end
    
 %% need just touches prior to lick and find average spks/touch there
tmp=find(U{rec}.S_ctk(9,:,:)==1)/4000;
    d = tmp-floor(tmp); %keep only decimals
    avail=round(min(d)*4000); %pole available time in ms
 touchIdx = cell(1,U{rec}.k);
 lickIdx = cell(1,U{rec}.k);
 spikes = cell(touchIdx);
 firstlick = zeros(1,U{rec}.k);
 
 for i = 1:U{rec}.k
    allLicks = find(U{rec}.S_ctk(16,:,i)==1);
    lickIdx{i} = allLicks(allLicks>avail);
    if isempty(lickIdx{i}) %index of all licks after pole availability 
        lickIdx{i}=0;
    end
    firstlick(i) = lickIdx{i}(1);
    alltouches = horzcat(find(U{rec}.S_ctk(9,:,i)==1),find(U{rec}.S_ctk(12,:,i)==1));
    touchIdx{i} = alltouches(alltouches<lickIdx{i}(1)); %all touches before first lick
    spikes{i}=zeros((numel(touchIdx{i})),45);
 end

 for k = 1:length(touchIdx)%within each trial 
     for m = 1:length(touchIdx{k})%for each touch before lick
         spikes{k}(m,:)=U{rec}.R_ntk(:,touchIdx{k}(m)+[6:50],k); %look at spikes 6:50ms after touch
     end
 end

%% window spikes - built for 197C (looking at window b/t 1000ms-1500ms
meanrxntime = round(mean(firstlick(find(firstlick>0))));
winrange = [1000:meanrxntime];
winspikes=cell(1,U{rec}.k);
for i = 1:length(winspikes);
    winspikes{i} = U{rec}.R_ntk(:,winrange,i);
end


%% Decision variable (dot product similarity of single trial to HIT trials -
% CR trials)
%need mean PSTH hit and mean PSTH CR

CR = find(U{rec}.meta.trialType==0 & U{rec}.meta.trialCorrect==1); %CR trials
hit = find(U{rec}.meta.trialType==1 & U{rec}.meta.trialCorrect==1); %HIT trials

hitPSTH=zeros(length(hit),numel(winrange));
CRPSTH = zeros(length(CR),numel(winrange));

for i = 1:length(hit)
    hitPSTH(i,:) = winspikes{hit(i)};
end

for i = 1:length(CR)
    CRPSTH(i,:) = winspikes{CR(i)};
end
%ASIDE:plotting hit vs CRPSTH of touches before decision 
% plot(smooth(sum(hitPSTH)));
% hold on; plot(smooth(sum(CRPSTH)),'r');

meanhitPSTH=mean(hitPSTH,2); %find mean of PSTH for each trial 
meanCRPSTH=mean(CRPSTH,2); 
DVhit=zeros(1,length(hit));
DVCR=zeros(1,length(CR));

%equation for "Decision Variable" for hit = ti (mean (hit - evaluating hit) - mean (CR)) 
%equation for "Decision Variable" for CR = ti (mean (hit) - mean (CR-evaluating CR)) 

for j=1:length(DVhit)
    tmphit=find(hit~=hit(j));%all hit trials - evaluating trial
    DVhit(j)=dot(meanhitPSTH(j),(mean(meanhitPSTH(tmphit)))-(mean(meanCRPSTH)));
end

for j=1:length(DVCR)
    tmpCR=find(CR~=CR(j));%all hit trials - evaluating trial
    DVCR(j)=dot(meanCRPSTH(j),(mean(meanhitPSTH))-(mean(meanCRPSTH(tmpCR))));
end

DVhit=DVhit*1000;%converting to spks/s^2 instead of spks/ms^2
DVCR=DVCR*1000;%converting to spks/s^2 instead of spks/ms^2
listDVhit=zeros(2,length(DVhit));
listDVCR=zeros(2,length(DVCR));

listDVhit=unique(DVhit);
for j=1:length(listDVhit)
    listDVhit(2,j)=numel(find(DVhit==listDVhit(1,j)));
end

listDVCR=unique(DVCR);
for j=1:length(listDVCR)
    listDVCR(2,j)=numel(find(DVCR==listDVCR(1,j)));
end

%plotting Decision variables 
% figure(9);clf;
% bar(listDVhit(1,:),listDVhit(2,:));
% hold on; bar(listDVCR(1,:),listDVCR(2,:),'r');
% ylabel('# of trials')
% xlabel('"Decision Variable" (spks/s^2)')
%% ROC Curve
%need probability CR 

% mu = mean(DVhit);
% sigma=std(DVhit);
% x=[0:1:45];
% pd=lognpdf(x,mu,sigma)*250;
% hold on; plot(pd)

lb=floor(min(horzcat(DVhit,DVCR)));
ub=ceil(max(horzcat(DVhit,DVCR)));
crit=lb:.0001:ub;
ROCCR=zeros(1,length(crit));
ROChit=zeros(1,length(crit));


for c=1:length(crit)
    ROCCR(c)=sum(DVCR>crit(c))/numel(DVCR);
    ROChit(c)=sum(DVhit>crit(c))/numel(DVhit);
end

 figure(11);scatter(ROCCR,ROChit);
 hold on;plot([0,1],[0,1],':b');
 ylabel('P(DVhit>crit)')
 xlabel('P(DVCR>crit)')

%% DV for hit trials
DVhit=zeros(1,length(hit));


for j = 1:length(DVhit)
    tmphit=hit(find(hit~=hit(j)));%all hit trials - evaluating trial
    hitPSTH=zeros(length(tmphit),45);
    CRPSTH = zeros(length(CR),45);
    evalPSTH = sum(spikes{hit(j)});
    for i = 1:length(tmphit)
        hitPSTH(i,:) = sum(spikes{tmphit(i)});
    end
    for i = 1:length(CR)
        CRPSTH(i,:) = sum(spikes{CR(i)});
    end
    DVhit(j)=evalPSTH.*(sum(hitPSTH)-sum(CRPSTH));
end
%plotting hit vs CRPSTH of touches before decision 
plot(smooth(sum(hitPSTH)));
hold on; plot(smooth(sum(CRPSTH)),'r');


