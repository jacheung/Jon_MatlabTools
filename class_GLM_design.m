%% load spikes data and convert to milliseconds
FFspkssparse=full(c.timeSeriesArrayHash.value{2}.valueMatrix(2,:));
dist = find(FFspkssparse==1);
binFFspkssparse=zeros(1,length(FFspkssparse)/10);
test2 = round(dist/10);
binFFspkssparse(test2)=1;

%% first touch onset
touchonsetIdx = find(c.timeSeriesArrayHash.value{1}.valueMatrix(7,:)==1);

trial=c.timeSeriesArrayHash.value{2}.trial;
trialIds=c.trialIds;
trialtimeidx=zeros(1,length(trialIds));%in milliseconds time points where each trial begins

for i=1:length(trialIds)
    tnum=trialIds(i);
    trialtimeidx(i)=find(trial==tnum,1);
end
trialtimeidx=trialtimeidx/10;
firsttouch=cell(numel(trialtimeidx),1);

for i=1:length(trialtimeidx)
    templogic =  find(touchonsetIdx>trialtimeidx(i) & touchonsetIdx<trialtimeidx(i)+4500,1);
    firsttouch{i}=touchonsetIdx(templogic);
end
firsttouch=firsttouch(~cellfun('isempty',firsttouch));
firsttouch=cell2mat(firsttouch);

spikesAlignedfirst = zeros(numel(firsttouch),100);
for i = 1:size(spikesAlignedfirst,1)
    spikesAlignedfirst(i,:) = binFFspkssparse(firsttouch(i)+[-50:49]);
end

firsttouchidx=zeros(1,length(trial));
firsttouchidx(firsttouch)=1; %index in time where first touch occurs

%% Normal dist fit
%need to feed it time of when spikes occurred
inds=find(spikesAlignedfirst==1); %finds index of all spikes
[~, col] = ind2sub(size(spikesAlignedfirst),inds); %gets column where spikes occurred
col(col<50)=NaN;%set all spikes before 0 to NaN and not be counted in distribution
dist=fitdist(col,'Normal');
dist.mu=dist.mu-50; %to account for index starting at 0 instead of -50

figure;
bar(-50:49,sum(spikesAlignedfirst)/numel(firsttouch),'k');
title('First Touch Onset')
xlabel('Time from touch (ms)')
ylabel('Spikes per touch')
hold on
plot([-50:49],pdf(dist,[-50:49]),'LineWidth',2)
%% First touch spike design matrix
ton=find(firsttouchidx==1);
epoch=cell(1,length(ton));
for i=1:length(ton)
    epoch{i}=binFFspkssparse(ton(i)-50:ton(i)+49); %using index of touch onset, pull spikes data out
end
epoch=(cell2mat(epoch))';

%tonfull = repmat([zeros(50,1);1;zeros(49,1)], [length(ton), 1]); %pretty much abov but in one line
%% First touch stim design matrix (-50ms:50ms)
hypparam=dist.mu

tft=repmat([(-50:49)'], [length(ton), 1]); %time from touch (50ms before to 50ms after)
newtft=exp(-(tft-hypparam).^2/(2*dist.sigma^2));%tft into gaussian distrubtion w/ mean at 10ms from touch

%% spike history design matrix
spkidx=find(binFFspkssparse==1); %find all periods of spikes
spkhistdesign=cell(1,length(spkidx));
for i =1:length(spkidx)
    spkhistdesign{i}=binFFspkssparse(spkidx(i)-5:spkidx(i)-1);%create vector of 4ms before spike
end
spkhistdesign=(cell2mat(spkhistdesign)'); %spks surrounding spkidx
spkhiststim=repmat([(-5:-1)'],length(spkidx),1); %stim for time scale of spkhistdesign
%% GLM Build
[gnew,dev,stats] = glmfit(newtft,epoch,'poisson');%glm fit w/ modeled Gaussian dist.
figure; plot(-50:49,exp(gnew(1)+gnew(2)*(newtft(1:100)))/0.001)
title('GLM with Gaussian Fit')
xlabel('Time from Touch Onset (ms)')
ylabel('Firing Rate (spks/s)')
[glin,dev,stats] = glmfit(tft,epoch,'poisson'); %glmfit using only linear fit
figure; plot(-50:49,exp(glin(1)+glin(2)*(tft(1:100)))/0.001)
title('GLM with Linear Fit')
xlabel('Time from Touch Onset (ms)')
ylabel('Firing Rate (spks/s)')
[gzero,dev,stats] = glmfit(ones(size(epoch)),epoch,'poisson','constant','off');%constant is normally vector of 1s to count for baseline frequency. Here we turn it off and add our own vector of 1s
[ghist,dev,stats] = glmfit(spkhiststim,spkhistdesign,'poisson');
%% AIC for model selection. Lower AIC = better model 
%use exp(a0)*Delta to calculate FR/ms as the results are in log. exp gets rid of
delt=.001;%for ms

ratefunc=(gnew(1)+gnew(2)*newtft);
params=length(gnew);
llfunc=sum(epoch.*(ratefunc)-delt*exp(ratefunc)); %log likelihood function
aic=(-2*llfunc)+(2*params); % aic = -2loglikelihod + 2x(params)

linratefunc=(glin(1)+glin(2)*tft);
linparams=length(glin);
linllfunc=sum(epoch.*(linratefunc)-delt*exp(linratefunc)); %log likelihood function
linaic=(-2*linllfunc)+(2*linparams); % aic = -2loglikelihod + 2x(params)

zratefunc=(gzero(1));
zparams=length(gzero);
zllfunc=sum(epoch.*zratefunc-delt*exp(zratefunc));
zaic=(-2*zllfunc)+(2*zparams);

allaic=[aic linaic zaic]
%% KS plots for computed GLM w/ hyperparams
glmratefunc=exp(gnew(1)+gnew(2)*newtft);

%glmnewfr=glmratefunc/delt; %spks/second

%1) computing time rescaled theorem zk (find index of spikes in each trial,
%sum up all points in rate function correlating to those points)
u=cell(1,length(ton));
for i=1:length(ton)
    base=1:100:length(epoch);
    tmpidx=find(epoch(base(i):base(i)+99)==1)+base(i)-1;
    empt=cell(1,length(tmpidx)-1);
    for j=1:length(tmpidx)-1 
        empt{j}=sum(glmratefunc(tmpidx(j):tmpidx(j+1)));%finds time between each spike in trial
    end
    u{i}=cell2mat(empt);
end
 
%2) compute 1-e^(-zk)
tmp2=1-exp(-(cell2mat(u)));

%3) sort tmp2=y axis
uk=sort(tmp2);

%4) plot against bk = (k-1/2/L)..... L = #uk
bk=([1:length(uk)]-.5)/length(uk);
figure; plot(bk,uk);
xlim([0 1])

%5) add confidence bounds (1.36/root(L))
hold on
bnd = 1.36/sqrt(length(tmp2));
plot([0 1-bnd],[bnd 1]);
plot([bnd 1],[0 1-bnd]);
title('Goodness of Fit for GLM')

%% cross validation training w/ first half, testing w/ 2nd half 
select=(round(length(ton)/2)/length(ton));%gives % of 1st touches rounding up in touches odd number

trainnewtft=newtft(1:select*length(newtft));%first half stim  for test 
trainepoch=epoch(1:select*length(newtft));%spikes
testnewtft=newtft(select*length(newtft)+1:length(newtft));
testepoch=epoch(select*length(newtft)+1:length(newtft));

[gnewtrain,dev,stats] = glmfit(trainnewtft,trainepoch,'poisson'); %using last half epoch to train
gnewratefunc=exp(gnewtrain(1)+gnewtrain(2)*trainnewtft);%rate function
%time rescaling 
ux=cell(1,length(testepoch)/100);
for i=1:length(ux)
    base=1:100:length(testepoch);
    tmpidx=find(testepoch(base(i):base(i)+99)==1)+base(i)-1;
    empt=cell(1,length(tmpidx)-1);
    for j=1:length(tmpidx)-1 
        empt{j}=sum(gnewratefunc(tmpidx(j):tmpidx(j+1)));%sums ratefunction = [p(firing) at each ms] between each spike in trial
    end
    ux{i}=cell2mat(empt);
end
%computing 1-e^(-zk)
tmpx=1-exp(-(cell2mat(ux)));

%sort
ukx=sort(tmpx);

%4) plot against bk = (k-1/2/L)..... L = #uk
bkx=([1:length(tmpx)]-.5)/length(tmpx);
figure; plot(bkx,ukx);
xlim([0 1])

%5) add confidence bounds (1.36/root(L))
hold on
bnd = 1.36/sqrt(length(tmpx));
plot([0 1-bnd],[bnd 1]);
plot([bnd 1],[0 1-bnd]);
title('Cross Validation of GLM')

%% Garbage: Fake spike data
% spkwindow = rand(15,1)>.3; %70% chance of generating 1's in a 15ms window
% spksmall=cell(1,length(ton));
% for s=1:length(spksmall);
%     spks=zeros(55,1);
%     spks=insertrows(spks,(rand(15,1)>.3),[55]);
%     spksmall{s}=insertrows(spks,zeros(30,1),[70])';
% end
% spkfull=(cell2mat(spksmall)');
