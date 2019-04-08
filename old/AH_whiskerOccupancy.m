%% JC edt

L=horzcat(L3,L5b);


%%
thetaAtTouch = {};
thetaAtTouchW = {};
goTrials = {};
goTimes = {};
noGoTrials = {};
noGoTimes = {};
whisking = {};
touchOn = {};
theta={};
phase={};
setpoint={};
touchpert={};
edges = -40:5:60;
radedges = -pi:pi/8:pi;
gotouchpert={};
Nogotouchpert={};
goLick={};
NogoLick={};
golickpert={};
Nogolickpert={};
golickprob={};
Nogolickprob={};
motorpos={};
trialcorr={};
sorted={};
sortedBy={};
pcorrect={};
normpcorrect={};

for i = 1:length(L)
    whisking{i} = L{i}.S_ctk(3,:,:)>2.5;
    nonwhisking{i} = L{i}.S_ctk(3,:,:)<1.25;
    %goTrials{i} = sort([find(L{i}.trialTypeMat(1,:)) find(L{i}.trialTypeMat(2,:))]);
    goTrials{i}=find(L{i}.meta.trialType==1);
    %noGoTrials{i} = sort([find(L{i}.trialTypeMat(3,:)) find(L{i}.trialTypeMat(4,:))]);
    noGoTrials{i}=find(L{i}.meta.trialType==0);
    
    motorpos{i}=L{i}.meta.motorPosition;
    trialcorr{i}=L{i}.meta.trialCorrect;
    
    goTimes{i} = zeros(size(whisking{i}));
    goTimes{i}(:,:,goTrials{i}) = 1;
    
    noGoTimes{i} = zeros(size(whisking{i}));
    noGoTimes{i}(:,:,noGoTrials{i}) = 1;
    
    touchOn{i} = ~isnan(L{i}.S_ctk(9,:,:)) + ~isnan(L{i}.S_ctk(12,:,:));
    
    
    %     touchGo = ~isnan(L{i}.S_ctk(9,:,goTrials{i})) + ~isnan(L{i}.S_ctk(12,:,goTrials));
    %     touchNoGo = ~isnan(L{i}.S_ctk(9,:,noGoTrials)) + ~isnan(L{i}.S_ctk(12,:,noGoTrials));
    %
    phase{i}= L{i}.S_ctk(5,:,:);
    theta{i} = L{i}.S_ctk(1,:,:);
    setpoint{i} = L{i}.S_ctk(4,:,:);
    
    %     thetaGo = L{i}.S_ctk(1,:,goTrials);
    %     thetaNoGo = L{i}.S_ctk(1,:,noGoTrials);

    preLickMask{i} = zeros(size(whisking{i}));
    for j = 1:L{i}.k
        firstLick = min([find(L{i}.S_ctk(16,:,j)==1,1,'first') find(L{i}.S_ctk(15,:,j)==1,1,'last') L{i}.t]);
        preLickMask{i}(1,1:firstLick,j) = 1;
    end
    %finding probability lick based on prelick touches
    for k = 1:L{i}.k
        touchpert{i}(1,1,k)=numel(find(touchOn{i}(:,:,k).*preLickMask{i}(:,:,k)));%finds prelick touches in each trial
        firstlick{i}(1,1,k)=numel(find(L{i}.S_ctk(16,:,k)==1,1,'first'));%finds time point of first lick
    end
    
    for d= 1:length(goTrials{i}) %sort trials by go v nogo
        gotouchpert{i}(d)=touchpert{i}(1,1,goTrials{i}(d)); %adds row under touches per trial that list trial type
        goLick{i}(d)=firstlick{i}(1,1,goTrials{i}(d)); %
    end
    
    for c= 1:length(noGoTrials{i}) %sort trials by go v nogo
        Nogotouchpert{i}(c)=touchpert{i}(1,1,noGoTrials{i}(c));
        NogoLick{i}(c)=firstlick{i}(1,1,noGoTrials{i}(c));
    end
    
    gotouchpert{i}(2,:)=goLick{i};
    for z=0:max(gotouchpert{i}(1,:))%cycling through each touch order
        [~,col]=find(gotouchpert{i}(1,:)==z);
        golickpert{i}(z+1)= sum(gotouchpert{i}(2,col))/numel(col);
        golickprob{z+1}(1,i)=golickpert{i}(z+1);%finding mean p(lick) based on touch number
    end
    Nogotouchpert{i}(2,:)=NogoLick{i};
    for y=0:max(Nogotouchpert{i}(1,:));
        [~,col]=find(Nogotouchpert{i}(1,:)==y);
        Nogolickpert{i}(y+1)= sum(Nogotouchpert{i}(2,col))/numel(col);
        Nogolickprob{y+1}(1,i)=Nogolickpert{i}(y+1);
    end
    
    % motor position vs accuracy 
    mpvsa{i}=vertcat(motorpos{i},trialcorr{i});   
    [sorted{i} sortedBy{i} binBounds]=binslin(mpvsa{i}(1,:),mpvsa{i}(2,:),'equalE',19,0,180000);
    pcorrect{i}=cellfun(@mean,sorted{i});
    normpcorrect{i}=pcorrect{i};
    normpcorrect{i}(isnan(normpcorrect{i}))=[];
    
    
    thetaAtTouchGo{i} = theta{i}(find(touchOn{i}.*preLickMask{i}.*goTimes{i}));
    thetaAtTouchGoW{i} = theta{i}(find(touchOn{i}.*preLickMask{i}.*goTimes{i}.*whisking{i}));
    thetaAtTouchNoGo{i} = theta{i}(find(touchOn{i}.*preLickMask{i}.*noGoTimes{i}));
    thetaAtTouchNoGoW{i} = theta{i}(find(touchOn{i}.*preLickMask{i}.*noGoTimes{i}.*whisking{i}));
    
    phaseAtTouchGo{i} = phase{i}(find(touchOn{i}.*preLickMask{i}.*goTimes{i}));
    phaseAtTouchNoGo{i} = phase{i}(find(touchOn{i}.*preLickMask{i}.*noGoTimes{i}));
    
    h_tgo{i} = histc(thetaAtTouchGo{i},edges);
    h_tgow{i} = histc(thetaAtTouchGoW{i},edges);
    
    h_tnogo{i} = histc(thetaAtTouchNoGo{i},edges);
    h_tnogow{i} = histc(thetaAtTouchNoGoW{i},edges);
    
    h_pgo{i} = histc(phaseAtTouchGo{i},radedges);
    h_pnogo{i} = histc(phaseAtTouchNoGo{i},radedges);
    
    gngTRatio{i} = (sum(h_tgow{i})-sum(h_tnogow{i})) / sum(sum(h_tgow{i})+sum(h_tnogow{i}));
    goRange{i} = [-5 5]+mean(thetaAtTouchGoW{i});
    noGoRange{i}= [-5 5]+mean(thetaAtTouchNoGoW{i});
    
    cropThetaGo{i} = theta{i}(find(whisking{i}.*preLickMask{i}.*goTimes{i}));
    cropThetaGo_noW{i} = theta{i}(find(whisking{i}.*preLickMask{i}.*goTimes{i}));

    cropThetaNoGo{i} = theta{i}(find(whisking{i}.*preLickMask{i}.*noGoTimes{i}));
    occupancyGo{i} = sum(cropThetaGo{i} > goRange{i}(1) & cropThetaGo{i} < goRange{i}(2))/numel(cropThetaGo{i});
    occupancyNoGo{i}= sum(cropThetaNoGo{i} > noGoRange{i}(1) & cropThetaNoGo{i} < noGoRange{i}(2))/numel(cropThetaNoGo{i});


    setpointAtTouchGo{i} = setpoint{i}(find(touchOn{i}.*goTimes{i}.*whisking{i}));
    setpointAtTouchNoGo{i} = setpoint{i}(find(touchOn{i}.*noGoTimes{i}.*whisking{i}));
end
goodSessions = find(cellfun(@(x)x.k,L)>100) %finds good sessions based on number of trials

%% plot theta and phase at each point of touch. p selects which cell
%cluster cells based on same mouse 
%ichi:AH0171 =1,2,3,4,5,15
%ni:AH0122 =6,16,17,18,19,20
%yon:AH0124 =7,8,9,10
%go:AH0112 =11,12,22
%roku:AH0287 =13,14,25
%nana:AH0283 =23,24
mouseNum='Grand Mean Continuous Whisker Occupancy (n=19)';
unonum=[6 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25];%fill in w/ cell nums you want to plot 

t_gouno=sum(horzcat(h_tgo{unonum}),2);
t_nogouno=sum(horzcat(h_tnogo{unonum}),2);
p_gouno=sum(horzcat(h_pgo{unonum}),2);
p_nogouno=sum(horzcat(h_pnogo{unonum}),2);

figure;
subplot(2,1,1);hold on
bar(edges(1:end)-.5,t_gouno,'b')
bar(edges(1:end)-.5,t_nogouno,'r')
xlabel('Theta at Touch')
ylabel('# of touches')
title(mouseNum)

subplot(2,1,2);hold on
bar(radedges(1:end),p_gouno,'b')
bar(radedges(1:end),p_nogouno,'r')
xlabel('Phase at Touch')
ylabel('# of touches')
set(gca,'xlim',[-pi pi],'xtick',pi*[-1:.5:1],'xticklabel',{'-pi','-pi/2',0,'pi/2','pi'})

whiskingTouchRatio = cellfun(@nume,thetaAtTouchW)./cellfun(@numel,thetaAtTouch);

%% Touch v Lick Prob
%nogolickprob/golickprob: each cell = touch number, each col = lick prob for a
%trial
accuracy=zeros(1,length(L));
for u=1:length(L)
accuracy(u)=sum(L{u}.meta.trialCorrect)/numel(L{u}.meta.trialCorrect);
end
a=find(accuracy>=.7);
discrim = {golickpert{a}};
figure(1);

totalgolprob=cellfun(@nanmean,golickprob); %all trials binned together regardless of discrim accuracy 
totalnogolprob=cellfun(@nanmean,Nogolickprob);

plot([0:length(totalgolprob)-1],totalgolprob,'b')
ylim([0 1])
hold on;plot([0:length(totalnogolprob)-1],totalnogolprob,'r')
title('p(lick) vs pre lick touches')

% look at single trial
o=5
hold on
plot([0:length(golickpert{o})-1],golickpert{o},'b')
xlim([0 length(golickpert{o})])
ylim([0 1])
hold on;plot([0:length(Nogolickpert{o})-1],Nogolickpert{o},'r')
title('p(lick) vs pre lick touches')

%% % accuracy by pole position binned @ 10000 microunits
%sorted holds all information of trial correct vs incorrect in a cell array
%pcorrect takes avg of correct answers at each position within each trial
%fill in w/ cell nums you want to plot 

motoredges=[0:1:17];
dosnum=[8 9 10]; %select which sessions you would like to cluster together
%8 9 10 good AH0124
%
dosmean=nanmean(horzcat(pcorrect{dosnum}),2);
plot(motoredges,dosmean);
split=(L{dosnum(1)}.meta.goPosition(1)-L{dosnum(1)}.meta.nogoPosition(1))/20000; %finds midpoint in sessions
xlabel('Anterior -> Posterior')
ylabel('p(correct)')
title('p(correct) vs motor position')
hold on
plot([split split],[0 1],'--r')

%normalizing using normpcorrect to try to plot ALL trials regardless of mouse
% This pretty much removes all nans from unused motor positions. Not sure how good this is since the ranges on some
%trials vary. Need to find better way to normalize (subtract by max and min
%of pole pos and then divide in two?)
motoredges2=[0:1:7];
tresnum=[3 4 5 15];
tresmean=nanmean(horzcat(normpcorrect{tresnum}),2);
plot(motoredges2,tresmean);
xlabel('Anterior -> Posterior')
ylabel('p(correct)')
title('p(correct) vs motor position')
hold on
plot([3.5 3.5],[0 1],'--r')


%%
for i = 1:52
    closeness(i) = noGoRange{i}(1)+5-mean(setpointAtTouchNoGo{i})
end


% figure(2)
% bar(edges-.5,sum([h_tw{:}],2),'k')
