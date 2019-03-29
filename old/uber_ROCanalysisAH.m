


L.meta.trialType = L.trialTypeMat(1,:,:)+L.trialTypeMat(2,:,:);
L.meta.trialCorrect =  L.trialTypeMat(1,:,:)+L.trialTypeMat(3,:,:);
AUClist=zeros(length(U),2);

 %% need just touches prior to lick and find average spks/touch there
 disp(rec)
tmp=find(L.S_ctk(9,:,:)==1)/L.t;
    d = tmp-floor(tmp); %keep only decimals
    avail=round(min(d)*L.t); %pole available time in ms
 touchIdx = cell(1,L.k);
 lickIdx = cell(1,L.k);
 spikes = cell(touchIdx);
 meanspikes = cell(touchIdx);
 firstlick = zeros(1,L.k);
 
 for i = 1:L.k
    allLicks = find(L.S_ctk(16,:,i)==1);
    lickIdx{i} = allLicks(allLicks>avail);
    if isempty(lickIdx{i}) %index of all licks after pole availability 
        lickIdx{i}=0;
    end
    firstlick(i) = lickIdx{i}(1);
    alltouches = horzcat(find(L.S_ctk(9,:,i)==1),find(L.S_ctk(12,:,i)==1));
    alltouches = alltouches(alltouches<L.t-100);
    %touchIdx{i} = alltouches(alltouches<lickIdx{i}(1)); %all touches before first lick
    touchIdx{i} = alltouches; %all touches after avail to end 
    spikes{i}=zeros((numel(touchIdx{i})),45);
    meanspikes{i}=zeros((numel(touchIdx{i})),1);
 end

 for k = 1:length(touchIdx)%within each trial 
     for m = 1:length(touchIdx{k})%for each touch before lick
         spikes{k}(m,:) = L.R_ntk(:,touchIdx{k}(m)+[6:50],k); %look at spikes 6:50ms after touch
         meanspikes{k}(m) = mean(spikes{k}(m,:),2);
     end
 end

%useful variables here are 
%spikes = cell array where each cell = 1 trial; rows in each cell =
%touches and columns correspond to 6:50ms after touch.
%meanspikes = cell array where each cell = 1 trial; rows in each cell = 
%touches and column = mean spike of the PSTH in spikes. 
 
%% Only looking at PSTH within window from pole availability to average RXN time
meanrxntime = round(mean(firstlick(find(firstlick>0))));
availtorxn = [avail:meanrxntime];
availtoend = [avail:L.t];

winrange=availtorxn;

winspikes=cell(1,L.k);
for i = 1:length(winspikes)
    winspikes{i} = L.R_ntk(:,winrange,i);
end

%useful variables here are
% availtorxn or availtoend which decides window with which to build
% decision variables. Choose by setting winrange to either one.

nogo = find(L.meta.trialType==0); %NogoTrials
go = find(L.meta.trialType==1); %GoTrials
nolick = [find(L.meta.trialType == 1 & L.meta.trialCorrect == 0) find(L.meta.trialType == 0 & L.meta.trialCorrect == 1)];
lick = [find(L.meta.trialType == 1 & L.meta.trialCorrect == 1) find(L.meta.trialType == 0 & L.meta.trialCorrect == 0)];


goPSTH=zeros(length(go),numel(winrange));
nogoPSTH = zeros(length(nogo),numel(winrange));
lickPSTH=zeros(length(lick),numel(winrange));
nolickPSTH=zeros(length(nolick),numel(winrange));

for i = 1:length(go)
    goPSTH(i,:) = winspikes{go(i)};
end
for i = 1:length(nogo)
    nogoPSTH(i,:) = winspikes{nogo(i)};
end
for i = 1:length(lick)
    lickPSTH(i,:) = winspikes{lick(i)};
end
for i = 1:length(nolick)
    nolickPSTH(i,:) = winspikes{nolick(i)};
end

meangoPSTH=mean(goPSTH,2); %find mean of PSTH for each trial 
meannogoPSTH=mean(nogoPSTH,2); 
meanlickPSTH=mean(lickPSTH,2); 
meannolickPSTH=mean(nolickPSTH,2); 


[listDVgo,listDVnogo,ROCtot,AUCtot] = assist_DVandROC(meangoPSTH,meannogoPSTH);
[listDVlick,listDVnolick,ROCdtot,AUCdtot] = assist_DVandROC(meanlickPSTH,meannolickPSTH);



%ASIDE:plotting hit vs CRPSTH of touches before decision 
%plot(smooth(sum(hitPSTH)));
%hold on; plot(smooth(sum(CRPSTH)),'r');

% %plotting Decision variables: Could use some work (R2015B much better 
% see code at bottom)
% figure(9);clf;
% bar(listDVhit(1,:),listDVhit(2,:));
% hold on; bar(listDVCR(1,:),listDVCR(2,:),'r');
% ylabel('# of trials')
% xlabel('"Decision Variable" (spks/s^2)')
% title('DVs for window')

figure(11);clf;scatter(ROCtot(:,1),ROCtot(:,2))
hold on;plot(ROCtot(:,1),ROCtot(:,2))
hold on;plot([0,1],[0,1],':b');
 ylabel('P(DVgo>crit)')
 xlabel('P(DVnogo>crit)')
 title('ROC Window (Pole Available:Mean Rxn Time)')
 if sum(ROCtot)>0
 AUCtot=abs(trapz(ROCtot(:,1),ROCtot(:,2)));
  else
     AUCtot = 0;
 end
 text(.7,.2,['AUC =' num2str(AUCtot)],'FontSize',12,'Color','black')
 
 
 
 figure(12);clf;scatter(ROCdtot(:,1),ROCdtot(:,2))
hold on;plot(ROCdtot(:,1),ROCdtot(:,2))
hold on;plot([0,1],[0,1],':b');
 ylabel('P(DVlick>crit)')
 xlabel('P(DVnolick>crit)')
 title('ROC Window (Pole Available:Mean Rxn Time)')
 if sum(ROCdtot)>0
 AUCdtot=abs(trapz(ROCdtot(:,1),ROCdtot(:,2)));
  else
     AUCdtot = 0;
 end
 text(.7,.2,['AUC =' num2str(AUCdtot)],'FontSize',12,'Color','black')
 
%% DV for TOUCH COUNT  NOT WORKING - MATRIX TOO BIG FOR DVANDROC
 
nogo = find(L.meta.trialType==0); %NogoTrials
go = find(L.meta.trialType==1); %GoTrials
nolick = [find(L.meta.trialType == 1 & L.meta.trialCorrect == 0) find(L.meta.trialType == 0 & L.meta.trialCorrect == 1)];
lick = [find(L.meta.trialType == 1 & L.meta.trialCorrect == 1) find(L.meta.trialType == 0 & L.meta.trialCorrect == 0)];

goPSTH = zeros(length(go),1);
nogoPSTH=zeros(length(nogo),1);
lickPSTH = zeros(length(lick),1);
nolickPSTH = zeros(length(nolick),1);
for i = 1:length(go)
    goPSTH(i)=numel(meanspikes{go(i)});
end   
for i = 1:length(nogo)
    nogoPSTH(i)= numel(meanspikes{nogo(i)});
end
for i = 1:length(lick)
    lickPSTH(i)=numel(meanspikes{lick(i)});
end   
for i = 1:length(nolick)
    nolickPSTH(i)=numel(meanspikes{nolick(i)});
end   

[DVgo,DVnogo,ROCgng,AUCgng] = assist_DVandROC(goPSTH,nogoPSTH);
[DVlick,DVnolick,ROClicks,AUClicks] = assist_DVandROC(lickPSTH,nolickPSTH); 

%% DV FOR MEAN OF TOUCHES PSTH WITHIN TRIAL
nogo = find(L.meta.trialType==0); %NogoTrials
go = find(L.meta.trialType==1); %GoTrials
closePSTH=[];
farPSTH=[];

for i = 1:length(go)
    closePSTH=vertcat(closePSTH,mean(meanspikes{go(i)}));
end
   
for i = 1:length(nogo)
    farPSTH=vertcat(farPSTH,mean(meanspikes{nogo(i)}));
end

[listDVclose,listDVfar,ROC,AUC] = assist_DVandROC(farPSTH,closePSTH);
% PLOT DECISION VARIABLES
% figure(10);clf;
% bar(listDVclose(1,:),listDVclose(2,:));
% hold on; bar(listDVfar(1,:),listDVfar(2,:),'r');
% ylabel('# of trials')
% xlabel('"Decision Variable" (spks/s^2)')
% title('DVs for each touch treated individually')

% PLOT ROC CURVE
figure(12);clf;scatter(ROC(:,1),ROC(:,2))
hold on;plot(ROC(:,1),ROC(:,2))
hold on;plot([0,1],[0,1],':b');
title('ROC Mean Touch PSTH/Trial')
 ylabel('P(DVFAR>crit)')
 xlabel('P(DVCLOSE>crit)')
 text(.7,.2,['AUC =' num2str(AUC)],'FontSize',12,'Color','black')
 
%% DV FOR ALL TOUCH PSTH TREATED AS INDIVIDUALS: optional if you wanted to treat each touch individually 
% far = find(L.meta.trialType==0); %NogoTrials
% close = find(L.meta.trialType==1); %GoTrials
% closePSTH=[];
% farPSTH=[];
% for i = 1:length(close)
%     closePSTH=vertcat(closePSTH,meanspikes{close(i)});
% end
%    
% for i = 1:length(far)
%     farPSTH=vertcat(farPSTH,meanspikes{far(i)});
% end
%% Putting AUC in a list for downstream plotting 
if numel(cell2mat(touchIdx(nogo)))<10 %if less than 10 PSTHs total at farPSTH then call that cell NaN
    AUClist(rec,:)=NaN;
else
 AUClist(rec,1)=AUCtot;
 AUClist(rec,2)=AUC;
end

% 
%  %% Plot scatter for AUC(rxn time) vs AUC (touches)
% AUClist(~any(~isnan(AUClist), 2),:)=[];
%  scatter(AUClist(:,1),AUClist(:,2),'k');
%   scatter(L5bAUClist(:,1),L5bAUClist(:,2),'r');
%  hold on;plot([0,1],[0,1],':k');
%  xlabel('AUC Reaction Window')
%  ylabel('AUC Touches')
%  title(['AUC for L3 vs L5b'])
%  axis('square')
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 2])
% %% Plot both ROCtot and ROC in the same window
% 
% figure(15);clf
% scatter(ROCtot(:,1),ROCtot(:,2),'r')
% hold on; plot(ROCtot(:,1),ROCtot(:,2),'r');
% hold on; scatter(ROC(:,1),ROC(:,2))
% hold on;plot(ROC(:,1),ROC(:,2))
% hold on;plot([0,1],[0,1],':k');
%  text(.7,.2,['AUC =' num2str(AUCtot)],'FontSize',12,'Color','red')
%  text(.7,.25,['AUC =' num2str(AUC)],'FontSize',12,'Color','blue')
% axis('square')
% title('ROC')
%  ylabel('P(DVFAR>crit)')
%  xlabel('P(DVCLOSE>crit)')
%  
 %% R2015b optimization for histogram function
 
 firstvar = listDVclose
 secondvar = listDVfar
     lb=floor(min(horzcat(firstvar,secondvar))); %find lower bound for ranges of criterion
    ub=ceil(max(horzcat(firstvar,secondvar)));
 
 figure(16);clf;histogram(firstvar,lb:.05:ub,'FaceColor','b','FaceAlpha',.5)
hold on; histogram(secondvar,lb:.05:ub,'FaceColor','r','FaceAlpha',.5)
ylabel('# of trials')
xlabel('"Decision Variable" (spks/s^2)')
title('DVs for mean touch spikes')

%      lb=floor(min(horzcat(DVclose,DVfar))); %find lower bound for ranges of criterion
%     ub=ceil(max(horzcat(DVclose,DVfar)));
%     edge=lb:.5:ub;
%  figure(17);clf;histogram(DVfar,edge,'FaceColor','r')
%  hold on;histogram(DVclose,edge,'FaceColor','b')
% ylabel('# of trials')
% xlabel('"Decision Variable" (spks/s^2)')
% title('DVs for Mean Touch Spikes')


hold on; histogram(DVfar,edge)






