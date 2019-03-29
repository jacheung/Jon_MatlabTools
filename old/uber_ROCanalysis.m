if strcmp(U{1}.meta.layer,'L5b')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','I23:J50');
    disp(U{1}.meta.layer);
    layer = 'L5b';
elseif strcmp(U{1}.meta.layer,'L3')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','B25:C44');
    disp(U{1}.meta.layer);
    layer = 'L3';
elseif strcmp(U{1}.meta.layer,'L4')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','O23:P28');
    disp(U{1}.meta.layer);
    layer = 'L4';
elseif strcmp(U{1}.meta.layer,'L3Out')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','B77:C82');
    disp(U{1}.meta.layer);
    layer = 'L3Out';
else strcmp(U{1}.meta.layer,'L5bOut')
    [ndata, txt, alldata] =xlsread('CellsConversionChart','I75:J80');
    disp(U{1}.meta.layer);
    layer = 'L5bOut';
end
txt=txt(~isnan(ndata));
ndata=ndata(~isnan(ndata));

AUClist=zeros(length(ndata),2);

for rec = 1:length(ndata)
 %% need just touches prior to lick and find average spks/touch there
 disp(rec)
tmp=find(U{rec}.S_ctk(9,:,:)==1)/4000;
    d = tmp-floor(tmp); %keep only decimals
    avail=round(min(d)*4000); %pole available time in ms
 touchIdx = cell(1,U{rec}.k);
 lickIdx = cell(1,U{rec}.k);
 spikes = cell(touchIdx);
 meanspikes = cell(touchIdx);
 firstlick = zeros(1,U{rec}.k);
 
 for i = 1:U{rec}.k
    allLicks = find(U{rec}.S_ctk(16,:,i)==1);
    lickIdx{i} = allLicks(allLicks>avail);
    if isempty(lickIdx{i}) %index of all licks after pole availability 
        lickIdx{i}=0;
    end
    firstlick(i) = lickIdx{i}(1);
    alltouches = horzcat(find(U{rec}.S_ctk(9,:,i)==1),find(U{rec}.S_ctk(12,:,i)==1));
    alltouches = alltouches(alltouches<3950);
    %touchIdx{i} = alltouches(alltouches<lickIdx{i}(1)); %all touches before first lick
    touchIdx{i} = alltouches; %all touches after avail to end 
    spikes{i}=zeros((numel(touchIdx{i})),45);
    meanspikes{i}=zeros((numel(touchIdx{i})),1);
 end

 for k = 1:length(touchIdx)%within each trial 
     for m = 1:length(touchIdx{k})%for each touch before lick
         spikes{k}(m,:) = U{rec}.R_ntk(:,touchIdx{k}(m)+[6:50],k); %look at spikes 6:50ms after touch
         meanspikes{k}(m) = mean(spikes{k}(m,:),2);
     end
 end

%useful variables here are 
%spikes = cell array where each cell = 1 trial; rows in each cell =
%touches and columns correspond to 6:50ms after touch.
%meanspikes = cell array where each cell = 1 trial; rows in each cell = 
%touches and column = mean spike of the PSTH in spikes. 
 
%% Only looking within window from pole availability to average RXN time
meanrxntime = round(mean(firstlick(find(firstlick>0))));
availtorxn = [avail:meanrxntime];
availtoend = [avail:U{rec}.t];

winrange=availtorxn;

winspikes=cell(1,U{rec}.k);
for i = 1:length(winspikes);
    winspikes{i} = U{rec}.R_ntk(:,winrange,i);
end

%useful variables here are
% availtorxn or availtoend which decides window with which to build
% decision variables. Choose by setting winrange to either one.

far = find(U{rec}.meta.trialType==0); %NogoTrials
close = find(U{rec}.meta.trialType==1); %GoTrials

hitPSTH=zeros(length(close),numel(winrange));
CRPSTH = zeros(length(far),numel(winrange));

for i = 1:length(close)
    hitPSTH(i,:) = winspikes{close(i)};
end

for i = 1:length(far)
    CRPSTH(i,:) = winspikes{far(i)};
end
meanhitPSTH=mean(hitPSTH,2); %find mean of PSTH for each trial 
meanCRPSTH=mean(CRPSTH,2); 

[listDVhit,listDVCR,ROCtot,AUCtot] = assist_DVandROC(meanhitPSTH,meanCRPSTH);

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
 ylabel('P(DVFAR>crit)')
 xlabel('P(DVCLOSE>crit)')
 title('ROC Window (Pole Available:Mean Rxn Time)')
 if sum(ROCtot)>0
 AUCtot=abs(trapz(ROCtot(:,1),ROCtot(:,2)));
  else
     AUCtot = 0;
 end
 text(.7,.2,['AUC =' num2str(AUCtot)],'FontSize',12,'Color','black')
 
%% DV FOR MEAN OF TOUCHES PSTH WITHIN TRIAL
far = find(U{rec}.meta.trialType==0); %NogoTrials
close = find(U{rec}.meta.trialType==1); %GoTrials
closePSTH=[];
farPSTH=[];

for i = 1:length(close)
    closePSTH=vertcat(closePSTH,mean(meanspikes{close(i)}));
end
   
for i = 1:length(far)
    farPSTH=vertcat(farPSTH,mean(meanspikes{far(i)}));
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
% far = find(U{rec}.meta.trialType==0); %NogoTrials
% close = find(U{rec}.meta.trialType==1); %GoTrials
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
if numel(cell2mat(touchIdx(far)))<10 %if less than 10 PSTHs total at farPSTH then call that cell NaN
    AUClist(rec,:)=NaN;
else
 AUClist(rec,1)=AUCtot;
 AUClist(rec,2)=AUC;
end

end
 %% Plot scatter for AUC(rxn time) vs AUC (touches)
AUClist(~any(~isnan(AUClist), 2),:)=[];
 scatter(AUClist(:,1),AUClist(:,2),'k');
  scatter(L5bAUClist(:,1),L5bAUClist(:,2),'r');
 hold on;plot([0,1],[0,1],':k');
 xlabel('AUC Reaction Window')
 ylabel('AUC Touches')
 title(['AUC for L3 vs L5b'])
 axis('square')
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 2])
%% Plot both ROCtot and ROC in the same window

figure(15);clf
scatter(ROCtot(:,1),ROCtot(:,2),'r')
hold on; plot(ROCtot(:,1),ROCtot(:,2),'r');
hold on; scatter(ROC(:,1),ROC(:,2))
hold on;plot(ROC(:,1),ROC(:,2))
hold on;plot([0,1],[0,1],':k');
 text(.7,.2,['AUC =' num2str(AUCtot)],'FontSize',12,'Color','red')
 text(.7,.25,['AUC =' num2str(AUC)],'FontSize',12,'Color','blue')
axis('square')
title('ROC')
 ylabel('P(DVFAR>crit)')
 xlabel('P(DVCLOSE>crit)')
 
 %% R2015b optimization for histogram function
     lb=floor(min(horzcat(DVhit,DVCR))); %find lower bound for ranges of criterion
    ub=ceil(max(horzcat(DVhit,DVCR)));
 
 figure(16);clf;histogram(DVhit,lb:.05:ub,'FaceColor','b','FaceAlpha',.5)
hold on; histogram(DVCR,lb:.05:ub,'FaceColor','r','FaceAlpha',.5)
ylabel('# of trials')
xlabel('"Decision Variable" (spks/s^2)')
title('DVs for Window')
print(gcf, '-opengl','-depsc','197C_DVWindow')

     lb=floor(min(horzcat(DVclose,DVfar))); %find lower bound for ranges of criterion
    ub=ceil(max(horzcat(DVclose,DVfar)));
    edge=lb:.5:ub;
 figure(17);clf;histogram(DVfar,edge,'FaceColor','r')
 hold on;histogram(DVclose,edge,'FaceColor','b')
ylabel('# of trials')
xlabel('"Decision Variable" (spks/s^2)')
title('DVs for Mean Touch Spikes')


hold on; histogram(DVfar,edge)






