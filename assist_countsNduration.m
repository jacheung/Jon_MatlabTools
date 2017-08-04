function [touchcounts, gosorted, nogosorted, godurfinal, nogodurfinal] = assist_countsNduration(array,prelixGo,prelixNoGo,duration);
go=find(array.meta.trialType==1);
nogo=find(array.meta.trialType==0);

%% Sorting touch count vs trial correct
gocount=[cellfun(@numel,prelixGo);array.meta.trialCorrect(go)==1];
nogocount=[cellfun(@numel,prelixNoGo);array.meta.trialCorrect(nogo)==1];
[tmpgosorted,~,~]=binslinTMP(gocount(1,:)',gocount','equalE',21,0,20);
[tmpnogosorted,~,~]=binslinTMP(nogocount(1,:)',nogocount','equalE',21,0,20);
gosorted=cell2mat(cellfun(@(x) mean(x,1),tmpgosorted,'uniformoutput',0));
nogosorted=cell2mat(cellfun(@(x) mean(x,1),tmpnogosorted,'uniformoutput',0));
% gosorted=gosorted(~isnan(gosorted(:, 1)), :);
% nogosorted=nogosorted(~isnan(nogosorted(:, 1)), :);

touchcounts=[histc(cellfun(@numel,prelixGo),0:20)/numel(prelixGo);histc(cellfun(@numel,prelixNoGo),0:20)/numel(prelixNoGo)]';
figure(123);clf;subplot(1,2,1)
bar(0:20,touchcounts);
set(gca,'xlim',[-1 20])
hold on;plot(gosorted(:,1),gosorted(:,2),'-o');
hold on;plot(nogosorted(:,1),1-nogosorted(:,2),'r-o');
xlabel('Number of contacts per trial');ylabel('P(lick) or Proportion of Trials')
title('Touch Count vs P(lick)')

% saved for future when we want to reanalyze touch count x P(lick) based on
% hit,miss,CR,FA
% hitcount=prelixGo(array.meta.trialCorrect(go)==1);
% misscount=prelixGo(array.meta.trialCorrect(go)==0);
% CRcount=prelixNoGo(array.meta.trialCorrect(nogo)==1);
% FAcount=prelixNoGo(array.meta.trialCorrect(nogo)==0);

%% Sorting touch duration vs trial correct
goduration=cell2mat(duration(go));nogoduration=cell2mat(duration(nogo));

[godursorted,~,gobounds]=binslinTMP(goduration(1,:)',goduration','equalE',26,1,101); %binning durations in 2-3ms bins
goboundsx=mean([gobounds;[gobounds(2:end) 0]]); %midway b/t points for bounds 
godurcounts=[cellfun('size',godursorted,1)'; goboundsx(1:end-1)];
goduracc=cell2mat(cellfun(@(x) mean(x,1),godursorted,'uniformoutput',0));%mean accuracy for each duration bin
godurfinal=[godurcounts;goduracc(:,2)'];%row 1) # samples, 2) duration, 3)mean accuracy
godurfinal(1,:)=godurfinal(1,:)/sum(goduration(1,:)>0);%dividing duration by number of touches (use sum b/c 0s count as no touches)

[nogodursorted,~,nogobounds]=binslinTMP(nogoduration(1,:)',nogoduration','equalE',26,1,101);
nogoboundsx=mean([nogobounds;[nogobounds(2:end) 0]]);
nogodurcounts=[cellfun('size',nogodursorted,1)'; nogoboundsx(1:end-1)];
nogoduracc=cell2mat(cellfun(@(x) mean(x,1),nogodursorted,'uniformoutput',0));
nogodurfinal=[nogodurcounts;nogoduracc(:,2)'];
nogodurfinal(1,:)=nogodurfinal(1,:)/sum(nogoduration(1,:)>0);%dividing duration by number of touches (use sum b/c 0s count as no touches)

figure(123);subplot(1,2,2);
bar(godurfinal(2,:),[godurfinal(1,:)' nogodurfinal(1,:)'])
hold on; plot(godurfinal(2,:),godurfinal(3,:),'-o')
hold on; plot(nogodurfinal(2,:),1-nogodurfinal(3,:),'r-o')
xlabel('Duration of Touches (ms)');ylabel('P(lick) or Proportion of touches per trial')
title('Duration of Touches vs P(Lick)')
set(gca,'xlim',[0 100])


