prob=cell(1,length(U));

for rec=1:length(U)
    gos=find(U{rec}.meta.trialType==1);nogos=find(U{rec}.meta.trialType==0);
        gotmp=U{rec}.meta.motorPosition(gos)-(U{rec}.meta.ranges(end)-50000);%x-min/max-min
            gonorm=gotmp/(U{rec}.meta.ranges(end)-(U{rec}.meta.ranges(end)-50000));
        nogotmp=U{rec}.meta.motorPosition(nogos)-(U{rec}.meta.ranges(1));
            nogonorm=(-1+(nogotmp/((U{rec}.meta.ranges(1)+50000)-U{rec}.meta.ranges(1))));
        
        [~,~,~,~,~,touchesPreD] = assist_predecisionVar(U,rec);
        normpos=zeros(length(U{rec}.k),3);
        
        normpos(gos,1)=gonorm;normpos(nogos,1)=nogonorm;
        normpos(:,2)=cellfun(@numel,touchesPreD);
        normpos(:,3)=U{rec}.meta.trialCorrect;
            
        clust=[]; %bin by motor positions
            for k=1:length(normpos)
                clust=binslin(normpos(:,1),normpos(:,1:3),'equalE',11,-1,1);
            end    
         
        range=linspace(-1,1,11);
        
        prob{rec}=cell2mat(cellfun(@mean,clust,'uniformoutput',0));
        prob{rec}(1:5,3)=1-prob{rec}(1:5,3); %since these are nogo trials and we want lick prob
end            


figure(10);clf;plot(1:10,prob{1}(:,3),'k')
hold on;plot(1:10,prob{2}(:,3),'k');
plot(mean([prob{1}(:,3), prob{2}(:,3)],2),'r','linewidth',4);
hold on; bar(1:10,(mean([prob{1}(:,2), prob{2}(:,2)],2)/50)+.25,.5,'FaceColor',[.7 .7 .7]);%mean number of touches per trial
xlabel('Decision Boundary');ylabel('Lick Probability')
title('Expert Session Full Continuous')
set(gca,'xlim',[1 10],'xtick',[1 5.5 10],'xticklabel',{-1 0 1},'ytick',[0:.25:1],'ylim',[.25 1])
a2 = axes('YAxisLocation', 'Right');
set(a2,'color','none');set(a2,'XTick',[]);set(a2,'YLim',[0 2],'YTick',[0:.1:.5],'yticklabel',[0:2:10]),ylabel('# Pre Decision Touches/Trial')
