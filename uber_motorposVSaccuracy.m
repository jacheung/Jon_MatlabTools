function [prob,raw] = uber_motorposVSaccuracy(U,varargin)
%varargin = 'semi' or 'cont'. 

%outputs a cell array. Each array = 1 session. 
%Column 1 = normalized pole positions. 2 = AVG number of touches. 3 =
%p(lick)

prob=cell(1,length(U));

for rec=1:length(U)
    gos=find(U{rec}.meta.trialType==1);nogos=find(U{rec}.meta.trialType==0);
    gotmp=U{rec}.meta.motorPosition(gos)-min(U{rec}.meta.motorPosition(gos));%x-min/max-min
    gonorm=gotmp/(max(U{rec}.meta.motorPosition(gos))-min(U{rec}.meta.motorPosition(gos)));
    if strcmp(varargin{1},'semi')
        nogonorm=zeros(1,numel(nogos))-1;
    else
        nogotmp=U{rec}.meta.motorPosition(nogos)-(U{rec}.meta.ranges(1));
        nogonorm=(-1+(nogotmp/((U{rec}.meta.ranges(1)+50000)-U{rec}.meta.ranges(1))));
    end
    
    [touchesPreD,~,~,~,~,~] =  assist_predecisionVar_vOLD(U,rec);

    normpos=zeros(length(U{rec}.k),3);
    
    normpos(gos,1)=gonorm;normpos(nogos,1)=nogonorm;
    normpos(:,2)=cellfun(@numel,touchesPreD);
    normpos(:,3)=U{rec}.meta.trialCorrect;
    
    clust=[]; %bin by motor positions
    for k=1:length(normpos)
        clust=binslin(normpos(:,1),normpos(:,1:3),'equalE',11,-1,1);
    end
    
    raw{rec} = cell2mat(clust);
    
    if strcmp(varargin{1},'semi')
        prob{rec}=cell2mat(cellfun(@(x) mean(x,1),clust([1 5:end]),'uniformoutput',0));
        prob{rec}(1,3)=1-prob{rec}(1,3); %since these are nogo trials and we want lick prob
    else
        prob{rec}=cell2mat(cellfun(@(x) mean(x,1),clust,'uniformoutput',0));
        prob{rec}(1:5,3)=1-prob{rec}(1:5,3);
    end
end



if strcmp(varargin{1},'semi')
    figure(9);clf;%plot(6:11,prob{1}(end-5:end,3),'k')
    for i = 1:length(prob)
        hold on;plot(6:11,prob{i}(end-5:end,3),'color',[.8 .8 .8]);
        plot(0,prob{i}(1,3),'ko')
    end
    tmp=cell2mat(prob);
    allavg = mean(tmp(:,3:3:end),2);
    numavg = mean(tmp(:,2:3:end),2);
    plot(6:11,allavg(end-5:end),'r','linewidth',4);
    plot(0,allavg(1),'ro','linewidth',4);
    hold on; bar(6:11,numavg(end-5:end)/50,1,'FaceColor',[.7 .7 .7]);%mean number of touches per trial
    hold on; bar(0,numavg(1)/50,1,'FaceColor',[.7 .7 .7]);%mean number of touches per trial
    xlabel('Decision Boundary');ylabel('Lick Probability')
    title(['Expert Session Semi Continuous n = ' num2str(length(U))])
    set(gca,'xlim',[0 11],'xtick',[0 5.5 11],'xticklabel',{-1 0 1},'ytick',[0:.25:1],'ylim',[0 1])
    a2 = axes('YAxisLocation', 'Right');
    set(a2,'color','none');set(a2,'XTick',[]);set(a2,'YLim',[0 2],'ytick',[0:.1:.4],'yticklabel',[0:2.5:10]);ylabel('# Pre Decision Touches/Trial')

else
    figure(10);clf;plot(1:10,prob{1}(:,3),'k')
    for i = 1:length(prob)
        hold on;plot(1:10,prob{i}(:,3),'color',[.8 .8 .8]);
    end
    tmp=cell2mat(prob);
    allavg = mean(tmp(:,3:3:end),2);
    numavg = mean(tmp(:,2:3:end),2);
    plot(1:10,allavg,'r','linewidth',4);
    hold on; bar(1:10,numavg/50,1,'FaceColor',[.7 .7 .7]);%mean number of touches per trial
    xlabel('Decision Boundary');ylabel('Lick Probability')
    title(['Expert Session Full Continuous n = ' num2str(length(U))])
    set(gca,'xlim',[1 10],'xtick',[1 5.5 10],'xticklabel',{-1 0 1},'ytick',[0:.25:1],'ylim',[0 1])
    a2 = axes('YAxisLocation', 'Right');
    set(a2,'color','none');set(a2,'XTick',[]);set(a2,'YLim',[0 2],'ytick',[0:.1:.4],'yticklabel',[0:2.5:10]);ylabel('# Pre Decision Touches/Trial')

end

