type = {D,SM,BV};
for d = 1:length(type) 
U=type{d} ;

clear g;clear ng
[V] = classifierWrapper(U);
for rec = 1:length(U)
    nogotouches = [V(rec).var.FA{1} V(rec).var.CR{1}]; %makes sure there is at least one touch before using array for building distribution
    
    if ~isempty(nogotouches) || numel(nogotouches)>1
        gos = [V(rec).touchNum.hit V(rec).touchNum.miss];
        nogos = [V(rec).touchNum.CR V(rec).touchNum.FA];
        totalT=numel([gos nogos]);
        goDist = histc(gos,0:10);
        nogoDist = histc(nogos,0:10);
        
        figure(50);subplot(2,5,rec);
        bar(0:10,goDist./numel(gos),'b');alpha(.35)
        hold on; bar(0:10,nogoDist./numel(nogos),'r');alpha(.35)
        set(gca,'xlim',[-.5 11],'xtick',0:5:10,'ylim',[0 .5],'ytick',[0 .25 .5])
        
        g{rec}=goDist;
        ng{rec}=nogoDist;
    end
end

GgoDist = sum(cell2mat(g'));
GngDist = sum(cell2mat(ng'));
gotouch = [];
nogotouch = [];
for i = 0:10
    gotouch = [gotouch ; repmat(i,GgoDist(i+1),1)];
    nogotouch = [nogotouch; repmat(i,GngDist(i+1),1)];
end

gos = [mean(gotouch) std(gotouch)]
nogos = [mean(nogotouch) std(nogotouch)]



figure(50+d);
bar(0:10,GgoDist./sum(GgoDist),'b');alpha(.35);
hold on; bar(0:10,GngDist./sum(GngDist),'r');alpha(.35);
set(gca,'xlim',[-.5 11],'xtick',0:5:10,'ylim',[0 .5],'ytick',[0 .25 .5])
% print(figure(51),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_countsDistribution' ])
% print(figure(50),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_countsDistributionINDIV' ])
% 

%% THETA DIST for one single mouse


% gos = [V(rec).var.hit{1} V(rec).var.miss{1}];
% nogos = [V(rec).var.CR{1} V(rec).var.FA{1}];
% thetas = histc([gos nogos],[-50:50]);
% figure(23);subplot(2,1,1)
% plot([-50:50],thetas./numel([gos nogos]),'color',[.5 .5 .5]);
% set(gca,'xtick',[-50:25:50],'xticklabel',[-50:25:50]);
% title('Theta Distribution')
% ylabel('Proportion of Touches')


gos = [];
nogos = [];
groupgorange = zeros(length(U),2);
groupnogorange = zeros(length(U),2);
for var = 1
    if var == 5
        range = linspace(-3.14,3.14,15);
    else
        range = [-25:50];
    end
    for rec = 1:length(U)
        tmpgos = [V(rec).var.hit{var} V(rec).var.miss{var}];
        tmpnogos = [V(rec).var.FA{var} V(rec).var.CR{var}];
        
        tmpgothetas = histc([tmpgos],range);
        tmpnogothetas = histc([tmpnogos],range);
        figure(32); subplot(2,5,rec); 
        plot(range,tmpgothetas,'b');alpha(.5)
        hold on; plot(range,tmpnogothetas,'r');alpha(.5)
        
        thetagovals = range(tmpgothetas>0);
        thetanogovals = range(tmpnogothetas>0);
        groupgorange(rec,1:2) = [min(thetagovals) max(thetagovals)];
        groupnogorange(rec,1:2) = [min(thetanogovals) max(thetanogovals)];
        
        
        if var == 5
            set(gca,'xlim',[min(range) max(range)],'xtick',linspace(-3.14,3.14,3),'xticklabel',{'-\pi','0','\pi'})
        end
        
        gos = [gos V(rec).var.hit{var} V(rec).var.miss{var}];
        nogos = [nogos V(rec).var.CR{var} V(rec).var.FA{var}];
    end
    gothetas = histc([gos],range);
    nogothetas = histc([nogos],range);
    figure(89+d); clf; 
    bar(range,gothetas,'b');alpha(.5)
    hold on;bar(range,nogothetas,'r');alpha(.5)
    if var == 5
        hold on; bar(range,nogothetas,'r');alpha(.5)
        set(gca,'xlim',[min(range) max(range)],'xtick',linspace(-3.14,3.14,5),'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    end
    ylabel('Number of Touches')
%     print(figure(31),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' V(1).varNames{var} 'DistPOP'])
%     print(figure(32),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' V(1).varNames{var} 'DistINDIV'])
% %     close all
end


gowidth = groupgorange(:,2)-groupgorange(:,1);
nogowidth = groupnogorange(:,2)-groupnogorange(:,1);
end