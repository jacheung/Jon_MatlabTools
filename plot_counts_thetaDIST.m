clear g;clear ng
for rec = 1:length(U)
gos = [V(rec).touchNum.hit V(rec).touchNum.miss];
nogos = [V(rec).touchNum.CR V(rec).touchNum.FA];
totalT=numel([gos nogos]);
goDist = histc(gos,0:10);
nogoDist = histc(nogos,0:10);
figure(50);%subplot(2,5,rec);
bar(0:10,goDist./numel(gos),'b');alpha(.35)
hold on; bar(0:10,nogoDist./numel(nogos),'r');alpha(.35)
set(gca,'xlim',[-.5 11],'xtick',0:5:10,'ylim',[0 .5])

g{rec}=goDist;
ng{rec}=nogoDist;
end

GgoDist = sum(cell2mat(g'));
GngDist = sum(cell2mat(ng'));

figure(51);
bar(0:10,GgoDist./sum(GgoDist),'b');alpha(.35);
hold on; bar(0:10,GngDist./sum(GngDist),'r');alpha(.35);
set(gca,'xlim',[-.5 11],'xtick',0:5:10,'ylim',[0 .5])



%% THETA DIST


gos = [V(rec).var.hit{1} V(rec).var.miss{1}];
nogos = [V(rec).var.CR{1} V(rec).var.FA{1}];
thetas = histc([gos nogos],[-50:50]);
figure(23);subplot(2,1,1)
plot([-50:50],thetas./numel(thetas),'color',[.5 .5 .5]);
set(gca,'xtick',[-50:25:50],'xticklabel',[-50:25:50]);
title('Theta Distribution')
ylabel('Proportion of Touches')
