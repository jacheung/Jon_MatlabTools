j=2
touchOnIdx = [find(array.S_ctk(9,:,go(j))==1); find(array.S_ctk(12,:,go(j))==1)'];
figure;plot(1:4000,array.S_ctk(8,:,j));
hold on;scatter(touchOnIdx,zeros(1,length(touchOnIdx)))

%%
rec=34;
varcell=cell(1,2);
[varspikes] = assist_varAtTouch(U{rec},-25:50);
[varForce] = assist_varAtTouchFORCE(U{rec},-25:50);
varcell{1}=varspikes;varcell{2}=varForce;
figure(39); clf


for d=1:numel(varcell)
[sorted sortedBy binBounds]=binslin(varcell{d}(:,1),varcell{d}(:,7:82),'equalE',41,-100,100);
thetarange=[-97.5:5:97.5];

thetatouchspks=zeros(length(sorted),79);
for j=1:length(sorted)
    thetatouchspks(j,2:77)=mean(sorted{j},1);
    thetatouchspks(j,78)=size(sorted{j},1);
end
thetatouchspks(:,1)=thetarange;    
thetatouchspks = thetatouchspks(all(~isnan(thetatouchspks),2),:); %remove NaN rows

for k=1:size(thetatouchspks,1)
    thetatouchspks(k,2:77)=smooth(thetatouchspks(k,2:77));
        if thetatouchspks(k,78)>normbinmin
        thetatouchspks(k,79)=max(thetatouchspks(k,26:26+normbinwin));
        end
end
[~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(U,rec);
[govar, nogovar, ~, ~] = assist_vardistribution(U{rec},1,prelixGo,prelixNoGo);

figure(39);h1 = subplot(numel(varcell),1,d);
imagesc(thetatouchspks(:,2:77))
hold on
colormap(gca,parula)
caxis([0 max(thetatouchspks(:,79))])
plot([25 25],[25 0],'w:') 
set(gca,'Ydir','normal','ytick',(1:length(thetatouchspks(:,1))),'yticklabel',[thetatouchspks(:,1)],'xtick',(0:25:75),'xticklabel',(-25:25:50));
for k=1:size(thetatouchspks,1)
    text(20,k,num2str(thetatouchspks(k,78)),'FontSize',8,'Color','white')
end
axis('square')
%xlabel('time from touch onset (ms)')
title('Theta at touch')

end
figure(39);xlabel('time from touch onset')