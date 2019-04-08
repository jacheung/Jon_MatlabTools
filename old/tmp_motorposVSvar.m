stdall=[];
for i = 1:length(U)
var=1; %1 for theta    
array=U{i};
go=find(array.meta.trialType==1);
nogo=find(array.meta.trialType==0);
gotheta=cell(1,length(go));pregotheta=cell(1,length(go));
nogotheta=cell(1,length(nogo));prenogotheta=cell(1,length(nogo));

[~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(U,i);

for j = 1:length(go)
    gotouchOnIdx= [];
    gotouchOnIdx = [find(array.S_ctk(9,:,go(j))==1); find(array.S_ctk(12,:,go(j))==1)'];
    gotheta{j}=array.S_ctk(var,gotouchOnIdx,go(j));
    pregotheta{j}=array.S_ctk(var,prelixGo{j},go(j));
end

for k = 1:length(nogo)
    nogotouchOnIdx=[];
    nogotouchOnIdx = [find(array.S_ctk(9,:,nogo(k))==1); find(array.S_ctk(12,:,nogo(k))==1)'];
    nogotheta{k}=array.S_ctk(var,nogotouchOnIdx,nogo(k));
    prenogotheta{k}=array.S_ctk(var,prelixNoGo{k},nogo(k));
end

gomotor=array.meta.motorPosition(go);nogomotor=array.meta.motorPosition(nogo);
gotmean=cellfun(@nanmean, gotheta);nogotmean=(cellfun(@nanmean,nogotheta));
gotstd=cellfun(@nanstd, gotheta);nogotstd=(cellfun(@nanstd,nogotheta));
figure(232);clf
errorbar(gomotor,gotmean,gotstd,'ko'); 
hold on; errorbar(nogomotor,nogotmean,nogotstd,'ko')
xlabel('Motor Position');ylabel('Theta')
stdall(i)=nanmean([gotstd nogotstd]);


% figure(57);clf
% plot(smooth(histc(gotheta,-75:75),5));
% hold on;plot(smooth(histc(pregotheta,-75:75),5),'c');
% hold on;plot(smooth(histc(nogotheta,-75:75),5),'r');
% hold on;plot(smooth(histc(prenogotheta,-75:75),5),'m');

pause

end

