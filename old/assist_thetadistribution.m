function [govar, nogovar, pregovar, prenogovar] = assist_vardistribution(array,variable,prelixGo,prelixNoGo);

go=find(array.meta.trialType==1);
nogo=find(array.meta.trialType==0);
govar=cell(1,length(go));pregovar=cell(1,length(go));
nogovar=cell(1,length(nogo));prenogovar=cell(1,length(nogo));

for j = 1:length(go)
    gotouchOnIdx= [];
    gotouchOnIdx = [find(array.S_ctk(9,:,go(j))==1); find(array.S_ctk(12,:,go(j))==1)'];
    govar{j}=array.S_ctk(variable,gotouchOnIdx,go(j));
    pregovar{j}=array.S_ctk(variable,prelixGo{j},go(j));
end

for k = 1:length(nogo)
    nogotouchOnIdx=[];
    nogotouchOnIdx = [find(array.S_ctk(9,:,nogo(k))==1); find(array.S_ctk(12,:,nogo(k))==1)'];
    nogovar{k}=array.S_ctk(variable,nogotouchOnIdx,nogo(k));
    prenogovar{k}=array.S_ctk(variable,prelixNoGo{k},nogo(k));
end
govar=cell2mat(govar);nogovar=cell2mat(nogovar);
pregovar=cell2mat(pregovar);prenogovar=cell2mat(prenogovar);

figure(57);clf
plot(smooth(histc(govar,-75:75),5));
hold on;plot(smooth(histc(pregovar,-75:75),5),'c');
hold on;plot(smooth(histc(nogovar,-75:75),5),'r');
hold on;plot(smooth(histc(prenogovar,-75:75),5),'m');
set(gca,'xtick',[0:25:150],'xticklabel',[-75:25:75])



