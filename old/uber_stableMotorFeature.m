for rec=1:length(U)
mp=U{rec}.meta.motorPosition;
filt= mp>min(mp)+10000; %farthest nogo filter
mp(filt)=NaN;
farnogoT=find(~isnan(mp)==1);

farnogoonset=U{rec}.meta.poleOnset(farnogoT);


thetas = NaN(length(farnogoT),U{rec}.t);
amp = NaN(length(farnogoT),U{rec}.t);
setpoint = NaN(length(farnogoT),U{rec}.t);
phase = NaN(length(farnogoT),U{rec}.t);

for d = 1:length(farnogoT)
    thetas(d,:)=U{rec}.S_ctk(1,:,farnogoT(d));
    amp(d,:)=U{rec}.S_ctk(3,:,farnogoT(d));
    setpoint(d,:)=U{rec}.S_ctk(4,:,farnogoT(d));
    phase(d,:)=U{rec}.S_ctk(5,:,farnogoT(d));
end

var = {thetas,amp,setpoint,phase};
varNames = {'Theta','Amp','Setpoint','Phase'} ;
varStd = [];
for varNum = 1:length(var) 
figure(20+varNum);clf
subplot(3,2,[1:4])
% plot(1:length(var{varNum}),var{varNum},'color',[.8 .8 .8])
% hold on; plot(1:length(var{varNum}),mean(var{varNum}),'r','linewidth',2)
% hold on; plot([mean(farnogoonset)*1000 mean(farnogoonset)*1000],[min(min(var{varNum}))*1.5 max(max(var{varNum}))*1.5],'-.k') 
tmpVals = [min(min(var{varNum})) max(max(var{varNum}))];
normVar = (var{varNum}-tmpVals(1))./(tmpVals(2)-tmpVals(1));

plot(1:length(normVar),normVar,'color',[.8 .8 .8])
hold on; plot(1:length(normVar),mean(normVar),'r','linewidth',2)
hold on; plot([mean(farnogoonset)*1000 mean(farnogoonset)*1000],[min(min(normVar))*1.5 max(max(normVar))*1.5],'-.k') 
ylabel(varNames{varNum})
set(gca,'xticklabel',[],'xtick',[])
subplot(3,2,[5:6])
% plot(1:length(var{var
Num}),std(var{varNum}));
plot(1:length(var{varNum}),std(normVar));
set(gca,'ylim',[0 .5])
ylabel('Variance')
xlabel('Time from Trial Start (ms)')

varStd(varNum,:) = std(normVar);
end

[~,idx] = min(varStd(2:4,:)); %get index of variable with lowest variance, thus most similar 

figure(26);scatter(1:length(idx),idx);
set(gca,'ylim',[1 3],'ytick',[1 2 3],'yticklabel',varNames(2:4))
title('Motor Feature with Min Variance');
pause
end

