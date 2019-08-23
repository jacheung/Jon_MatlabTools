

numCell = 2;
gos=find(U{numCell}.meta.trialType==1);
nogos = find(U{numCell}.meta.trialType==0);

figure(480);clf
for i = 1:length(gos)
    figure(480); subplot(4,1,1)
    spks =find( U{numCell}.R_ntk(:,:,gos(i))==1);
    hold on;scatter(spks,ones(length(spks),1)*i,7,'k','filled')
    set(gca,'xtick',[],'xticklabel',[])
end



for k = 1:length(nogos)
    figure(480); subplot(4,1,2)
    spks =find(U{numCell}.R_ntk(:,:,nogos(k))==1);
    hold on;scatter(spks,ones(length(spks),1)*k,7,'k','filled')
end

gospks=mean(squeeze(U{numCell}.R_ntk(:,:,gos)),2);
nogospks=mean(squeeze(U{numCell}.R_ntk(:,:,nogos)),2);

figure(480);subplot(4,1,[3 4]);
hold on;
plot(1:length(gospks),smooth(gospks,50)*1000,'b')
plot(1:length(nogospks),smooth(nogospks,50)*1000,'r');
ylabel('spks/s');xlabel('Time from trial start (ms)');

% ponset = mean(U{numCell}.meta.poleOnset)*10000;
% 
% tmp=squeeze(find(U{numCell}.S_ctk(9,:,:)==1))/U{numCell}.t;
% firsttouch = min(tmp-floor(tmp))*10000;

% hold on; 
% plot([ponset ponset],[0 10],'-.k')
% plot([firsttouch firsttouch],[0 10],'-.k')
